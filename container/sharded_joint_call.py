#!/usr/bin/env python3

"""Run joint calling using a waterfall approach"""

from __future__ import annotations

import argparse
import asyncio
import asyncio.subprocess
import logging
import os
import pathlib
import shutil
import sys

from typing import Any, Optional

DOWNLOAD_CMD = (
    "bcftools view --no-version -r {shard} --threads 2 -o {gvcf_dest} {gvcf} && "
    "bcftools index --threads 2 -t {gvcf_dest}"
)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ref",
        required=True,
        type=argparse.FileType("r"),
        help="The reference fasta",
    )
    parser.add_argument(
        "--gvcf_list",
        required=True,
        type=argparse.FileType("r"),
        help="The list of gVCF input files",
    )
    parser.add_argument(
        "--shards",
        required=True,
        nargs="+",
        help="Shards to process",
    )
    parser.add_argument(
        "--dbsnp",
        type=argparse.FileType("r"),
        help="A dbsnp VCF",
    )
    parser.add_argument(
        "--driver_xargs",
        help="Extra arguments for the sentieon driver",
    )
    parser.add_argument(
        "--gvcftyper_xargs",
        help="Extra arguments for GVCFtyper",
    )
    parser.add_argument(
        "--output_basename",
        required=True,
        help="The basename for the output VCFs",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        default=logging.WARNING,
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Print debugging info",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
    )
    parser.add_argument(
        "--concurrent_downloads",
        help="The number of gVCFs to ingest concurrently",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--download_retries",
        help="The number of times to retry a failing download",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--ramdisk",
        help="Use a ramdisk to store partial gVCFs",
        action="store_true",
    )
    return parser.parse_args(argv)


async def download_gvcf(
    shard: str,
    gvcf: str,
    gvcf_dest: str,
    current_try: int = 0,
) -> tuple[asyncio.subprocess.Process, str, str, int]:
    """Download one gVCF"""
    if current_try > 0:
        logging.debug(
            "Re-running download for shard %s with gvcf %s",
            shard,
            gvcf,
        )
        await asyncio.sleep(current_try * 2)
    download_cmd = DOWNLOAD_CMD.format(
        shard=shard,
        gvcf=gvcf,
        gvcf_dest=gvcf_dest,
    )
    logging.debug("Running: %s", download_cmd)
    if logging.root.level <= logging.DEBUG:
        proc = await asyncio.create_subprocess_shell(download_cmd, shell=True)
    else:
        proc = await asyncio.create_subprocess_shell(
            download_cmd,
            stderr=asyncio.subprocess.DEVNULL,
            stdout=asyncio.subprocess.DEVNULL,
        )
    res = (proc, gvcf, gvcf_dest, current_try + 1)
    return res


async def download_shard(
    gvcf_list: list[str],
    shard: str,
    shard_idx: int,
    n_concurrent: int,
    download_retries: int = 5,
    ramdisk: bool = False,
) -> tuple[int, int, str]:
    """
    Async download of the next shard
    return: (returncode, shard_idx, shard)
    """
    logging.info("Downloading shard index '%s' and shard: %s", shard_idx, shard)

    base_dir = f"sharded_inputs_{shard_idx}"
    if ramdisk:
        base_dir = f"/dev/shm/sharded_inputs_{shard_idx}"
    os.makedirs(base_dir, exist_ok=True)
    running_downloads: list[tuple[asyncio.subprocess.Process, str, str, int]] = []
    running_tasks: set[asyncio.Task[int]] = set()
    finalized_downloads: list[tuple[str, str]] = []
    cur_gvcf = 0
    while True:
        logging.debug("Download looping shard %s with gvcf index %s", shard, cur_gvcf)
        if len(running_downloads) >= n_concurrent:
            _done, running_tasks = await asyncio.wait(
                running_tasks, return_when=asyncio.FIRST_COMPLETED
            )

            finished_idxs = [
                i
                for i, j in enumerate(running_downloads)
                if j[0].returncode is not None
            ]
            finished_downloads = [
                running_downloads.pop(i) for i in reversed(finished_idxs)
            ]
            for download in finished_downloads:
                proc, gvcf, gvcf_dest, current_try = download
                if proc.returncode == 0:
                    finalized_downloads.append((gvcf, gvcf_dest))
                    continue

                if current_try >= download_retries:
                    logging.error(
                        "Download failed for shard, %s with gvcf: %s",
                        shard,
                        gvcf,
                    )
                    return (-1, shard_idx, shard)

                download_ret = await download_gvcf(
                    shard=shard,
                    gvcf=gvcf,
                    gvcf_dest=gvcf_dest,
                    current_try=current_try,
                )
                running_tasks.add(asyncio.create_task(download_ret[0].wait()))
                running_downloads.append(download_ret)
            continue

        if cur_gvcf >= len(gvcf_list):
            break
        gvcf = gvcf_list[cur_gvcf]
        gvcf_dest = f"{base_dir}/sample_{cur_gvcf}.g.vcf.gz"

        download_ret = await download_gvcf(
            shard=shard,
            gvcf=gvcf,
            gvcf_dest=gvcf_dest,
        )
        running_tasks.add(asyncio.create_task(download_ret[0].wait()))
        running_downloads.append(download_ret)
        cur_gvcf += 1

    while running_downloads:
        if running_tasks:
            _done, running_tasks = await asyncio.wait(
                running_tasks, return_when=asyncio.FIRST_COMPLETED
            )

        finished_idxs = [
            i for i, j in enumerate(running_downloads) if j[0].returncode is not None
        ]
        finished_downloads = [running_downloads.pop(i) for i in reversed(finished_idxs)]

        for download in finished_downloads:
            proc, gvcf, gvcf_dest, current_try = download
            if proc.returncode == 0:
                finalized_downloads.append((gvcf, gvcf_dest))
                continue

            if current_try >= download_retries:
                logging.error(
                    "Download failed for shard, %s with gvcf: %s",
                    shard,
                    gvcf,
                )
                return (-1, shard_idx, shard)

            download_ret = await download_gvcf(
                shard=shard,
                gvcf=gvcf,
                gvcf_dest=gvcf_dest,
                current_try=current_try,
            )
            running_tasks.add(asyncio.create_task(download_ret[0].wait()))
            running_downloads.append(download_ret)

    # Check the dest files
    for gvcf, gvcf_dest in finalized_downloads:
        gvcf_path = pathlib.Path(gvcf_dest)
        if gvcf_path.exists() and (gvcf_path.stat().st_size > 0):
            continue

        logging.warning(
            "Empty destination file for shard, %s with gvcf: %s",
            shard,
            gvcf,
        )
        download_ret = await download_gvcf(shard, gvcf, gvcf_dest)
        proc = download_ret[0]
        ret = await proc.wait()
        if ret != 0 or (not gvcf_path.exists()) or (gvcf_path.stat().st_size < 1):
            logging.error(
                "Final download failed for shard, %s with gvcf: %s",
                shard,
                gvcf,
            )
            return (-1, shard_idx, shard)

    logging.info(
        "Download finished for shard index '%s' and shard: %s", shard_idx, shard
    )
    return (0, shard_idx, shard)


async def run_shard(
    basename: str,
    ref: str,
    n_samples: int,
    shard: str,
    shard_idx: int,
    driver_xargs: str = "",
    algo_xargs: str = "",
    dbsnp: Optional[str] = None,
    ramdisk: bool = False,
) -> int:
    """Async run the GVCFtyper for the shard"""
    base_dir = ""
    if ramdisk:
        base_dir = "/dev/shm/"
    input_gvcfs = [
        f"{base_dir}sharded_inputs_{shard_idx}/sample_{i}.g.vcf.gz"
        for i in range(n_samples)
    ]
    input_gvcfs = str.encode("\n".join(input_gvcfs))

    dbsnp_arg = "--dbsnp {dbsnp}" if dbsnp else ""
    run_cmd = (
        f"sentieon driver -r {ref} --shard {shard} {driver_xargs} "
        f"--algo GVCFtyper {algo_xargs} {dbsnp_arg} "
        f"{basename}_shard-{shard_idx}.vcf.gz -"
    )
    logging.info("Running: %s", run_cmd)
    p = await asyncio.create_subprocess_shell(
        run_cmd,
        stdin=asyncio.subprocess.PIPE,
    )
    await p.communicate(input=input_gvcfs)
    ret = await p.wait()
    if ret != 0:
        logging.error("Joint calling failed for shard %s", shard)
        return -1

    # Remove the input gVCFs to save space
    shutil.rmtree(f"{base_dir}sharded_inputs_{shard_idx}")

    return 0


async def main(argv: argparse.Namespace) -> int:
    """Main function"""

    # Read the list of gVCFs
    gvcf_list = argv.gvcf_list.read().rstrip().split("\n")

    shards_to_process: list[str] = argv.shards
    downloaded_shards: list[tuple[int, str]] = []

    running_downloads: list[asyncio.Task[tuple[int, int, str]]] = []
    running_shards: list[asyncio.Task[int]] = []
    shard_idx = 0
    while shards_to_process or downloaded_shards or running_downloads or running_shards:
        logging.debug(
            (
                "Main loop: iteration with\nshards_to_process: %s\n"
                "downloaded_shards: %s\nrunning_downloads: %s\n"
                "running_shards: %s\nshard_idx: %s\n\n"
            ),
            shards_to_process,
            downloaded_shards,
            running_downloads,
            running_shards,
            shard_idx,
        )

        # Check if the download is finished
        finished: list[int] = []
        for i, task in enumerate(running_downloads):
            if task.done():
                finished.append(i)
                res = task.result()
                if res[0] != 0:
                    logging.error("Exit: download failed for shard %i", res[1])
                    return -1
                downloaded_shards.append(res[1:])
        for i in reversed(finished):
            running_downloads.pop(i)

        if shards_to_process and not running_downloads and not downloaded_shards:
            # Download the next shard to process
            shard_idx += 1
            cur_shard = shards_to_process.pop(0)
            running_downloads.append(
                asyncio.create_task(
                    download_shard(
                        gvcf_list,
                        cur_shard,
                        shard_idx,
                        argv.concurrent_downloads,
                        download_retries=argv.download_retries,
                        ramdisk=argv.ramdisk,
                    )
                )
            )

        # Check if the currently running shard is finished
        finished: list[int] = []
        for i, task in enumerate(running_shards):
            if task.done():
                finished.append(i)
                res = task.result()
                if res != 0:
                    logging.error("Exit after failed run")
                    return -1
        for i in reversed(finished):
            running_shards.pop(i)

        if not running_shards and downloaded_shards:
            # Run the next shard
            shard_idx, shard = downloaded_shards.pop()
            dbsnp_arg = None
            if hasattr(argv, "dbsnp") and hasattr(argv.dbsnp, "name"):
                dbsnp_arg = argv.dbsnp.name
            running_shards.append(
                asyncio.create_task(
                    run_shard(
                        argv.output_basename,
                        argv.ref.name,
                        len(gvcf_list),
                        shard,
                        shard_idx,
                        driver_xargs=argv.driver_xargs,
                        algo_xargs=argv.gvcftyper_xargs,
                        dbsnp=dbsnp_arg,
                        ramdisk=argv.ramdisk,
                    )
                )
            )

        logging.debug("Main loop: checking to wait")
        if (
            (running_shards and running_downloads)
            or (running_shards and downloaded_shards)
            or (running_downloads and not downloaded_shards)
            or (running_shards and not shards_to_process)
        ):
            logging.debug("Main loop: wait")
            running_tasks: list[asyncio.Task[Any]] = [
                *running_downloads,
                *running_shards,
            ]
            await asyncio.wait(running_tasks, return_when=asyncio.FIRST_COMPLETED)
            logging.debug("Main loop: wake")

    return os.EX_OK


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=args.loglevel,
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    sys.exit(asyncio.run(main(args)))
