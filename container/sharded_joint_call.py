#!/usr/bin/env python3

"""Run joint calling using a waterfall approach"""

from __future__ import annotations

import argparse
import asyncio
import asyncio.subprocess
import logging
import os
import shutil
import sys

from typing import Any, Optional


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
    return parser.parse_args(argv)


async def download_shard(
    gvcf_list: list[str],
    shard: str,
    shard_idx: int,
    n_concurrent: int,
) -> tuple[int, int, str]:
    """
    Async download of the next shard
    return: (returncode, shard_idx, shard)
    """
    logging.info("Downloading shard index '%s' and shard: %s", shard_idx, shard)
    os.makedirs(f"sharded_inputs_{shard_idx}", exist_ok=True)
    running_downloads: list[tuple[asyncio.subprocess.Process, str]] = []
    running_tasks: set[asyncio.Task[int]] = set()
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
                if download[0].returncode != 0:
                    logging.error(
                        "Download failed for shard, %s with gvcf: %s",
                        shard,
                        download[1],
                    )
                    return (-1, shard_idx, shard)

        if cur_gvcf >= len(gvcf_list):
            break
        gvcf = gvcf_list[cur_gvcf]
        download_cmd = (
            f"bcftools view --no-version -r {shard} -o - {gvcf} | "
            f"sentieon util vcfconvert - sharded_inputs_{shard_idx}/sample_{cur_gvcf}.g.vcf.gz"
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
        running_tasks.add(asyncio.create_task(proc.wait()))
        running_downloads.append((proc, gvcf))
        cur_gvcf += 1

    logging.debug("Download waiting for shard %s", shard)
    await asyncio.wait(running_tasks, return_when=asyncio.ALL_COMPLETED)
    for download in running_downloads:
        if download[0].returncode != 0:
            logging.error(
                "Download failed for shard, %s with gvcf: %s", shard, download[1]
            )
            return (-1, shard_idx, shard)
    logging.info("Download finished for shard index '%s' and shard: %s", shard_idx, shard)

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
) -> int:
    """Async run the GVCFtyper for the shard"""
    input_gvcfs = [
        f"sharded_inputs_{shard_idx}/sample_{i}.g.vcf.gz" for i in range(n_samples)
    ]
    input_gvcfs = str.encode("\n".join(input_gvcfs))

    dbsnp_arg = "--dbsnp {dbsnp}" if dbsnp else ""
    run_cmd = (
        f"sentieon driver -r {ref} --shard {shard} {driver_xargs} "
        f"--algo GVCFtyper {algo_xargs} {dbsnp_arg} "
        f"{basename}_shard-{shard_idx}.vcf.gz -"
    )
    logging.info("Running: %s", run_cmd)
    p = await asyncio.subprocess.create_subprocess_shell(
        run_cmd,
        stdin=asyncio.subprocess.PIPE,
    )
    await p.communicate(input=input_gvcfs)
    ret = await p.wait()
    if ret != 0:
        logging.error("Joint calling failed for shard %s", shard)
        return -1

    # Remove the input gVCFs to save space
    shutil.rmtree(f"sharded_inputs_{shard_idx}")

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
                        gvcf_list, cur_shard, shard_idx, argv.concurrent_downloads
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
    logging.basicConfig(level=args.loglevel)
    sys.exit(asyncio.run(main(args)))
