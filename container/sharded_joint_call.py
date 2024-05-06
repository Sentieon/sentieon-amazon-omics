#!/usr/bin/env python3

"""Run joint calling using a waterfall approach"""

from __future__ import annotations

import argparse
import asyncio
import asyncio.subprocess
import logging
import os
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
    return parser.parse_args(argv)


async def run_shard(
    basename: str,
    ref: str,
    gvcf_list: str,
    shard: str,
    shard_idx: int,
    driver_xargs: str = "",
    algo_xargs: str = "",
    dbsnp: Optional[str] = None,
) -> int:
    """Async run the GVCFtyper for the shard"""

    dbsnp_arg = "--dbsnp {dbsnp}" if dbsnp else ""

    with open(gvcf_list, mode="rb") as gvcf_list_fh:
        run_cmd = (
            f"sentieon driver -r {ref} --shard {shard} {driver_xargs} "
            f"--algo GVCFtyper {algo_xargs} {dbsnp_arg} "
            f"{basename}_shard-{shard_idx}.vcf.gz -"
        )
        logging.info("Running: %s", run_cmd)
        p = await asyncio.create_subprocess_shell(
            run_cmd,
            stdin=gvcf_list_fh,
        )
        ret = await p.wait()
        if ret != 0:
            logging.error("Joint calling failed for shard %s", shard)
            return -1
    return 0


async def main(argv: argparse.Namespace) -> int:
    """Main function"""

    shards_to_process: list[str] = argv.shards
    running_shards: list[asyncio.Task[int]] = []
    shard_idx = 0
    while shards_to_process or running_shards:
        logging.debug(
            (
                "Main loop: iteration with\nshards_to_process: %s\n"
                "running_shards: %s\nshard_idx: %s\n\n"
            ),
            shards_to_process,
            running_shards,
            shard_idx,
        )

        if shards_to_process and not running_shards:
            # Run the next shard
            shard_idx += 1
            shard = shards_to_process.pop(0)
            dbsnp_arg = None
            if hasattr(argv, "dbsnp") and hasattr(argv.dbsnp, "name"):
                dbsnp_arg = argv.dbsnp.name
            running_shards.append(
                asyncio.create_task(
                    run_shard(
                        argv.output_basename,
                        argv.ref.name,
                        argv.gvcf_list,
                        shard,
                        shard_idx,
                        driver_xargs=argv.driver_xargs,
                        algo_xargs=argv.gvcftyper_xargs,
                        dbsnp=dbsnp_arg,
                    )
                )
            )

        logging.debug("Main loop: checking to wait")
        if running_shards:
            logging.debug("Main loop: wait")
            running_tasks: list[asyncio.Task[Any]] = [
                *running_shards,
            ]
            await asyncio.wait(
                running_tasks, return_when=asyncio.FIRST_COMPLETED
            )
            logging.debug("Main loop: wake")

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

    return os.EX_OK


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=args.loglevel,
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    sys.exit(asyncio.run(main(args)))
