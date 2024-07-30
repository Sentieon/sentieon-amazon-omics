#!/usr/bin/env python3

"""Merge sharded VCF files"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile

from typing import Optional


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vcfs",
        required=True,
        nargs="+",
        help="The sharded VCFs, in order",
    )
    parser.add_argument(
        "--region",
        help="Output the VCF over a specific region",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The output VCF file",
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


def main(argv: argparse.Namespace) -> int:
    """Main function"""

    vcf_list = argv.vcfs
    vcf_str = " ".join(vcf_list)
    output = argv.output

    tmpdir = tempfile.mkdtemp()
    vcf_fifo = os.path.join(tmpdir, "merged_tmp.vcf.gz")
    os.mkfifo(vcf_fifo)

    region_tmp_file = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".bed",
        delete=False,
    )
    region_tmp_file_name = region_tmp_file.name
    for region in argv.region.split(','):
        contig, coords = region.rsplit(':', 1)
        start, stop = coords.split('-', 1)
        print("\t".join((contig, start, stop)), file=region_tmp_file)
    region_tmp_file.close()

    subset_cmd = (
        f"cat {vcf_fifo} "
        f"| bcftools view -T {region_tmp_file_name} "
        f"| sentieon util vcfconvert - {output}"
    )
    logging.debug("Running: %s", subset_cmd)
    subset_p = subprocess.Popen(subset_cmd, shell=True)

    merge_cmd = (
        f"sentieon driver --passthru --algo GVCFtyper --merge {vcf_fifo} {vcf_str}"
    )
    logging.debug("Running: %s", merge_cmd)
    merge_p = subprocess.Popen(merge_cmd, shell=True)

    merge_p.wait()
    subset_p.wait()

    shutil.rmtree(tmpdir)

    if merge_p.returncode != 0 or subset_p.returncode != 0:
        logging.error("Merge failed")
        return -1
    return 0


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(level=args.loglevel)
    sys.exit(main(args))
