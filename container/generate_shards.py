#!/usr/bin/env python3

"""Generate shards for Sentieon"""

from __future__ import annotations

import argparse
import copy
import io
import json
import os
import pathlib
import re
import sys

from collections import OrderedDict
from typing import Generator, Optional

interval_pat = re.compile(r"(?P<chrom>.*):(?P<start>\d+)-(?P<stop>\d+)")


class Faidx:
    """Manipulate genomic regions from a samtools faidx index"""

    def __init__(self, fh: io.TextIOWrapper):
        contigs: OrderedDict[str, int] = OrderedDict()
        for line in fh:
            line = line.rstrip().split("\t")
            chrom, length = line[:2]
            contigs[chrom] = int(length)
        self.contigs = contigs


class IntervalList(object):
    """A list of intervals"""

    def __init__(
        self,
        path: Optional[pathlib.Path] = None,
        faidx: Optional[Faidx] = None,
        interval_str: Optional[str] = None,
    ):
        if path:
            regions = self.load_bed(path)
        elif faidx and interval_str:
            regions = self.load_interval_str(faidx, interval_str)
        elif faidx:
            regions = self.load_faidx(faidx)
        else:
            regions = {}
        self.regions = self.merge_flatten(regions)

    def __str__(self) -> str:
        res: list[str] = []
        for chrom, cur_intervals in self.regions.items():
            for i in range(0, len(cur_intervals), 2):
                s, e = cur_intervals[i : i + 2]
                res.append(f"{chrom}:{s}-{e}")
        return ",".join(res)

    @staticmethod
    def load_faidx(faidx: Faidx) -> dict[str, list[tuple[int, int]]]:
        """Load all of the regions in a faidx"""
        regions: dict[str, list[tuple[int, int]]] = {}
        for chrom, stop in faidx.contigs.items():
            regions[chrom] = [(1, stop)]
        return regions

    @staticmethod
    def load_interval_str(
        faidx: Faidx, interval_str: str
    ) -> dict[str, list[tuple[int, int]]]:
        """Load a comma-delimited list of intervals"""
        regions: dict[str, list[tuple[int, int]]] = {}
        for interval in interval_str.split(","):
            m = interval_pat.match(interval)
            if m:
                d = m.groupdict()
                chrom, start, stop = d["chrom"], d["start"], d["stop"]
                regions.setdefault(chrom, []).append((int(start), int(stop)))
            else:
                stop = faidx.contigs[interval]
                regions.setdefault(interval, []).append((1, stop))
        return regions

    @staticmethod
    def merge_flatten(
        regions: dict[str, list[tuple[int, int]]],
        dist: int = 1,
    ) -> dict[str, list[int]]:
        """Merge overlapping regions"""
        out_regions: dict[str, list[int]] = {}
        for chrom, intervals in regions.items():
            v: list[int] = []
            s0 = e0 = None
            for s, e in sorted(intervals):
                if e0 is None or s0 is None:
                    s0 = s
                    e0 = e
                elif s - dist > e0:
                    v.extend((s0, e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = e
            if e0 is not None and s0 is not None:
                v.extend((s0, e0))
            out_regions[chrom] = v
        return out_regions

    @staticmethod
    def load_bed(path: pathlib.Path) -> dict[str, list[tuple[int, int]]]:
        """Load a bed file"""
        regions: dict[str, list[tuple[int, int]]] = {}
        fp = io.open(path, "rb")
        for line in fp:
            if line.startswith(b"#") or line.startswith(b"track"):
                continue
            cols = line.rstrip().decode().split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[:3]
            regions.setdefault(chrom, []).append((int(start), int(end)))
        fp.close()
        return regions

    def split_into_n(self, n: int) -> list[IntervalList]:
        """Split the interval into n equal parts"""
        size = 0
        for se_list in self.regions.values():
            for i in range(0, len(se_list), 2):
                s, e = se_list[i : i + 2]
                size += e - s
        split_size = size // n
        carry = size % n

        return self.split_by_bases(split_size, carry)

    def split_by_bases(self, n_bases: int, carry: int = 0) -> list[IntervalList]:
        """Split an interval into parts of n_bases each"""

        res: list[IntervalList] = []
        chroms = list(self.regions.keys())
        chrom_idx, intv_idx, offset = 0, 0, 0
        while chrom_idx < len(chroms):
            cur_size = 0
            regions: dict[str, list[int]] = {}
            target_size = (n_bases + 1) if carry else n_bases
            carry = max(0, carry - 1)

            while True:
                if chrom_idx >= len(chroms):
                    break

                cur_chrom = chroms[chrom_idx]
                cur_intervals = self.regions[cur_chrom]

                if intv_idx + 1 >= len(cur_intervals):
                    chrom_idx += 1
                    intv_idx = 0
                    continue

                s, e = cur_intervals[intv_idx : intv_idx + 2]
                size = e - s - offset
                if size + cur_size >= target_size:
                    diff = size + cur_size - target_size
                    regions.setdefault(cur_chrom, []).extend((s + offset, e - diff))
                    offset = (e - diff) - s + 1
                    break

                regions.setdefault(cur_chrom, []).extend((s + offset, e))
                cur_size += size
                offset = 0
                intv_idx += 2

            if regions:
                _res = IntervalList()
                _res.regions = regions
                res.append(_res)

        return res

    def head(self, n_bases: int) -> IntervalList:
        """Create a new IntervalList with only the first `n_bases`"""
        cur_size = 0
        out_regions: dict[str, list[int]] = {}
        for chrom, regions in self.regions.items():
            for i in range(0, len(regions), 2):
                s, e = regions[i : i + 2]
                size = e - s
                if size + cur_size >= n_bases:
                    diff = size + cur_size - n_bases
                    out_regions.setdefault(chrom, []).extend((s, e - diff))
                    cur_size += size
                    break

                out_regions.setdefault(chrom, []).extend((s, e))
                cur_size += size

            if cur_size >= n_bases:
                break

        ret = IntervalList()
        ret.regions = out_regions
        return ret

    def after(self, n_bases: int) -> IntervalList:
        """Create a new IntervalList without the first `n_bases`"""
        regions: dict[str, list[int]] = {}
        found = False
        cur_size = 0
        chrom_idx, intv_idx = 0, 0
        chroms = list(self.regions.keys())
        chrom = chroms[chrom_idx]
        intv_list = self.regions[chrom]

        while True:
            while intv_idx < len(intv_list):
                s, e = intv_list[intv_idx : intv_idx + 2]
                size = e - s
                intv_idx += 2
                if size + cur_size >= n_bases:
                    diff = size + cur_size - n_bases
                    if diff > 0:
                        regions.setdefault(chrom, []).extend((e - diff, e))
                    found = True
                    break
                cur_size += size

            if found:
                break

            chrom_idx += 1
            intv_idx = 0
            if chrom_idx >= len(chroms):
                break

            chrom = chroms[chrom_idx]
            intv_list = self.regions[chrom]

        if found:
            while intv_idx < len(intv_list):
                s, e = intv_list[intv_idx : intv_idx + 2]
                regions.setdefault(chrom, []).extend((s, e))
                intv_idx += 2

            chrom_idx += 1
            while chrom_idx < len(chroms):
                chrom = chroms[chrom_idx]
                regions[chrom] = copy.copy(self.regions[chrom])
                chrom_idx += 1

        ret = IntervalList()
        ret.regions = regions
        return ret


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fai",
        required=True,
        type=argparse.FileType("r"),
        help="The reference fasta index file",
    )
    parser.add_argument(
        "--region",
        type=str,
        help="A comma-delimited list of regions to limit the joint call to",
    )
    parser.add_argument(
        "--n_workers",
        required=True,
        type=int,
        help="The number of GVCFtyper workers",
    )
    parser.add_argument(
        "--shard_size",
        required=True,
        type=int,
        help="The shard size",
    )
    parser.add_argument(
        "--n_merge_splits",
        required=True,
        type=int,
        help="The number of splits for the output VCF",
    )
    parser.add_argument(
        "--shards_json",
        required=True,
        type=argparse.FileType("w"),
        help="The output file containing a list of lists of shard coordinates",
    )
    parser.add_argument(
        "--merge_info",
        required=True,
        type=argparse.FileType("w"),
        help="The output file containing the merge information",
    )
    return parser.parse_args(argv)


def split_int(size: int, parts: int) -> Generator[tuple[int, int], None, None]:
    """
    Split `size` into `parts` of aprox equal size
    yield: (offset, split_size)
    """
    split_size = size // parts
    carry = size % parts
    idx = 0
    while idx < size:
        _split_size = (split_size + 1) if carry else split_size
        yield (idx, _split_size)
        idx += _split_size
        carry = max(0, carry - 1)


def main(argv: argparse.Namespace) -> int:
    """Main function"""
    faidx = Faidx(argv.fai)
    intervals = IntervalList(faidx=faidx, interval_str=argv.region)

    # Handle worker intervals
    shards = intervals.split_by_bases(argv.shard_size)
    out_shards = [
        shards[idx : idx + split_size]
        for idx, split_size in split_int(len(shards), argv.n_workers)
    ]

    res: list[list[str]] = []
    for worker in out_shards:
        res.append([str(x) for x in worker])
    json.dump(res, argv.shards_json)

    # Handle merge intervals
    merge_info = {}
    shards_per_merge = int(len(shards) / argv.n_merge_splits + 0.999)
    merge_info["size"] = shards_per_merge + 1

    out_merge_shards: list[tuple[int, str]] = []
    for i in range(argv.n_merge_splits):
        # Merge the covered shards + 1kb into the next one
        beg = i * shards_per_merge
        end = i * shards_per_merge + shards_per_merge
        merge_shards = shards[beg:end]
        if not merge_shards:
            continue

        if i > 0:  # Not the first merge
            merge_shards[0] = merge_shards[0].after(1001)

        if end < len(shards):  # Not the last merge
            merge_shards.append(shards[end].head(1000))
        merge_shards = IntervalList(
            faidx=faidx,
            interval_str=",".join([str(x) for x in merge_shards]),
        )
        out_merge_shards.append((beg, str(merge_shards)))
    merge_info["shards"] = out_merge_shards
    json.dump(merge_info, argv.merge_info)

    return os.EX_OK


if __name__ == "__main__":
    args = parse_args()
    sys.exit(main(args))
