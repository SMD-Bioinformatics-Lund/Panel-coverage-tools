#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

from src.classes import GTFEntry
from src.parse_utils import parse_gtf


def main(in_gtf: Path, panel_names: Path, out_gtf: Path, keep_chr: bool):

    panel_hgnc_names = panel_names.read_text().split("\n")
    gtf_entries = parse_gtf(in_gtf, keep_chr, with_progress=True)

    print(f"{len(gtf_entries)} unfiltered entries found")

    filtered_entries: List[GTFEntry] = []
    for entry in gtf_entries:
        if entry.hgnc_name in panel_hgnc_names:
            if entry.mol_type == "gene" or entry.has_mane_tag():
                filtered_entries.append(entry)

    print(f"Writing {len(filtered_entries)} filtered rows to {out_gtf}")

    with out_gtf.open("w") as out_fh:
        for entry in filtered_entries:
            print(entry.raw_line, file=out_fh)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_gtf", required=True)
    parser.add_argument("--panel_genes", required=True)
    parser.add_argument("--out_gtf", required=True)
    parser.add_argument(
        "--keep_chr", action="store_true", help="Assign if you want to keep the chr prefix"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(Path(args.in_gtf), Path(args.panel_genes), Path(args.out_gtf), args.keep_chr)
