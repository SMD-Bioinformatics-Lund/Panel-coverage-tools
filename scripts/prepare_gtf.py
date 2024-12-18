#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Set


def main(in_gtf: Path, panel_json: Path, out_gtf: Path):

    gtf_entries = []

    with in_gtf.open() as in_fh:
        for row in in_fh:
            if row.startswith("#"):
                continue



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_gtf", required=True)
    parser.add_argument("--panel_genes", required=True)
    parser.add_argument("--out_gtf", required=True)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(Path(args.in_gtf), Path(args.panel_genes), Path(args.out_gtf))
