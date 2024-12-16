#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

from src.parse_utils import parse_mim2gene, parse_panel_json, parse_panel_text

__version_info__ = ("1", "0", "1")
__version__ = ".".join(__version_info__)

description = "Placeholder description"


def main(panel_genes: List[str]):
    pass

def parse_arguments():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--panel_json", help="Either --panel_json or --panel_txt needs to be supplied")
    parser.add_argument("--panel_text", help="Either --panel_json or --panel_txt needs to be supplied")
    parser.add_argument("--mim2gene", help="OMIM map of gene names to ENSEMBL gene IDs", required=True)
    parser.add_argument("--mane_gtf", help="GTF file with MANE transcripts", required=True)
    parser.add_argument("--d4format", help="Path to d4tools executable", required=True)

    args = parser.parse_args()

    if args.panel_json is None and args.panel_text is None:
        raise ValueError("Either --panel_json or --panel_text needs to be supplied")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    if args.panel_text is None and args.panel_json is not None:
        panel_genes = parse_panel_json(Path(args.panel_json))
    else:
        panel_genes = parse_panel_text(Path(args.panel_text))

    hgnc_to_ensembl = parse_mim2gene(args.mim2gene)

    main(
        panel_genes,
    )