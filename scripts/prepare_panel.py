#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from typing import Set


def main(panel_json: Path, out_tsv: Path):
    hgnc_symbols: Set[str] = set()
    with panel_json.open() as json_fh:
        json_data = json.load(json_fh)
        hgnc_symbols = set([gene_info["symbol"] for gene_info in json_data["genes"]])

    out_tsv.write_text("\n".join(list(hgnc_symbols)))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel_json", required=True)
    parser.add_argument("--out_tsv", required=True)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(Path(args.panel_json), Path(args.out_tsv))
