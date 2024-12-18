#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List, Dict, Set

from src.output import write_output
from src.exons_utils import calculate_exon_coverage
from src.d4tools_utils import (
    Coverage,
    collect_d4_coverages,
)
from src.parse_utils import Gene, parse_mane_gtf, parse_panel_text

__version_info__ = ("1", "0", "1")
__version__ = ".".join(__version_info__)

description = "Calculate coverage for gene, MANE transcript and MANE transcript exons"


def main(
    d4tools_command: str,
    d4_file: Path,
    panel_hgnc_symbols: Set[str],
    all_genes: List[Gene],
    out_dir: Path,
    thresholds: List[int],
):

    panel_genes = [gene for gene in all_genes if gene.hgnc_symbol in panel_hgnc_symbols]
    found_symbols = [gene.hgnc_symbol for gene in panel_genes]

    print(f"Now I have {len(panel_genes)} omim genes")

    missed_hgnc_symbols = [
        hgnc_symbol for hgnc_symbol in panel_hgnc_symbols if hgnc_symbol not in found_symbols
    ]

    print(f"Missed {len(missed_hgnc_symbols)} HGNC symbols: {missed_hgnc_symbols}")
    out_dir.mkdir(parents=True, exist_ok=True)

    d4tools_out_dir = out_dir / "d4tools_results"
    d4tools_out_dir.mkdir(parents=True, exist_ok=True)

    gene_bed_rows: List[str] = []
    mane_bed_rows: List[str] = []
    for gene in panel_genes:
        gene_bed = gene.get_bed_row()
        gene_bed_rows.append(gene_bed)
        if gene.mane_transcript:
            mane_bed = gene.get_mane_transcript_bed_row()
            mane_bed_rows.append(mane_bed)

    gene_coverages = collect_d4_coverages(
        gene_bed_rows, out_dir, "gene", d4tools_command, d4_file, thresholds
    )

    mane_coverages = collect_d4_coverages(
        mane_bed_rows, out_dir, "mane", d4tools_command, d4_file, thresholds
    )

    all_exon_bed_rows: List[str] = []
    for gene in panel_genes:
        if gene.mane_transcript:
            for bed_row in gene.get_bed_exons():
                all_exon_bed_rows.append(bed_row)
    exon_coverages = collect_d4_coverages(
        all_exon_bed_rows, out_dir, "exons", d4tools_command, d4_file, thresholds
    )

    gene_mane_exons_coverage: Dict[str, Coverage] = {}
    for gene in panel_genes:
        if gene.mane_transcript:
            coverage = calculate_exon_coverage(gene, exon_coverages)
            gene_mane_exons_coverage[gene.get_gene_loc()] = coverage

    out_path = out_dir / "results.tsv"
    write_output(
        panel_genes, gene_coverages, gene_mane_exons_coverage, mane_coverages, out_path, thresholds
    )


def parse_arguments():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--panel_tsv", required=True, help="One column with gene panel HGNC symbols"
    )
    parser.add_argument("--mane_gtf", help="GTF file with MANE transcripts", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--d4tools", help="Path to d4tools executable", required=True)
    parser.add_argument("--d4_file", help="Path to d4tools executable", required=True)
    parser.add_argument(
        "--cov_thresholds",
        help="d4tools calculates perc of reads passing these coverage thresholds. Comma separated integers.",
        default="10,30",
    )
    parser.add_argument(
        "--keep_chr", action="store_true", help="Trim leading 'chr' from contig names"
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    if not args.keep_chr:
        print("Trimming leading 'chr', use option '--keep_chr' to keep them")

    print("Loading panel genes")
    panel_hgnc_symbols = parse_panel_text(Path(args.panel_tsv))

    print(f"{len(panel_hgnc_symbols)} HGNC symbols loaded")

    print("Parsing MANE transcripts ...")
    genes = parse_mane_gtf(Path(args.mane_gtf), args.keep_chr)

    thresholds = [int(thres) for thres in args.cov_thresholds.split(",")]

    main(args.d4tools, Path(args.d4_file), panel_hgnc_symbols, genes, Path(args.outdir), thresholds)
