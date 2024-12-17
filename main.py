#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Callable, List, Dict, Set

from src.exons_utils import calculate_exon_coverage, parse_exon_coverage
from src.d4tools_utils import calculate_coverage, calculate_perc_at_thres
from src.parse_utils import Gene, parse_mane_gtf, parse_mim2gene, parse_panel_json, parse_panel_text

__version_info__ = ("1", "0", "1")
__version__ = ".".join(__version_info__)

description = "Placeholder description"


def main(
    d4tools_command: str,
    d4_file: Path,
    panel_hgnc_symbols: Set[str],
    all_genes: List[Gene],
    out_dir: Path,
    force: bool,
):

    panel_genes = [gene for gene in all_genes if gene.hgnc_symbol in panel_hgnc_symbols]
    found_symbols = [gene.hgnc_symbol for gene in panel_genes]

    print(f"Now I have {len(panel_genes)} omim genes")

    missed_hgnc_symbols = [
        hgnc_symbol for hgnc_symbol in panel_hgnc_symbols if hgnc_symbol not in found_symbols
    ]

    print(f"Missed {len(missed_hgnc_symbols)} HGNC symbols: {missed_hgnc_symbols}")
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_bed = out_dir / "genes.bed"
    d4tools_out_dir = out_dir / "d4tools_results"
    d4tools_out_dir.mkdir(parents=True, exist_ok=True)
    thresholds = [10, 30]

    # FIXME: Refactor this

    with gene_bed.open("w") as out_fh:
        for gene in panel_genes:
            print(gene.get_bed_row(), file=out_fh)
    gene_cov_out = d4tools_out_dir / "gene_coverage.tsv"
    calculate_coverage(d4tools_command, d4_file, gene_bed, gene_cov_out)
    gene_thres_out = d4tools_out_dir / "gene_cov_at_thres.tsv"
    calculate_perc_at_thres(d4tools_command, d4_file, gene_bed, thresholds, gene_thres_out)

    mane_transcripts_bed = out_dir / "mane_transcripts.bed"
    with mane_transcripts_bed.open("w") as out_fh:
        for gene in panel_genes:
            print(gene.get_mane_transcript_bed_row(), file=out_fh)
    mane_cov_out = d4tools_out_dir / "mane_coverage.tsv"
    calculate_coverage(d4tools_command, d4_file, mane_transcripts_bed, mane_cov_out)
    mane_thres_out = d4tools_out_dir / "mane_cov_at_thres.tsv"
    calculate_perc_at_thres(d4tools_command, d4_file, gene_bed, thresholds, mane_thres_out)

    mane_exons_bed = out_dir / "mane_exons.bed"
    with mane_exons_bed.open("w") as out_fh:
        all_exons = set()
        for gene in panel_genes:
            for bed_row in gene.get_bed_exons():
                all_exons.add(bed_row)
        # Unique rows - only calculate coverage once
        for bed_row in all_exons:
            print(bed_row, file=out_fh)
    exons_cov_out = d4tools_out_dir / "exons_coverage.tsv"
    calculate_coverage(d4tools_command, d4_file, mane_exons_bed, exons_cov_out)
    exons_thres_out = d4tools_out_dir / "mane_cov_at_thres.tsv"
    calculate_perc_at_thres(d4tools_command, d4_file, mane_exons_bed, thresholds, exons_thres_out)

    # FIXME: Parsing back of the exons
    calculate_exon_coverage(panel_genes, exons_cov_out, exons_thres_out, thresholds)


def parse_arguments():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--panel_json", help="Either --panel_json or --panel_txt needs to be supplied"
    )
    parser.add_argument(
        "--panel_text", help="Either --panel_json or --panel_txt needs to be supplied"
    )
    parser.add_argument(
        "--mim2gene", help="OMIM map of gene names to ENSEMBL gene IDs", required=True
    )
    parser.add_argument("--mane_gtf", help="GTF file with MANE transcripts", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--d4tools", help="Path to d4tools executable", required=True)
    parser.add_argument("--d4_file", help="Path to d4tools executable", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--keep_chr", action="store_true", help="Trim leading 'chr' from contig names"
    )

    args = parser.parse_args()

    if args.panel_json is None and args.panel_text is None:
        raise ValueError("Either --panel_json or --panel_text needs to be supplied")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    if not args.keep_chr:
        print("Trimming leading 'chr', use option '--keep_chr' to keep them")

    print("Loading panel genes")
    if args.panel_text is None and args.panel_json is not None:
        panel_hgnc_symbols = parse_panel_json(Path(args.panel_json))
    else:
        panel_hgnc_symbols = parse_panel_text(Path(args.panel_text))

    print(f"{len(panel_hgnc_symbols)} HGNC symbols loaded")

    hgnc_to_ensembl = parse_mim2gene(args.mim2gene)

    print("Parsing MANE transcripts ...")
    genes = parse_mane_gtf(Path(args.mane_gtf), args.keep_chr, args.verbose)

    main(args.d4tools, Path(args.d4_file), panel_hgnc_symbols, genes, Path(args.outdir), args.force)
