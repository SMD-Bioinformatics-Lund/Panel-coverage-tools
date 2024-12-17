#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List, Dict, Set

from src.exons_utils import calculate_exon_coverage
from src.d4tools_utils import (
    Coverage,
    collect_d4_coverages,
)
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
        for bed_row in gene.get_bed_exons():
            all_exon_bed_rows.append(bed_row)
    exon_coverages = collect_d4_coverages(
        all_exon_bed_rows, out_dir, "exons", d4tools_command, d4_file, thresholds
    )

    gene_mane_exons_coverage: Dict[str, Coverage] = {}
    for gene in panel_genes:
        coverage = calculate_exon_coverage(gene, exon_coverages)
        gene_mane_exons_coverage[gene.get_gene_loc()] = coverage

    headers = ["hgnc_symbol"]

    headers.append("gene_cov")
    for thres in thresholds:
        headers.append(f"gene_${thres}x")

    headers.append("mane_cov")
    for thres in thresholds:
        headers.append(f"mane_${thres}x")

    headers.append("exon_cov")
    for thres in thresholds:
        headers.append(f"exon_${thres}x")

    out_path = out_dir / "results.tsv"
    with out_path.open("w") as out_fh:
        print("\t".join(headers), file=out_fh)
        # Collect and print
        # output_rows: List[str] = []
        for gene in panel_genes:
            gene_loc = gene.get_gene_loc()
            gene_cov = gene_coverages[gene_loc]
            gene_exon_cov = gene_mane_exons_coverage[gene_loc]

            mane_loc = gene.get_mane_loc()
            mane_cov = mane_coverages[mane_loc]

            output_row: List[str] = [gene.hgnc_symbol]

            output_row.append(str(gene_cov.cov))
            for cov_at_thres in gene_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))

            output_row.append(str(mane_cov.cov))
            for cov_at_thres in mane_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))

            output_row.append(str(gene_exon_cov.cov))
            for cov_at_thres in gene_exon_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))

            output_string = "\t".join(output_row)
            print(output_string, file=out_fh)
        # print(output_rows)


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
    parser.add_argument(
        "--cov_thresholds",
        help="d4tools calculates perc of reads passing these coverage thresholds. Comma separated integers.",
        default="10,30",
    )
    parser.add_argument("--verbose", action="store_true")
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

    thresholds = [int(thres) for thres in args.cov_thresholds.split(",")]

    main(args.d4tools, Path(args.d4_file), panel_hgnc_symbols, genes, Path(args.outdir), thresholds)
