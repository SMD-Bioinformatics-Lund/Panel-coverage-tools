from pathlib import Path
from typing import Callable, Dict, List

from src.d4tools_utils import calculate_coverage, calculate_perc_at_thres
from src.exons_utils import GeneCoverage


def assign_single_coverage(
    gene_cov_entries: List[GeneCoverage],
    loc_to_gene_cov: Dict[str, GeneCoverage],
    out_dir: Path,
    d4tools_command: str,
    d4_file: Path,
    get_bed_row: Callable[[GeneCoverage], str],
    label: str,
):

    # FIXME
    thresholds = [10, 30]

    regions_bed = out_dir / f"{label}.bed"
    thres_out = out_dir / f"{label}_cov_at_thres.tsv"
    gene_cov_out = out_dir / f"{label}_coverage.tsv"
    with regions_bed.open("w") as out_fh:
        for gene_entry in gene_cov_entries:
            # gene_entry = gene_entry.gene
            print(get_bed_row(gene_entry), file=out_fh)
            # print(gene_entry.get_bed_row(), file=out_fh)

    cov_results = calculate_coverage(d4tools_command, d4_file, regions_bed, gene_cov_out)
    at_thres_results = calculate_perc_at_thres(
        d4tools_command, d4_file, regions_bed, thresholds, thres_out
    )

    for cov_row in cov_results:
        chr, start, end, cov = cov_row.split("\t")
        loc = f"{chr}_{start}_{end}"
        if label == "gene":
            loc_to_gene_cov[loc].gene_cov = float(cov)
        elif label == "mane":
            loc_to_gene_cov[loc].mane_cov = float(cov)
        else:
            raise ValueError(f"Unexpected situation, found label: {label}")

    for cov_at_thres in at_thres_results:
        fields = cov_at_thres.split("\t")
        chr = fields[0]
        start = fields[1]
        end = fields[2]
        # FIXME: Generalize
        cov_10x = fields[3]
        cov_30x = fields[4]
        cov_at_thres = {10: float(cov_10x), 30: float(cov_30x)}
        # chr, start, end, cov = gene_cov_row
        loc = f"{chr}_{start}_{end}"
        if label == "gene":
            loc_to_gene_cov[loc].gene_cov_at_thres = cov_at_thres
        elif label == "mane":
            loc_to_gene_cov[loc].mane_cov_at_thres = cov_at_thres
        else:
            raise ValueError(f"Unexpected situation, found label: {label}")
