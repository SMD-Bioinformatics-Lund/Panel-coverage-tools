from pathlib import Path
from typing import Dict, List

from d4tools_utils import calculate_coverage, calculate_perc_at_thres
from exons_utils import GeneCoverage


def assign_single_coverage(
    gene_covs: List[GeneCoverage],
    loc_to_gene_cov: Dict[str, GeneCoverage],
    out_dir: Path,
    d4tools_command: str,
    d4_file: Path,
):

    # FIXME
    thresholds = [10, 30]

    gene_bed = out_dir / "genes.bed"
    gene_thres_out = out_dir / "gene_cov_at_thres.tsv"
    gene_cov_out = out_dir / "gene_coverage.tsv"
    with gene_bed.open("w") as out_fh:
        for gene_cov in gene_covs:
            gene = gene_cov.gene
            print(gene.get_bed_row(), file=out_fh)

    gene_cov_results = calculate_coverage(d4tools_command, d4_file, gene_bed, gene_cov_out)
    gene_at_thres_results = calculate_perc_at_thres(
        d4tools_command, d4_file, gene_bed, thresholds, gene_thres_out
    )

    for gene_cov_row in gene_cov_results:
        chr, start, end, cov = gene_cov_row.split("\t")
        loc = f"{chr}_{start}_{end}"
        loc_to_gene_cov[loc].gene_cov = float(cov)

    for gene_cov_at_thres in gene_at_thres_results:
        fields = gene_cov_at_thres.split("\t")
        chr = fields[0]
        start = fields[1]
        end = fields[2]
        # FIXME: Generalize
        cov_10x = fields[3]
        cov_30x = fields[4]
        cov_at_thres = {10: float(cov_10x), 30: float(cov_30x)}
        # chr, start, end, cov = gene_cov_row
        loc = f"{chr}_{start}_{end}"
        loc_to_gene_cov[loc].gene_cov_at_thres = cov_at_thres
