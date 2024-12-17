from collections import defaultdict
from typing import Dict

from src.d4tools_utils import Coverage
from src.parse_utils import Gene


# Mutates panel_genes by adding exon coverage
def calculate_exon_coverage(gene: Gene, loc_to_exon_cov: Dict[str, Coverage]) -> Coverage:

    exons = gene.mane_transcript.exons

    total_length = 0
    total_length_cov_products = 0
    total_length_perc_at_thres_products: Dict[int, float] = defaultdict(float)
    for exon in exons:
        loc = f"{exon.chr}_{exon.start}_{exon.end}"
        exon_cov = loc_to_exon_cov[loc]
        exon_len = exon.end - exon.start + 1
        total_length += exon_len

        for thres, thres_cov in exon_cov.perc_at_thres.items():
            total_length_cov_products += exon_len * thres_cov
            total_length_perc_at_thres_products[thres] += total_length_cov_products

    coverage = total_length_cov_products / total_length

    cov_at_thres: Dict[int, float] = {}
    for thres, thres_cov in total_length_perc_at_thres_products.items():
        cov_at_thres[thres] = thres_cov / total_length

    return Coverage(
        gene.gene_entry.chr, gene.gene_entry.start, gene.gene_entry.end, coverage, cov_at_thres
    )
