from collections import defaultdict
from typing import Dict, List

from src.d4tools_utils import Coverage
from src.parse_utils import Gene


# Mutates panel_genes by adding exon coverage
def calculate_exon_coverage(gene: Gene, loc_to_exon_cov: Dict[str, Coverage]) -> Coverage:

    exons = gene.mane_transcript.exons

    exon_lengths: List[int] = []
    exon_covs: List[float] = []
    exon_cov_at_threS: Dict[int, List[float]] = defaultdict(list)
    # total_length_perc_at_thres_products: Dict[int, float] = defaultdict(float)
    exon_ranges: List[str] = []

    for exon in exons:
        loc = f"{exon.chr}_{exon.start}_{exon.end}"
        exon_cov = loc_to_exon_cov[loc]
        exon_covs.append(exon_cov.cov)

        exon_len = exon.end - exon.start + 1
        exon_lengths.append(exon_len)

        exon_range = f"{exon.start}-{exon.end}"
        exon_ranges.append(exon_range)

        for thres, thres_cov in exon_cov.perc_at_thres.items():
            exon_cov_at_threS[thres].append(thres_cov)

    cov_x_len = sum([exon_covs[i] * exon_lengths[i] for i in range(len(exons))])
    coverage = cov_x_len / sum(exon_lengths)

    cov_at_thres: Dict[int, float] = {}
    for thres, thres_covs in exon_cov_at_threS.items():
        thres_sum = sum([thres_covs[i] * exon_lengths[i] for i in range(len(exons))])
        cov_at_thres[thres] = thres_sum / sum(exon_lengths)

    location = f"{exons[0].chr}:{','.join(exon_ranges)}"

    return Coverage(
        gene.gene_entry.chr,
        gene.gene_entry.start,
        gene.gene_entry.end,
        coverage,
        cov_at_thres,
        location,
    )
