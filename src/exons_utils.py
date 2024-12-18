from collections import defaultdict
from typing import Dict, List

from src.d4tools_utils import Coverage
from src.parse_utils import Gene


def calculate_exon_coverage(gene: Gene, loc_to_exon_cov: Dict[str, Coverage]) -> Coverage:
    """
    Calculate average coverage and % cov at threshold for exons
    Coverage is weighted by length of exons
    """

    exons = gene.mane_transcript.exons

    exon_coverages: List[Coverage] = []
    exon_lengths: List[int] = []
    exon_cov_at_threS: Dict[int, List[float]] = defaultdict(list)

    for exon in exons:
        loc = f"{exon.chr}:{exon.start}-{exon.end}"
        exon_cov = loc_to_exon_cov[loc]
        exon_coverages.append(exon_cov)
        exon_len = exon.end - exon.start + 1
        exon_lengths.append(exon_len)

        for thres, thres_cov in exon_cov.perc_at_thres.items():
            exon_cov_at_threS[thres].append(thres_cov)

    cov_x_len = sum([exon_coverages[i].cov * exon_lengths[i] for i in range(len(exons))])
    coverage = cov_x_len / sum(exon_lengths)

    cov_at_thres: Dict[int, float] = {}
    for thres, thres_covs in exon_cov_at_threS.items():
        thres_sum = sum([thres_covs[i] * exon_lengths[i] for i in range(len(exons))])
        cov_at_thres[thres] = thres_sum / sum(exon_lengths)

    exon_ranges = [f"{exon.start}-{exon.end}" for exon in exon_coverages]
    location = f"{exons[0].chr}:{','.join(exon_ranges)}"

    return Coverage(
        gene.gtf_entry.chr,
        gene.gtf_entry.start,
        gene.gtf_entry.end,
        coverage,
        cov_at_thres,
        location,
        exon_coverages,
    )
