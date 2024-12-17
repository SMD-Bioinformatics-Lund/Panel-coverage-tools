from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from d4tools_utils import Coverage
from src.parse_utils import Gene


class GeneCoverage:
    def __init__(self, gene: Gene):
        self.gene = gene
        self.exons = gene.mane_transcript.exons
        self.exon_coverage = {}
        self.exon_cov_at_thres: Dict[str, Dict[int, float]] = {}

        self.gene_cov: Optional[float] = None
        self.gene_cov_at_thres: Optional[Dict[int, float]] = None
        self.mane_cov: Optional[float] = None
        self.mane_cov_at_thres: Optional[Dict[int, float]] = None

    def get_gene_loc(self) -> str:
        gene_entry = self.gene.gene_entry
        return f"{gene_entry.chr}_{gene_entry.start}_{gene_entry.end}"

    def get_mane_loc(self) -> str:
        mane_transcript = self.gene.mane_transcript
        return f"{mane_transcript.chr}_{mane_transcript.start}_{mane_transcript.end}"

    def add_coverage(self, loc: str, coverage: float):
        if self.exon_coverage.get(loc) and self.exon_coverage[loc] != coverage:
            raise ValueError("Exon coverage already loaded with different coverage")
        self.exon_coverage[loc] = coverage

    def add_cov_at_thres(self, loc: str, cov_at_thres: Dict[int, float]):
        if self.exon_cov_at_thres.get(loc):
            for thres, cov in cov_at_thres.items():
                if self.exon_cov_at_thres[loc][thres] != cov:
                    raise ValueError("Exon coverage already loaded")
        self.exon_cov_at_thres[loc] = cov_at_thres

    def get_weighted_coverage(self) -> float:
        if len(self.exons) != len(self.exon_coverage):
            import pdb

            pdb.set_trace()
            raise ValueError("Unexpected difference in length in number of exons")

        tot_length = 0
        cov_x_lens = []
        for exon in self.exons:
            exon_len = exon.end - exon.start + 1
            loc = f"{exon.chr}_{exon.start}_{exon.end}"

            cov = self.exon_coverage[loc]
            cov_x_len = cov * exon_len
            tot_length += exon_len
            cov_x_lens.append(cov_x_len)

        return sum(cov_x_lens) / tot_length

    def get_weighted_cov_at_thres(self, thres: int) -> float:
        assert len(self.exons) == len(self.exon_cov_at_thres)

        tot_length = 0
        cov_at_thres_x_lens = []
        for exon in self.exons:
            exon_len = exon.end - exon.start + 1
            loc = f"{exon.chr}_{exon.start}_{exon.end}"

            cov_at_thres = self.exon_cov_at_thres[loc][thres]
            cov_x_len = cov_at_thres * exon_len
            tot_length += exon_len
            cov_at_thres_x_lens.append(cov_x_len)
        return sum(cov_at_thres_x_lens) / tot_length


# Mutates panel_genes by adding exon coverage
def calculate_exon_coverage(
    gene: GeneCoverage, loc_to_exon_cov: Dict[str, Coverage]
) -> Dict[str, float]:

    exons = gene.exons

    total_length = 0
    total_length_cov_products = 0
    total_length_perc_at_thres_products = defaultdict(float)
    for exon in exons:
        loc = f"{exon.chr}_{exon.start}_{exon.end}"
        exon_cov = loc_to_exon_cov[loc]
        exon_len = exon.end - exon.start + 1
        total_length += exon_len

        for thres, thres_cov in exon_cov.perc_at_thres.items():
            total_length_cov_products += total_length * thres_cov
            total_length_perc_at_thres_products[thres] += total_length_cov_products

    results = {"coverage": total_length_cov_products / total_length}

    for thres, thres_cov in total_length_perc_at_thres_products:
        results[f"{thres}x"] = thres_cov / total_length

    return results

    # total_length = [ex.end - ex.start + 1 for ex in exon_coverages]

    # parse_exon_coverage(gene_entries, exons_coverage, exons_cov_at_thres, thresholds)

    # for gene in gene_entries:
    #     exon_cov = gene.get_weighted_coverage()
    #     cov_at_thres = []
    #     for thres in thresholds:
    #         exon_cov_at_thres = gene.get_weighted_cov_at_thres(thres)
    #         cov_at_thres.append(exon_cov_at_thres)
    #     print(f"{gene.gene.hgnc_symbol} {exon_cov} {cov_at_thres}")


def parse_exon_coverage(
    panel_genes: List[GeneCoverage], exon_coverage: Path, exon_thres: Path, thresholds: List[int]
):
    # genes_only: List[GeneCoverage] = []
    # Same exon can be used in multiple genes
    exon_loc_to_genes: Dict[str, List[GeneCoverage]] = defaultdict(list)
    for gene_exon_coverage in panel_genes:
        # gene_exon_coverage = GeneCoverage(gene)
        gene = gene_exon_coverage.gene
        for exon in gene.mane_transcript.exons:
            loc = f"{exon.chr}_{exon.start}_{exon.end}"
            exon_loc_to_genes[loc].append(gene_exon_coverage)
        # genes_only.append(gene_exon_coverage)

    coverage_header = None
    with exon_coverage.open() as in_fh:
        for line in in_fh:
            line = line.rstrip()

            if line == "":
                continue

            if not coverage_header:
                coverage_header = line.split("\t") + ["Coverages"]
                continue

            chr, start, end, cov = line.split("\t")

            loc = f"{chr}_{start}_{end}"

            for gene in exon_loc_to_genes[loc]:
                gene.add_coverage(loc, float(cov))

    thres_header = None
    with exon_thres.open() as in_fh:
        for line in in_fh:
            line = line.rstrip()

            if line == "":
                continue

            if not thres_header:
                thres_header = line.split("\t")
                continue

            fields = line.split("\t")
            chr = fields[0]
            start = fields[1]
            end = fields[2]
            loc = f"{chr}_{start}_{end}"
            cov_at_thresholds = fields[3:]
            cov_at_thres_dict: Dict[int, float] = {}
            for i, thres in enumerate(thresholds):
                cov_at_thres = float(cov_at_thresholds[i])
                cov_at_thres_dict[thres] = cov_at_thres

            for gene in exon_loc_to_genes[loc]:
                gene.add_cov_at_thres(loc, cov_at_thres_dict)
