from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from src.parse_utils import Gene


class GeneExonCoverage:
    def __init__(self, gene: Gene):
        self.gene = gene
        self.exons = gene.mane_transcript.exons
        self.exon_coverage = {}
        self.exon_cov_at_thres: Dict[str, Dict[int, float]] = {}

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


def calculate_exon_coverage(
    panel_genes: List[Gene], exons_coverage: Path, exons_cov_at_thres: Path, thresholds: List[int]
):
    genes_with_exon_coverage = parse_exon_coverage(
        panel_genes, exons_coverage, exons_cov_at_thres, thresholds
    )

    for gene in genes_with_exon_coverage:
        exon_cov = gene.get_weighted_coverage()
        cov_at_thres = []
        for thres in thresholds:
            exon_cov_at_thres = gene.get_weighted_cov_at_thres(thres)
            cov_at_thres.append(exon_cov_at_thres)
        # print(f"{gene.gene.hgnc_symbol} {exon_cov} {cov_at_thres}")


def parse_exon_coverage(
    panel_genes: List[Gene], exon_coverage: Path, exon_thres: Path, thresholds: List[int]
) -> List[GeneExonCoverage]:
    genes_only: List[GeneExonCoverage] = []
    # Same exon can be used in multiple genes
    exon_loc_to_genes: Dict[str, List[GeneExonCoverage]] = defaultdict(list)
    for gene in panel_genes:
        gene_exon_coverage = GeneExonCoverage(gene)
        for exon in gene.mane_transcript.exons:
            loc = f"{exon.chr}_{exon.start}_{exon.end}"
            exon_loc_to_genes[loc].append(gene_exon_coverage)
        genes_only.append(gene_exon_coverage)

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
    return genes_only
