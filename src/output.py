from pathlib import Path
from typing import Dict, List

from src.d4tools_utils import Coverage
from src.parse_utils import Gene


def write_output(
    panel_genes: List[Gene],
    gene_coverages: Dict[str, Coverage],
    gene_mane_exons_coverage: Dict[str, Coverage],
    mane_coverages: Dict[str, Coverage],
    out_path: Path,
    thresholds: List[int],
):

    # FIXME: OK, let's think about how to refactor this

    headers = ["hgnc_symbol"]

    headers.append("gene_cov")
    for thres in thresholds:
        headers.append(f"gene_{thres}x")
    headers.append("gene")

    headers.append("mane_cov")
    for thres in thresholds:
        headers.append(f"mane_{thres}x")
    headers.append("mane")

    headers.append("exon_cov")
    for thres in thresholds:
        headers.append(f"exon_{thres}x")
    headers.append("exon_all_cov")
    for thres in thresholds:
        headers.append(f"exon_all_{thres}x")
    headers.append("exons")

    with out_path.open("w") as out_fh:
        print("\t".join(headers), file=out_fh)
        for gene in panel_genes:

            gene_loc = gene.gtf_entry.get_loc()
            gene_cov = gene_coverages[gene_loc]
            gene_exon_cov = gene_mane_exons_coverage[gene_loc]

            mane_loc = gene.mane_transcript.get_loc()
            mane_cov = mane_coverages[mane_loc]

            output_row: List[str] = [gene.hgnc_symbol]

            output_row.append(str(gene_cov.cov))
            for cov_at_thres in gene_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))
            output_row.append(gene_cov.location_str)

            output_row.append(str(mane_cov.cov))
            for cov_at_thres in mane_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))
            output_row.append(mane_cov.location_str)

            output_row.append(str(gene_exon_cov.cov))
            for cov_at_thres in gene_exon_cov.perc_at_thres.values():
                output_row.append(str(cov_at_thres))

            if not gene_exon_cov.exons:
                raise ValueError(
                    "Unexpected situation, exon coverage should have individual exon coverage"
                )

            exon_covs = ",".join([str(exon_cov.cov) for exon_cov in gene_exon_cov.exons])
            output_row.append(exon_covs)
            for thres in thresholds:
                exon_perc_at_thres = ",".join(
                    [str(exon_cov.perc_at_thres[thres]) for exon_cov in gene_exon_cov.exons]
                )
                output_row.append(str(exon_perc_at_thres))
            output_row.append(gene_exon_cov.location_str)

            output_string = "\t".join(output_row)
            print(output_string, file=out_fh)
