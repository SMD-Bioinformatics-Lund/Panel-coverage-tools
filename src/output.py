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

    headers = ["hgnc_symbol"]

    for mol_type in ["gene", "mane", "exon"]:
        headers.append(f"{mol_type}_cov")
        for thres in thresholds:
            headers.append(f"{mol_type}_{thres}x")
        headers.append(f"{mol_type}_loc")

    headers.append("exon_all_cov")
    for thres in thresholds:
        headers.append(f"exon_all_{thres}x")

    with out_path.open("w") as out_fh:
        print("\t".join(headers), file=out_fh)
        for gene in panel_genes:

            gene_loc = gene.gtf_entry.get_loc()
            gene_cov = gene_coverages[gene_loc]

            output_row: List[str] = [gene.hgnc_symbol]

            output_row.extend(gene_cov.get_output_vals(thresholds))

            if gene.mane_transcript:
                gene_exon_cov = gene_mane_exons_coverage[gene_loc]
                mane_loc = gene.mane_transcript.get_loc()
                mane_cov = mane_coverages[mane_loc]

                output_row.extend(mane_cov.get_output_vals(thresholds))
                output_row.extend(gene_exon_cov.get_output_vals(thresholds, skip_loc=True))
                output_row.append(gene_exon_cov.location_str)

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
            else:
                output_row.extend(["-"] * (len(headers) - len(output_row)))

            output_string = "\t".join(output_row)
            print(output_string, file=out_fh)
