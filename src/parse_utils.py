from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Set
import gzip

from src.classes import GTFEntry, Gene


def parse_panel_text(panel_text: Path) -> Set[str]:
    hgnc_symbols: Set[str] = set()
    with panel_text.open() as txt_fh:
        for line in txt_fh:
            line = line.rstrip()
            hgnc_symbols.add(line)
    return hgnc_symbols


def parse_mim2gene(mim2gene: Path) -> Dict[str, str]:
    hgnc_to_ensembl: Dict[str, str] = {}

    with open(mim2gene) as in_fh:
        nbr_skipped = 0
        for line in in_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            fields = line.split("\t")
            if len(fields) < 5:
                nbr_skipped += 1
                continue
            _, _, _, hgnc_symbol, ensembl_gene = fields
            hgnc_to_ensembl[hgnc_symbol] = ensembl_gene
        print(f"Number skipped: {nbr_skipped} nbr stored: {len(hgnc_to_ensembl)}")

    return hgnc_to_ensembl


def parse_gtf(gtf: Path, keep_chr: bool) -> List[GTFEntry]:

    gtf_entries: List[GTFEntry] = []

    if gtf.suffix == ".gz":
        print(f"Found suffix: {gtf.suffix}, parsing as gzip file")
        fh = gzip.open(gtf, "rt")
    else:
        fh = gtf.open()

    for line in fh:
        line = line.rstrip()
        if line.startswith("#!"):
            print(f"Skipping header line: {line}")
            continue

        if not keep_chr and line.startswith("chr"):
            line = line.replace("chr", "", 1)

        gtf_entry = GTFEntry.parse_row(line)
        gtf_entries.append(gtf_entry)

    fh.close()

    return gtf_entries


def parse_mane_gtf(mane_gtf: Path, keep_chr: bool) -> List[Gene]:

    gtf_entries = parse_gtf(mane_gtf, keep_chr)

    genes: List[Gene] = []

    gene_to_gtf_entries: Dict[str, List[GTFEntry]] = defaultdict(list)

    for gtf_entry in gtf_entries:
        gene_to_gtf_entries[gtf_entry.gene_id].append(gtf_entry)

    for gene_entries in gene_to_gtf_entries.values():
        gene = Gene.parse_gtf_entries(gene_entries)
        genes.append(gene)

    return genes
