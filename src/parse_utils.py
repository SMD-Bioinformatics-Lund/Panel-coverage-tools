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


def parse_gtf(gtf: Path, keep_chr: bool, with_progress: bool = False) -> List[GTFEntry]:

    gtf_entries: List[GTFEntry] = []

    if gtf.suffix == ".gz":
        print(f"Found suffix: {gtf.suffix}, parsing as gzip file")
        fh = gzip.open(gtf, "rt")
    else:
        fh = gtf.open()

    nbr_parsed = 0
    for line in fh:
        line = line.rstrip()
        if line.startswith("#!"):
            print(f"Skipping header line: {line}")
            continue

        if not keep_chr and line.startswith("chr"):
            line = line.replace("chr", "", 1)

        gtf_entry = GTFEntry.parse_row(line)
        gtf_entries.append(gtf_entry)

        if with_progress:
            nbr_parsed += 1
            if nbr_parsed % 100000 == 0:
                print(f"Number of rows parsed: {nbr_parsed}")

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
