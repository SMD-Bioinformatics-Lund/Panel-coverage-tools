import json
from pathlib import Path
from typing import List, Dict
import gzip

def parse_panel_json(panel_json: Path) -> List[str]:
    hgnc_symbols = []
    with panel_json.open() as json_fh:
        json_data = json.load(json_fh)
        hgnc_symbols = json_data['genes']
    return hgnc_symbols

def parse_panel_text(panel_text: Path) -> List[str]:
    hgnc_symbols = []
    with panel_text.open() as txt_fh:
        for line in txt_fh:
            line = line.rstrip()
            hgnc_symbols.append(line)
    return hgnc_symbols

def parse_mim2gene(mim2gene: Path) -> Dict[str, str]:
    hgnc_to_ensembl = {}

    with open(mim2gene) as in_fh:
        for line in in_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            fields = line.split("\t")
            _mim_number, _mim_entry_type, _entrez_gene_id, hgnc_symbol, ensembl_gene = fields
            hgnc_to_ensembl[hgnc_symbol] = ensembl_gene
            
    return hgnc_to_ensembl


class Gene:
    def __init__(self, gene_id: str):
        self.gene_id = gene_id
        self.entries = []

        self.gene_entry = None
        self.transcript_entry = None
        self.exons = []
    
    def add_info_entry(self, entry_type: str, chr: str, start: int, end: int, info_fields: Dict[str, str]):
        self.entries.append(info_fields)


def parse_gtf_info(info_str: str) -> Dict[str, str]:
    info_fields = [field.strip() for field in info_str.split(";")]
    gtf_info = {}
    for field in info_fields:
        key, value_raw = field.split(" ")
        value_clean = value_raw.strip('"')
        gtf_info[key] = value_clean
    return gtf_info


def parse_mane_gtf(mane_gtf: Path) -> Dict[str, str]:

    if mane_gtf.suffix == "gz":
        fh = gzip.open(mane_gtf, "rt") 
    else:
        fh = mane_gtf.open()
    
    mane_transcripts = {}

    for line in fh:
        line = line.rstrip()
        fields = line.split("\t")
        chr, evidence, mol_type, start, end, _, strand, _, info_str = fields
        parsed_info = parse_gtf_info(info_str)

        gene_id = parsed_info.get("gene_id")
        if not gene_id:
            raise ValueError("FIXME: Should not gene_id be present in all?")
        
        if not mane_transcripts.get(gene_id):
            mane_transcripts[gene_id] = Gene(gene_id)
        
        mane_transcripts[gene_id].add_info_entry(chr, start, end, parsed_info)


    fh.close()


