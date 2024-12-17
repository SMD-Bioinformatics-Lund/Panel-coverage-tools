import json
from pathlib import Path
from typing import List, Dict, Optional, Set
import gzip
from collections import defaultdict

MANE_SELECT_TAG = "MANE_Select"
MANE_SELECT_PLUS_CLINICAL_TAG = "MANE_Plus_Clinical"


def parse_panel_json(panel_json: Path) -> Set[str]:
    hgnc_symbols: Set[str] = set()
    with panel_json.open() as json_fh:
        json_data = json.load(json_fh)
        hgnc_symbols = set([gene_info["symbol"] for gene_info in json_data["genes"]])
    return hgnc_symbols


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


class GTFEntry:
    def __init__(
        self,
        mol_type: str,
        chr: str,
        start: int,
        end: int,
        gene_id: str,
        hgnc_name: str,
        tags: List[str],
        transcript_id: Optional[str],
        db_xref: Optional[str],
    ):
        self.mol_type = mol_type
        self.chr = chr
        self.start = start
        self.end = end

        self.gene_id = gene_id
        self.hgnc_name = hgnc_name
        self.transcript_id = transcript_id
        self.db_xref = db_xref
        self.tags = tags

    @staticmethod
    def parse_row(row: str, keep_chr: bool) -> "GTFEntry":
        fields = row.split("\t")
        line = row.rstrip()
        fields = line.split("\t")
        raw_chr, evidence_, mol_type, start, end, _, strand_, _, info_str = fields
        if not keep_chr and raw_chr.find("chr") == 0:
            chr = raw_chr.replace("chr", "", 1)
        else:
            chr = raw_chr
        parsed_info = GTFEntry._parse_gtf_info(info_str)
        gene_id = parsed_info["gene_id"][0]
        hgnc_name = parsed_info["gene_name"][0]

        tags = parsed_info.get("tag")
        transcript_id = parsed_info.get("transcript_id")
        db_xref = parsed_info.get("db_xref")

        return GTFEntry(
            mol_type,
            chr,
            int(start),
            int(end),
            gene_id,
            hgnc_name,
            tags or [],
            transcript_id[0] if transcript_id else None,
            db_xref[0] if db_xref else None,
        )

    @staticmethod
    def _parse_gtf_info(info_str: str) -> Dict[str, List[str]]:
        info_fields = [field.strip() for field in info_str.strip(";").split(";")]
        gtf_info = defaultdict(list)
        for field in info_fields:
            key, value_raw = field.split(" ")
            value_clean = value_raw.strip('"')
            gtf_info[key].append(value_clean)
        return gtf_info

    def __str__(self) -> str:
        return f"{self.gene_id} {self.hgnc_name} {self.mol_type} {self.chr} {self.start} {self.end}"


class Transcript:
    def __init__(
        self,
        ensembl_transcript_id: str,
        chr: str,
        start: int,
        end: int,
        db_xref: str,
        mane_tag: str,
        exon_entries: List[GTFEntry],
    ):
        self.ensembl_id = ensembl_transcript_id.split(".")[0]
        self.chr = chr
        self.start = start
        self.end = end
        self.db_xref = db_xref
        self.mane_tag = mane_tag
        self.exons = exon_entries

    def get_bed_exons(self) -> List[str]:
        bed_exons = []
        for exon in self.exons:
            bed_exon = "\t".join([exon.chr, str(exon.start), str(exon.end)])
            bed_exons.append(bed_exon)
        return bed_exons

    @staticmethod
    def parse_gtf_entries(gtf_entries: List[GTFEntry]) -> "Transcript":
        transcript_entry = gtf_entries[0]

        if transcript_entry.mol_type != "transcript":
            raise ValueError("First entry expected to be a transcript entry")

        exon_entries = []

        for gtf_entry in gtf_entries[1:]:
            if gtf_entry.mol_type == "transcript":
                raise ValueError("Transcript type already found, something is wrong")

            if gtf_entry.mol_type == "exon":
                exon_entries.append(gtf_entry)

        transcript_id = transcript_entry.transcript_id
        db_xref = transcript_entry.db_xref

        mane_tag_patterns = {MANE_SELECT_TAG, MANE_SELECT_PLUS_CLINICAL_TAG}

        mane_tags = [tag for tag in transcript_entry.tags if tag in mane_tag_patterns]

        if transcript_id is None:
            raise ValueError("Expected transcript_id missing")

        if db_xref is None:
            raise ValueError("Expected RefSeq ID missing")

        if mane_tags is None:
            raise ValueError("Expected mane_tag missing")

        if len(mane_tags) > 1:
            raise ValueError(f"No more than one mane tag expected, found: {mane_tags}")

        return Transcript(
            transcript_id,
            transcript_entry.chr,
            transcript_entry.start,
            transcript_entry.end,
            db_xref,
            mane_tags[0],
            exon_entries,
        )


class Gene:
    def __init__(
        self,
        ensembl_gene_id: str,
        hgnc_symbol: str,
        gene: GTFEntry,
        mane_transcript: Transcript,
        mane_plus_clinical_transcript: Optional[Transcript],
    ):
        self.ensembl_gene_id = ensembl_gene_id.split(".")[0]
        self.hgnc_symbol = hgnc_symbol
        self.gene_entry = gene
        self.mane_transcript = mane_transcript
        self.mane_plus_clinical_transcript = mane_plus_clinical_transcript

    def get_bed_row(self) -> str:
        return "\t".join(
            [self.gene_entry.chr, str(self.gene_entry.start), str(self.gene_entry.end)]
        )

    def get_mane_transcript_bed_row(self) -> str:
        return "\t".join(
            [
                self.mane_transcript.chr,
                str(self.mane_transcript.start),
                str(self.mane_transcript.end),
            ]
        )

    def get_gene_loc(self) -> str:
        gene = self.gene_entry
        return f"{gene.chr}_{gene.start}_{gene.end}"

    def get_mane_loc(self) -> str:
        mane = self.mane_transcript
        return f"{mane.chr}_{mane.start}_{mane.end}"

    def get_bed_exons(self) -> List[str]:
        return self.mane_transcript.get_bed_exons()

    @staticmethod
    def parse_gtf_entries(gtf_entries: List[GTFEntry], verbose: bool) -> "Gene":

        gene_entry = None
        mane_transcript: Optional[Transcript] = None
        mane_plus_transcript: Optional[Transcript] = None

        gene_entry = gtf_entries[0]
        if gene_entry.mol_type != "gene":
            raise ValueError(f"First entry expected to be a gene entry")

        transcript_entries = []
        for gtf_entry in gtf_entries[1:]:
            if gtf_entry.mol_type == "gene":
                raise ValueError("Gene type already found, something is wrong")
            if gtf_entry.mol_type == "transcript":
                if len(transcript_entries) > 0:
                    transcript = Transcript.parse_gtf_entries(transcript_entries)
                    if transcript.mane_tag == MANE_SELECT_TAG:
                        mane_transcript = transcript
                    elif transcript.mane_tag == MANE_SELECT_PLUS_CLINICAL_TAG:
                        mane_plus_transcript = transcript
                    else:
                        raise ValueError(f"Unknown MANE tag: {transcript.mane_tag}")
                    transcript_entries = []
            transcript_entries.append(gtf_entry)

        # FIXME: Refactor dual blocks here
        if len(transcript_entries) > 0:
            transcript = Transcript.parse_gtf_entries(transcript_entries)
            if transcript.mane_tag == MANE_SELECT_TAG:
                mane_transcript = transcript
            elif transcript.mane_tag == MANE_SELECT_PLUS_CLINICAL_TAG:
                mane_plus_transcript = transcript
            else:
                raise ValueError(f"Unknown MANE tag: {transcript.mane_tag}")

        if gene_entry is None:
            raise ValueError("No gene entry found")

        if mane_transcript is None:
            raise ValueError("No MANE transcript entry found")

        return Gene(
            gene_entry.gene_id,
            gene_entry.hgnc_name,
            gene_entry,
            mane_transcript,
            mane_plus_transcript,
        )


def parse_mane_gtf(mane_gtf: Path, keep_chr: bool, verbose: bool) -> List[Gene]:

    print(f"Found suffix: {mane_gtf.suffix}")
    if mane_gtf.suffix == ".gz":
        fh = gzip.open(mane_gtf, "rt")
    else:
        fh = mane_gtf.open()

    genes: List[Gene] = []
    curr_gene_entries: List[GTFEntry] = []
    for line in fh:
        line = line.rstrip()

        gtf_entry = GTFEntry.parse_row(line, keep_chr)

        if gtf_entry.mol_type == "gene":
            if len(curr_gene_entries) > 0:
                if verbose:
                    print(f"Found gene {gtf_entry.hgnc_name} with {len(curr_gene_entries)} entries")
                gene = Gene.parse_gtf_entries(curr_gene_entries, verbose)
                if gene.gene_entry.chr.find("_") == -1:
                    genes.append(gene)
                else:
                    print(
                        f"Skipping gene with non-normal chr: {gene.hgnc_symbol} {gene.gene_entry.chr}"
                    )
                curr_gene_entries = []
        curr_gene_entries.append(gtf_entry)

    gene = Gene.parse_gtf_entries(curr_gene_entries, verbose)
    genes.append(gene)

    fh.close()

    return genes
