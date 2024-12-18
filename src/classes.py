from collections import defaultdict
from typing import List, Dict, Optional


MANE_SELECT_TAG = "MANE_Select"
MANE_SELECT_PLUS_CLINICAL_TAG = "MANE_Plus_Clinical"
MANE_TAGS = [MANE_SELECT_TAG, MANE_SELECT_PLUS_CLINICAL_TAG]


class GTFEntry:
    def __init__(
        self,
        raw_line: str,
        mol_type: str,
        chr: str,
        start: int,
        end: int,
        gene_id: str,
        hgnc_name: Optional[str],
        tags: List[str],
        transcript_id: Optional[str],
        db_xref: Optional[str],
    ):
        self.raw_line = raw_line
        self.mol_type = mol_type
        self.chr = chr
        self.start = start
        self.end = end

        self.gene_id = gene_id
        self.hgnc_name = hgnc_name
        self.transcript_id = transcript_id
        self.db_xref = db_xref
        self.tags = tags

    def has_mane_tag(self) -> bool:
        mane_tags = [tag for tag in self.tags if tag in MANE_TAGS]
        return len(mane_tags) > 0

    def get_loc(self) -> str:
        return f"{self.chr}:{self.start}-{self.end}"

    @staticmethod
    def parse_row(row: str) -> "GTFEntry":
        fields = row.split("\t")
        line = row.rstrip()
        fields = line.split("\t")
        raw_chr, _, mol_type, start, end, _, _, _, info_str = fields
        chr = raw_chr
        parsed_info = GTFEntry._parse_gtf_info(info_str)

        gene_id = parsed_info["gene_id"][0]
        hgnc_name = parsed_info["gene_name"][0] if parsed_info.get("gene_name") else None

        tags = parsed_info.get("tag")
        transcript_id = parsed_info.get("transcript_id")
        db_xref = parsed_info.get("db_xref")

        return GTFEntry(
            line,
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
        gtf_info: Dict[str, List[str]] = defaultdict(list)
        for field in info_fields:
            try:
                key, value_raw = field.split(" ", 1)
            except:
                raise ValueError(f"Failing to parse for fields: {info_fields}")
            value_clean = value_raw.strip('"')
            gtf_info[key].append(value_clean)
        return gtf_info

    def __str__(self) -> str:
        return f"{self.gene_id} {self.hgnc_name} {self.mol_type} {self.chr} {self.start} {self.end}"


class Transcript:
    def __init__(
        self,
        ensembl_transcript_id: str,
        transcript_entry: GTFEntry,
        chr: str,
        start: int,
        end: int,
        db_xref: Optional[str],
        mane_tag: str,
        exon_entries: List[GTFEntry],
    ):
        self.ensembl_id = ensembl_transcript_id.split(".")[0]
        self.gtf_entry = transcript_entry
        self.chr = chr
        self.start = start
        self.end = end
        self.db_xref = db_xref
        self.mane_tag = mane_tag
        self.exons = exon_entries

    def get_loc(self) -> str:
        return self.gtf_entry.get_loc()

    def get_bed_exons(self) -> List[str]:
        bed_exons: List[str] = []
        for exon in self.exons:
            bed_exon = "\t".join([exon.chr, str(exon.start), str(exon.end)])
            bed_exons.append(bed_exon)
        return bed_exons

    @staticmethod
    def parse_gtf_entries(gtf_entries: List[GTFEntry]) -> "Transcript":
        transcript_entry = gtf_entries[0]

        if transcript_entry.mol_type != "transcript":
            raise ValueError("First entry expected to be a transcript entry")

        exon_entries: List[GTFEntry] = []

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

        if len(mane_tags) == 0:
            raise ValueError("Expected mane_tag missing")

        if len(mane_tags) > 1:
            raise ValueError(f"No more than one mane tag expected, found: {mane_tags}")

        return Transcript(
            transcript_id,
            transcript_entry,
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
        mane_transcript: Optional[Transcript],
        mane_plus_clinical_transcript: Optional[Transcript],
    ):
        self.ensembl_gene_id = ensembl_gene_id.split(".")[0]
        self.hgnc_symbol = hgnc_symbol
        self.gtf_entry = gene
        self.mane_transcript = mane_transcript
        self.mane_plus_clinical_transcript = mane_plus_clinical_transcript

    def get_bed_row(self) -> str:
        return "\t".join([self.gtf_entry.chr, str(self.gtf_entry.start), str(self.gtf_entry.end)])

    def get_mane_transcript_bed_row(self) -> str:
        if not self.mane_transcript:
            raise ValueError("No MANE transcript available")
        return "\t".join(
            [
                self.mane_transcript.chr,
                str(self.mane_transcript.start),
                str(self.mane_transcript.end),
            ]
        )

    def get_gene_loc(self) -> str:
        # gene = self.gene_entry
        return self.gtf_entry.get_loc()

    def get_mane_loc(self) -> str:
        if not self.mane_transcript:
            raise ValueError("No MANE transcript available")
        return self.mane_transcript.get_loc()

    def get_bed_exons(self) -> List[str]:
        if not self.mane_transcript:
            raise ValueError("No MANE transcript available")
        return self.mane_transcript.get_bed_exons()

    @staticmethod
    def parse_gtf_entries(gtf_entries: List[GTFEntry]) -> "Gene":

        gene_entry = None
        mane_transcript: Optional[Transcript] = None
        mane_plus_transcript: Optional[Transcript] = None

        gene_entry = gtf_entries[0]
        if gene_entry.mol_type != "gene":
            raise ValueError(f"First entry expected to be a gene entry")

        transcript_entries: Dict[str, List[GTFEntry]] = defaultdict(list)

        # transcript_entries: Dict[str, List[GTFEntry]] = defaultdict(list)
        for gtf_entry in gtf_entries[1:]:
            if not gtf_entry.transcript_id:
                raise ValueError(
                    f"Expected (but did not find) transcript_id for GTF entry: {gtf_entry}"
                )
            transcript_entries[gtf_entry.transcript_id].append(gtf_entry)

        for entries in transcript_entries.values():

            if entries[0].mol_type != "transcript":
                raise ValueError(
                    f"First entry expected to be transcript, found: {entries[0].mol_type}"
                )

            for entry in entries[1:]:
                if entry.mol_type == "gene":
                    raise ValueError("Gene already found, something is wrong")
                if entry.mol_type == "transcript":
                    raise ValueError("Transcript already found, something is wrong")

            transcript = Transcript.parse_gtf_entries(entries)
            if transcript.mane_tag == MANE_SELECT_TAG:
                mane_transcript = transcript
            elif transcript.mane_tag == MANE_SELECT_PLUS_CLINICAL_TAG:
                mane_plus_transcript = transcript
            else:
                raise ValueError(f"Unknown MANE tag: {transcript.mane_tag}")

        if mane_transcript is None:
            print(f"No MANE found for {gene_entry.hgnc_name}")
            # raise ValueError(f"No MANE transcript entry found for {gene_entry.hgnc_name}")

        if gene_entry.hgnc_name is None:
            raise ValueError(f"No HGNC entry name found for gene entry: {gene_entry}")

        return Gene(
            gene_entry.gene_id,
            gene_entry.hgnc_name,
            gene_entry,
            mane_transcript,
            mane_plus_transcript,
        )

    def __str__(self) -> str:
        return f"{self.hgnc_symbol} {self.get_gene_loc()} {self.get_mane_loc()}"
