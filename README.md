Tools for working with coverage and panels in d4 files.

# Preparing input data

## Panel data

The input panel is expected to be a flat text file with HGNC names for each gene of interest.

If using Scout, this can be retrieved from the mongo database as such.

```
mongoexport \
    --uri="mongodb://localhost:27017/scout" \
    --collection=gene_panel \
    --query='{"panel_name": "OMIM-AUTO"}' \
    --sort='{ "version": -1 }' \
    --limit=1 \
    --out=omim_latest.json
```

This can subsequently be converted using the util script `prepare_panel.py`. Note the `python -m scripts.prepare_panel` syntax is needed for it to import other modules correctly.

```
python3 -m scripts.prepare_panel --panel_json omim_latest.json --out_tsv omim_latest.tsv
```

## Preparing the input GTF

A GTF file with the genes, their annotations and mane transcripts + exons is required.

You can either download a MANE GTF directly, but you might miss out on some transcripts which have no MANE-transcripts. If included, you will still get gene-level coverage.

Or, you download a full GTF and optionally trim it down to panel genes. This will speed up analysis.

### Preparing a GTF from full ENSEMBL data

In the moment of writing, the latest ENSEMBL GTF can be downloaded from:

https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

This can be parsed:

```
python3 -m scripts.prepare_gtf \
    --in_gtf Homo_sapiens.GRCh38.113.gtf.gz \
    --panel_genes omim_latest.tsv \
    --out_gtf parsed_mane_from_full.gtf
```

### Downloading MANE only

The latest MANE version can be downloaded from: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/

Download the genomic GTF file. In the moment of writing, it is this one:

https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz

# Running the main script

```
python3 main.py --help
```
