Tools for working with coverage and panels in d4 files.

# Preparing input data

## Panel data

The input panel is expected to be a flat text file with HGNC names for each gene of interest.

If using Scout, this can be retrieved from the mongo database as such:

```
mongoexport \
    --uri="mongodb://localhost:27017/scout" \
    --collection=gene_panel \
    --query='{"panel_name": "OMIM-AUTO"}' \
    --sort='{ "version": -1 }' \
    --out=omim_latest.json
```

You can run the coverage calculation either using this import or by directly supplying a flat text file.

## Mapping HGNC names to ENSEMBL gene IDs

The officially supported map by OMIM can be downloaded from here: https://www.omim.org/static/omim/data/mim2gene.txt

## MANE and exon data

The latest MANE version can be downloaded from: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/

Download the genomic GTF file. In the moment of writing, it is this one:

https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz

This file also contains relevant ranges for MANE transcripts and their exons.