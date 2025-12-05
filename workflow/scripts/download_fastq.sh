#!/bin/bash

# Directory to store FASTQ files
# TODO: Read the output directory from command line argument
OUTDIR="/home/nobackup/molinerislab/mosquito/fastq_jenkins"

# Create directory if it doesn't exist
mkdir -p "$OUTDIR"

# List of SRA accessions to download
#TODO: Read the SRA accessions from a file, so that you don't have to hardcode them and make the script executable every time
SRA_LIST=(
    ERR440788
    ERR440789
    ERR440790
    ERR440791
    ERR440792
    ERR440793
)

for SRA in "${SRA_LIST[@]}"; do
    echo "Downloading $SRA ..."
    fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-3 --clip --outdir "$OUTDIR" "$SRA" #split-3: split into files for paired-end reads | gzip compression | other flags for quality control
done

echo "Download complete."

#TODO: Add verification of downloaded files and error handling
