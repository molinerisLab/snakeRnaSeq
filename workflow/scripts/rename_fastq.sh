#!/bin/bash

# Check for required argument
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 sample_map.txt"
    exit 1
fi

mapfile="$1"

# Ensure mapping file exists
if [[ ! -f "$mapfile" ]]; then
    echo "Error: Mapping file '$mapfile' not found"
    exit 1
fi

# Path to the directory containing the FASTQ files
FASTQ_DIR="/home/nobackup/molinerislab/mosquito/fastq_saizonou"

# Read mapping and rename files
while read -r sra sample || [[ -n "$sra" ]]; do

    # Skip empty lines
    # [[ -z "$sample" || -z "$sra" ]] && continue
    [[ -z "$sra" || -z "$sample" ]] && continue

    R1="${FASTQ_DIR}/${sra}_pass_1.fastq.gz"
    R2="${FASTQ_DIR}/${sra}_pass_2.fastq.gz"

    # Check if files exist and rename
    if [[ -f "$R1" ]]; then
        echo "Renaming $R1 → ${FASTQ_DIR}/${sample}_R1.fastq.gz"
        mv "$R1" "${FASTQ_DIR}/${sample}_R1.fastq.gz"
    else
        echo "Warning: $R1 not found"
    fi

    if [[ -f "$R2" ]]; then
        echo "Renaming $R2 → ${FASTQ_DIR}/${sample}_R2.fastq.gz"
        mv "$R2" "${FASTQ_DIR}/${sample}_R2.fastq.gz"
    else
        echo "Warning: $R2 not found"
    fi

done < "$mapfile"
