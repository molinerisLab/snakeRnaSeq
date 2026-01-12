#!/bin/bash

# Usage 
usage() {
    echo "Usage: $0 -f SAMPLE_MAP_FILE -d FASTQ_DIR"
    echo "  -f SAMPLE_MAP_FILE   Mapping file: <SRA_ID> <SampleName>"
    echo "  -d FASTQ_DIR         Directory containing FASTQ files"
    exit 1
}

# Argument parsing
MAPFILE=""
FASTQ_DIR=""

while getopts ":f:d:" opt; do
    case "${opt}" in
        f) MAPFILE="${OPTARG}" ;;
        d) FASTQ_DIR="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "${MAPFILE}" || -z "${FASTQ_DIR}" ]]; then
    usage
fi

# Validate inputs
if [[ ! -f "${MAPFILE}" ]]; then
    echo "ERROR: Mapping file not found: ${MAPFILE}" >&2
    exit 1
fi

if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}" >&2
    exit 1
fi

# Read mapping and rename files
while read -r sra sample || [[ -n "$sra" ]]; do

    # Skip empty lines
    [[ -z "$sra" || -z "$sample" ]] && continue

    R1="${FASTQ_DIR}/${sra}_1.fastq.gz"
    R2="${FASTQ_DIR}/${sra}_2.fastq.gz"

    # Rename Read 1
    if [[ -f "$R1" ]]; then
        echo "Renaming $R1 → ${FASTQ_DIR}/${sample}_R1.fastq.gz"
        mv "$R1" "${FASTQ_DIR}/${sample}_R1.fastq.gz"
    else
        echo "Warning: $R1 not found"
    fi

    # Rename Read 2
    if [[ -f "$R2" ]]; then
        echo "Renaming $R2 → ${FASTQ_DIR}/${sample}_R2.fastq.gz"
        mv "$R2" "${FASTQ_DIR}/${sample}_R2.fastq.gz"
    else
        echo "Warning: $R2 not found"
    fi

done < "${MAPFILE}"
