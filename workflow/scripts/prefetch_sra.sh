#!/bin/bash
set -euo pipefail

# Usage and argument parsing
usage() {
    echo "Usage: $0 -o OUTDIR -s SRA_LIST_FILE"
    echo "  -o OUTDIR          Output directory for FASTQ files"
    echo "  -s SRA_LIST_FILE   File containing one SRA accession per line"
    exit 1
}

# Default values
OUTDIR=""
SRA_FILE=""

# Parse command-line options
while getopts "o:s:" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        s) SRA_FILE="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$OUTDIR" || -z "$SRA_FILE" ]]; then
    usage
fi

# Validate SRA file
if [[ ! -f "$SRA_FILE" ]]; then
    echo "ERROR: SRA list file '$SRA_FILE' not found."
    exit 1
fi

if [[ ! -s "$SRA_FILE" ]]; then
    echo "ERROR: SRA list file '$SRA_FILE' is empty."
    exit 1
fi

# Create output directory
mkdir -p "$OUTDIR"

# Log file
LOGFILE="$OUTDIR/download_log.txt"
touch "$LOGFILE"

# Download loop
echo "Starting downloads..."
echo "Log file: $LOGFILE"

while read -r SRA; do
    # Skip empty or commented lines
    [[ -z "$SRA" || "$SRA" =~ ^# ]] && continue

    echo "Downloading $SRA ..."
    if fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-3 --clip --outdir "$OUTDIR" "$SRA" >>"$LOGFILE" 2>&1; then

        # Verify files were created
        if ls "$OUTDIR/${SRA}"*.fastq.gz >/dev/null 2>&1; then
            echo "Successfully downloaded $SRA"
        else
            echo "ERROR: No FASTQ files found for $SRA" | tee -a "$LOGFILE"
        fi
    else
        echo "ERROR: fastq-dump failed for $SRA" | tee -a "$LOGFILE"
    fi
done < "$SRA_FILE"

echo "Download complete."
echo "Check the log file for details: $LOGFILE"