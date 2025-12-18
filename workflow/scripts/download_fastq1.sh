#!/bin/bash
set -Eeuo pipefail


# Usage --------------------------------------
usage() {
    echo "Usage: $0 -o OUTDIR -s SRA_LIST_FILE"
    echo "  -o OUTDIR          Output directory for FASTQ files"
    echo "  -s SRA_LIST_FILE   File containing one SRA accession per line"
    exit 1
}

# Argument parsing ---------------------------
OUTDIR=""
SRA_FILE=""

while getopts "o:s:" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        s) SRA_FILE="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$OUTDIR" || -z "$SRA_FILE" ]] && usage
[[ ! -f "$SRA_FILE" ]] && { echo "ERROR: '$SRA_FILE' not found"; exit 1; }
[[ ! -s "$SRA_FILE" ]] && { echo "ERROR: '$SRA_FILE' is empty"; exit 1; }

# TMPDIR handling ---------------------------
BASE_TMPDIR="${TMPDIR:-/tmp}"

WORKDIR="$(mktemp -d "$BASE_TMPDIR/fastq_dump.XXXXXX")"
ERRLOG="$WORKDIR/fastq_dump.err.log"

cleanup() {
    rm -rf "$WORKDIR"
}
trap cleanup EXIT INT TERM

echo "Using temporary directory: $WORKDIR"

# Output directory ---------------------------
mkdir -p "$OUTDIR"

# Download loop ------------------------------- 
echo "Starting downloads..."

while read -r SRA; do
    [[ -z "$SRA" || "$SRA" =~ ^# ]] && continue

    echo "Downloading $SRA ..."

    if fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-3 --clip --outdir "$WORKDIR" "$SRA" 2>>"$ERRLOG"
    then
        if ls "$WORKDIR/${SRA}"*.fastq.gz >/dev/null 2>&1; then
            echo "✓ $SRA downloaded successfully"

            mv "$WORKDIR/${SRA}"*.fastq.gz "$OUTDIR/"
        else
            echo "ERROR: No FASTQ files produced for $SRA"
        fi
    else
        echo "ERROR: fastq-dump failed for $SRA (see error log)"
    fi

done < "$SRA_FILE"

# Final status --------------------------
if [[ -s "$ERRLOG" ]]; then
    echo
    echo "Some errors occurred. Error log:"
    echo "$ERRLOG"
else
    echo
    echo "All downloads completed successfully."
fi
