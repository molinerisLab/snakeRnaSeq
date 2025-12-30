#!/bin/bash
set -euo pipefail

# Usage -----------------------------------------------
usage() {
    echo "Usage: $0 -o OUTDIR -s SRA_LIST_FILE"
    echo "  -o OUTDIR          Output directory for FASTQ files"
    echo "  -s SRA_LIST_FILE   File containing one SRA accession per line"
    exit 1
}

# Argument parsing ------------------------------------
OUTDIR=""
SRA_LIST=""

while getopts ":o:s:" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        s) SRA_LIST="${OPTARG}" ;;
        *) usage ;;
    esac
done

if [[ -z "${OUTDIR}" || -z "${SRA_LIST}" ]]; then
    usage
fi

if [[ ! -f "${SRA_LIST}" ]]; then
    echo "ERROR: SRA list file not found: ${SRA_LIST}" >&2
    exit 1
fi

# TMPDIR handling --------------------------------------
TMPBASE="${TMPDIR:-/mnt/nvme_raid0/tmp}"
WORKTMP="${TMPBASE}/sra_work"
mkdir -p "${WORKTMP}"

# Setup output + logging ------------------------------
mkdir -p "${OUTDIR}"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${LOGDIR}"

LOGFILE="${LOGDIR}/download_$(date +%Y%m%d_%H%M%S).log"

{
    echo "[$(date)] Starting SRA pipeline"
    echo "OUTDIR: ${OUTDIR}"
    echo "SRA LIST: ${SRA_LIST}"
    echo "WORKTMP: ${WORKTMP}"
    echo "fasterq-dump: $(command -v fasterq-dump)"
    fasterq-dump --version
} | tee -a "${LOGFILE}"

# Prefetch SRA files ---------------------------------
echo "[$(date)] Checking SRA files" | tee -a "${LOGFILE}"

while read -r ACC; do
    [[ -z "${ACC}" ]] && continue
    ACC="${ACC%%[$'\r\t ']*}"

    SRA_PATH="${WORKTMP}/${ACC}/${ACC}.sra"

    if [[ -f "${SRA_PATH}" ]]; then
        echo "[$(date)] ${ACC}: SRA exists, skipping prefetch" | tee -a "${LOGFILE}"
    else
        echo "[$(date)] ${ACC}: downloading SRA" | tee -a "${LOGFILE}"
        prefetch "${ACC}" --output-directory "${WORKTMP}" --max-size 200G --progress 2>&1 | tee -a "${LOGFILE}"
    fi
done < "${SRA_LIST}"

# Convert SRA → FASTQ --------------------------------
while read -r ACC; do
    [[ -z "${ACC}" ]] && continue
    ACC="${ACC%%[$'\r\t ']*}"

    SRA_PATH="${WORKTMP}/${ACC}/${ACC}.sra"
    FASTQ1="${OUTDIR}/${ACC}_1.fastq"
    FASTQ2="${OUTDIR}/${ACC}_2.fastq"
    FASTQSE="${OUTDIR}/${ACC}.fastq"

    if [[ -f "${FASTQ1}" || -f "${FASTQSE}" ]]; then
        echo "[$(date)] ${ACC}: FASTQ exists, skipping conversion" | tee -a "${LOGFILE}"
        continue
    fi

    echo "[$(date)] ${ACC}: converting to FASTQ" | tee -a "${LOGFILE}"

    fasterq-dump "${SRA_PATH}" --outdir "${OUTDIR}" -t "${WORKTMP}" --threads 6 --split-files --skip-technical 2>&1 | tee -a "${LOGFILE}"

done < "${SRA_LIST}"

# Compress FASTQ -------------------------------------
echo "[$(date)] Compressing FASTQ files" | tee -a "${LOGFILE}"

find "${OUTDIR}" -maxdepth 1 -name "*.fastq" ! -name "*.gz" \
    -print -exec pigz -p 6 {} \; | tee -a "${LOGFILE}"

echo "[$(date)] Pipeline completed successfully" | tee -a "${LOGFILE}"
