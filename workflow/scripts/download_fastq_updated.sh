#!/bin/bash
set -euo pipefail

#---------------------------------------------#
# Usage:
#   ./download_fastq.sh /path/to/config.yaml
#---------------------------------------------#

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/config.yaml"
    exit 1
fi

CONFIG="$1"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file '$CONFIG' not found."
    exit 1
fi

# --- Check if Python can import PyYAML ---
if ! python3 -c "import yaml" &>/dev/null; then
    echo "ERROR: Python module 'yaml' not found."
    echo "Please install PyYAML, e.g.:"
    echo "  conda install -c conda-forge pyyaml"
    echo "or"
    echo "  pip install --user pyyaml"
    exit 1
fi

# --- Read OUTDIR and SRA_LIST from config.yaml ---
eval "$(python3 - <<EOF
import yaml, shlex
config = yaml.safe_load(open("$CONFIG"))

outdir = config['FASTQ']['OUTDIR']
sras = config['FASTQ']['SRA']

# Print Bash assignments
print(f"OUTDIR={shlex.quote(outdir)}")
print(f"SRA_LIST=({ ' '.join(shlex.quote(s) for s in sras) })")
EOF
)"

# --- Validate ---
if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: FASTQ.OUTDIR not defined in $CONFIG"
    exit 1
fi

if [[ ${#SRA_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No SRA accessions found in FASTQ.SRA"
    exit 1
fi

# --- Create output directory ---
mkdir -p "$OUTDIR"

# --- Download loop with validation ---
for SRA in "${SRA_LIST[@]}"; do
    echo "Downloading $SRA ..."
    fastq-dump --skip-technical --gzip --readids --read-filter pass --dumpbase --split-3 --clip --outdir "$OUTDIR" "$SRA"

    # --- Validation: check that at least one FASTQ file exists ---
    if ! ls "$OUTDIR" | grep -q "$SRA"; then
        echo "ERROR: No FASTQ files found for $SRA — download may have failed."
        exit 1
    fi

    echo "Finished downloading: $SRA"
    echo "--------------------------------------------"
done

echo "All downloads completed successfully."
