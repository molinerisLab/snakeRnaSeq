#!/bin/bash

# Set the folder containing BAM files
BAM_FOLDER="${1:-.}"  # Use first argument or current directory

OUTPUT_FILE="out_spladder/alignment.txt" 

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Find all .bam files and write their full paths to the output file
find "$BAM_FOLDER" -name "*.bam"  | while read -r bam_file; do
    realpath "$bam_file" >> "$OUTPUT_FILE"
done
chmod g+rw "$OUTPUT_FILE"
echo "BAM file list created: $OUTPUT_FILE"
echo "Total BAM files found: $(wc -l < "$OUTPUT_FILE")"