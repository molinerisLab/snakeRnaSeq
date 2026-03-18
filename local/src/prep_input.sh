#!/bin/bash

OUT="../../dataset/v1/stringtie/prepDE_input.txt"

> "$OUT"

for gtf in ../../dataset/v1/stringtie/*/quant_merged.gtf; do
    sample_name=$(basename "$(dirname "$gtf")")
    real_gtf=$(realpath "$gtf")
    echo -e "${sample_name}\t${real_gtf}" >> "$OUT"
done

echo "Creato $OUT"
