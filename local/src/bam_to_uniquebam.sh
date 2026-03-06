#!/bin/bash
THREADS=8
while read s; do
    IN="$SPL/${s}.bam"
    OUT="$UNIQ/${s}.unique.bam"
    # se già fatto, salta
    if [[ -s "$OUT" && -s "${OUT}.bai" ]]; then
        echo "[SKIP] $s"
        continue
    fi
    echo "=== UNIQUE $s ==="
    samtools view -@ "$THREADS" -b -q 255 "$IN" | samtools sort -@ "$THREADS" -o "$OUT" -
    samtools index -@ "$THREADS" "$OUT"
done < /tmp/spladder_345.samples