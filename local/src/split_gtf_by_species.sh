#!/usr/bin/env bash

set -euo pipefail

GTF="$1"

DIR=$(dirname "$GTF")
BASE=$(basename "$GTF" .gtf)

HUMAN_OUT="${DIR}/splitted_human.gtf"
MOUSE_OUT="${DIR}/splitted_mouse.gtf"

echo "Input GTF: $GTF"
echo "Human output: $HUMAN_OUT"
echo "Mouse output: $MOUSE_OUT"

awk '
/^#/ {print > human; print > mouse; next}
$1 ~ /^H/ {print > human}
$1 ~ /^M/ {print > mouse}
' human="$HUMAN_OUT" mouse="$MOUSE_OUT" "$GTF"

echo "Done."
echo "Line counts:"
wc -l "$HUMAN_OUT" "$MOUSE_OUT"
