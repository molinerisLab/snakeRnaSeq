STAR_DIR="/mnt/oncog/software/snakeRnaSeq/dataset/v1/Results/pass2"
LINK_DIR="bam_pass2_renamed"

for SAMPLE_DIR in "$STAR_DIR"/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    BAM="$SAMPLE_DIR/Aligned.sortedByCoord.out.bam"
    BAI="$BAM.bai"

    if [[ -f "$BAM" && -f "$BAI" ]]; then
        ln -sf "$BAM" "$LINK_DIR/${SAMPLE}.bam"
        ln -sf "$BAI" "$LINK_DIR/${SAMPLE}.bam.bai"
        echo "✅ Link creati per $SAMPLE"
    else
        echo "⚠️ BAM o BAI mancanti per $SAMPLE"
    fi
done
