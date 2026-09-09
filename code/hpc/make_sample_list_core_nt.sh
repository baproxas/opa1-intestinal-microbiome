#!/bin/bash

PROJECT_DIR="/groups/vkv/baproxas/opa1_drp"
FASTQ_DIR="$PROJECT_DIR/data/fastq"
SAMPLE_LIST="$FASTQ_DIR/sample_list.txt"

echo "Scanning for FASTQ files in: $FASTQ_DIR"

find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq" \
    | xargs -I{} basename {} .fastq \
    | sort \
    > "$SAMPLE_LIST"

N=$(wc -l < "$SAMPLE_LIST")
echo "Found $N samples. Written to: $SAMPLE_LIST"
echo ""
echo "Sample list:"
cat "$SAMPLE_LIST"
echo ""
echo "Update 02_kraken2_array_core_nt.slurm --array=0-$((N-1)) if needed"
