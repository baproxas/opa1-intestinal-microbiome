#!/bin/bash
# ============================================================================== 
# 02_kraken2_array_core_nt.slurm
# Classify all samples with Kraken2 using core_nt DB (two-pass host decontamination)
#
# IMPORTANT: core_nt index is ~316 GB. Standard Puma nodes have only 239 GB RAM.
# You MUST use high-memory nodes (3008 GB RAM) via --constraint=hi_mem.
# There are only 5 hi_mem nodes on Puma — jobs may queue.
#
# Two-pass strategy:
#   Pass 1: k2_host_db  (~20 GB RAM, ~2-5 min)  -> extract non-host reads
#   Pass 2: core_nt DB (~316 GB RAM, ~45-90 min) -> final taxonomy report
# ============================================================================== 

#SBATCH --job-name=kraken2_core_nt_classify
#SBATCH --account=vkv
#SBATCH --partition=standard
#SBATCH --constraint=hi_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=08:00:00
#SBATCH --array=0-34
#SBATCH --output=/groups/vkv/baproxas/opa1_drp/logs/kraken2_core_nt_%A_%a.out
#SBATCH --error=/groups/vkv/baproxas/opa1_drp/logs/kraken2_core_nt_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=baproxas@arizona.edu

set -euo pipefail

PROJECT_DIR="/groups/vkv/baproxas/opa1_drp"
CONDA_ENV="/groups/vkv/baproxas/conda_envs/metagenomics"
DB_PATH="/xdisk/vkv/baproxas/opa1_drp/data/db/k2_core_nt_20251015"
HOST_DB_PATH="/xdisk/vkv/baproxas/opa1_drp/data/db/k2_host_db"
FASTQ_DIR="$PROJECT_DIR/data/fastq"
SAMPLE_LIST="$FASTQ_DIR/sample_list.txt"
REPORT_DIR="$PROJECT_DIR/temp/kraken_reports_core_nt"
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Activate conda environment
source /opt/lmod/lmod/init/bash 2>/dev/null || true
module --ignore_cache load micromamba/2.0.2-2 2>/dev/null || true
export PATH="$CONDA_ENV/bin:$PATH"

mkdir -p "$REPORT_DIR" "$PROJECT_DIR/logs"

# Get sample name for this array task
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

FASTQ="$FASTQ_DIR/${SAMPLE}.fastq"
REPORT="$REPORT_DIR/${SAMPLE}.k2report"
HOST_REPORT="$REPORT_DIR/${SAMPLE}_host_filter.k2report"
FILTERED_FASTQ="$REPORT_DIR/${SAMPLE}_nonhost.fastq"

echo "============================================="
echo "  Kraken2 Classification (core_nt)"
echo "  Array task: $SLURM_ARRAY_TASK_ID"
echo "  Sample:     $SAMPLE"
echo "  FASTQ:      $FASTQ"
echo "  core_nt DB: $DB_PATH"
echo "  Host DB:    $HOST_DB_PATH"
echo "  Report:     $REPORT"
echo "  Threads:    $THREADS"
echo "  Node:       $(hostname)"
echo "  Date:       $(date)"
echo "============================================="

# Skip if already classified
if [ -f "$REPORT" ]; then
    echo "[SKIP] Report already exists: $REPORT"
    exit 0
fi

if [ ! -f "$FASTQ" ]; then
    echo "ERROR: FASTQ file not found: $FASTQ"
    exit 1
fi

# Pass 1: Host decontamination
kraken2 \
    --db "$HOST_DB_PATH" \
    --threads "$THREADS" \
    --unclassified-out "$FILTERED_FASTQ" \
    --report "$HOST_REPORT" \
    --output /dev/null \
    "$FASTQ"

echo "  Pass 1 complete. Non-host reads: $FILTERED_FASTQ"

# Pass 2: core_nt classification
kraken2 \
    --db "$DB_PATH" \
    --threads "$THREADS" \
    --output /dev/null \
    --report "$REPORT" \
    --confidence 0.01 \
    "$FILTERED_FASTQ"

echo "[Pass 2] Done. Report: $REPORT"

# Clean up temp file
rm -f "$FILTERED_FASTQ"
echo "  Temp non-host FASTQ removed."
echo "Completed: $(date)"
