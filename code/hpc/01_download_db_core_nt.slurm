#!/bin/bash
# ============================================================================== 
# 01_download_db_core_nt.slurm
# Download Kraken2 core_nt database to /xdisk/vkv/baproxas/opa1_drp/data/db
# ============================================================================== 

#SBATCH --job-name=kraken2_db_core_nt
#SBATCH --account=vkv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/groups/vkv/baproxas/opa1_drp/logs/%x_%j.out
#SBATCH --error=/groups/vkv/baproxas/opa1_drp/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=baproxas@arizona.edu

set -euo pipefail

PROJECT_DIR="/groups/vkv/baproxas/opa1_drp"
DB_DIR="/xdisk/vkv/baproxas/opa1_drp/data/db"
DB_NAME="k2_core_nt_20251015"
DB_PATH="$DB_DIR/$DB_NAME"
TAR_PATH="$DB_DIR/${DB_NAME}.tar.gz"
DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz"

# Reject tiny/partial files (the real archive is many GB, so 1 GB is a safe floor)
MIN_BYTES=1073741824

mkdir -p "$DB_DIR" "$DB_PATH" "$PROJECT_DIR/logs"

echo "============================================="
echo "  Kraken2 core_nt Database Download"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  CPUs: ${SLURM_CPUS_PER_TASK:-8}"
echo "  DB path: $DB_PATH"
echo "============================================="

echo "[Step 0] Checking archive URL..."
if ! wget --spider -q "$DB_URL"; then
    echo "ERROR: URL not reachable: $DB_URL"
    exit 4
fi

echo "[Step 1] Downloading core_nt database..."
if [ -f "$TAR_PATH" ]; then
    FILE_SIZE=$(stat -c%s "$TAR_PATH" 2>/dev/null || echo 0)
    if [ "$FILE_SIZE" -lt "$MIN_BYTES" ]; then
        echo "  Existing archive is too small (${FILE_SIZE} bytes). Removing and re-downloading."
        rm -f "$TAR_PATH"
    elif ! gzip -t "$TAR_PATH" 2>/dev/null; then
        echo "  Existing archive failed gzip integrity test. Removing and re-downloading."
        rm -f "$TAR_PATH"
    else
        echo "  Existing archive looks valid ($(du -sh "$TAR_PATH" | cut -f1)); skipping download."
    fi
fi

if [ ! -f "$TAR_PATH" ]; then
    wget --continue --progress=dot:giga -O "$TAR_PATH" "$DB_URL"
    FILE_SIZE=$(stat -c%s "$TAR_PATH" 2>/dev/null || echo 0)
    if [ "$FILE_SIZE" -lt "$MIN_BYTES" ]; then
        echo "ERROR: Downloaded file is unexpectedly small (${FILE_SIZE} bytes)."
        exit 2
    fi
    if ! gzip -t "$TAR_PATH" 2>/dev/null; then
        echo "ERROR: Downloaded file failed gzip integrity test."
        exit 2
    fi
    echo "  Download complete. $(du -sh "$TAR_PATH" | cut -f1) on disk."
fi

echo "[Step 2] Extracting archive..."
if [ -f "$DB_PATH/hash.k2d" ]; then
    echo "  Skipping extraction - $DB_PATH/hash.k2d already exists."
else
    tar -xzf "$TAR_PATH" -C "$DB_PATH"
    for f in hash.k2d opts.k2d taxo.k2d; do
        if [ ! -f "$DB_PATH/$f" ]; then
            echo "ERROR: Missing expected file after extraction: $DB_PATH/$f"
            exit 3
        fi
    done
    echo "  Extraction complete. Core files verified."
fi

echo "[Done] core_nt DB ready at $DB_PATH"
echo "  Size: $(du -sh "$DB_PATH" | cut -f1)"
echo "Completed: $(date)"
