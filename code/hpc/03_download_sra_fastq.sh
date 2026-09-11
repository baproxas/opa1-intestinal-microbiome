#!/bin/bash
# Download FASTQ files for BioProject PRJNA1523974 once it is public in NCBI SRA.
#
# Prerequisites:
#   - SRA Toolkit with prefetch and fasterq-dump on PATH
#   - NCBI Entrez Direct with esearch and efetch on PATH
#
# This script will fail until PRJNA1523974 is released by NCBI. It is safe to
# rerun: prefetch resumes downloads and fasterq-dump skips complete outputs.
#
# Usage:
#   ./03_download_sra_fastq.sh [OUTPUT_DIR]
#   ./03_download_sra_fastq.sh --output-dir OUTPUT_DIR [--cleanup]
#
# Interactive example:
#   ./03_download_sra_fastq.sh data/raw
#
# SLURM example:
#   sbatch --cpus-per-task=8 --mem=64G --time=72:00:00 \
#       --wrap='bash code/hpc/03_download_sra_fastq.sh data/raw'
#
# Optional SLURM header for users who want to save this file as a job script:
#SBATCH --job-name=sra_fastq_download
#SBATCH --account=vkv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

BIOPROJECT="PRJNA1523974"
OUTPUT_DIR="$(pwd)"
THREADS=8
CLEANUP=false
PREFETCH_DIR=""

usage() {
	sed -n '1,28p' "$0"
}

while [[ $# -gt 0 ]]; do
	case "$1" in
		--cleanup)
			CLEANUP=true
			shift
			;;
		--output-dir)
			[[ $# -ge 2 ]] || { echo "ERROR: --output-dir requires a path." >&2; exit 2; }
			OUTPUT_DIR="$2"
			shift 2
			;;
		--threads)
			[[ $# -ge 2 ]] || { echo "ERROR: --threads requires an integer." >&2; exit 2; }
			THREADS="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		--*)
			echo "ERROR: Unknown option: $1" >&2
			usage >&2
			exit 2
			;;
		*)
			if [[ "$OUTPUT_DIR" != "$(pwd)" ]]; then
				echo "ERROR: Output directory supplied more than once." >&2
				exit 2
			fi
			OUTPUT_DIR="$1"
			shift
			;;
	esac
done

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
	echo "ERROR: THREADS must be a positive integer; received '$THREADS'." >&2
	exit 2
fi

if ! mkdir -p "$OUTPUT_DIR"; then
	echo "ERROR: Could not create output directory: $OUTPUT_DIR" >&2
	exit 2
fi
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

PREFETCH_DIR="$OUTPUT_DIR/sra_cache"
SAMPLE_LIST="$OUTPUT_DIR/sample_list.txt"
ACCESSION_LIST="$OUTPUT_DIR/srr_accessions.txt"
mkdir -p "$PREFETCH_DIR"

echo "============================================="
echo "  NCBI SRA FASTQ Download"
echo "  BioProject: $BIOPROJECT"
echo "  Output: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo "============================================="

echo "[Step 1] Checking SRA Toolkit and Entrez Direct availability..."
if ! command -v prefetch >/dev/null 2>&1 || ! command -v fasterq-dump >/dev/null 2>&1; then
	echo "ERROR: SRA Toolkit is required (prefetch and fasterq-dump)." >&2
	exit 1
fi
if ! command -v esearch >/dev/null 2>&1 || ! command -v efetch >/dev/null 2>&1; then
	echo "ERROR: Entrez Direct is required (esearch and efetch)." >&2
	exit 1
fi
echo "  prefetch: $(command -v prefetch)"
echo "  fasterq-dump: $(command -v fasterq-dump)"

echo "[Step 2] Checking whether $BIOPROJECT is public..."
BIOPROJECT_QUERY="$(mktemp)"
trap 'rm -f "$BIOPROJECT_QUERY"' EXIT
if ! esearch -db bioproject -query "${BIOPROJECT}[Accession]" >"$BIOPROJECT_QUERY" 2>/dev/null; then
	echo "WARNING: NCBI could not query $BIOPROJECT. It may not yet be public." >&2
	exit 1
fi
if ! grep -q '<Id>' "$BIOPROJECT_QUERY"; then
	echo "WARNING: $BIOPROJECT is not yet public in NCBI." >&2
	echo "The download will fail until NCBI releases the BioProject." >&2
	exit 1
fi
echo "  BioProject is visible in NCBI."

echo "[Step 3] Retrieving SRR accessions from $BIOPROJECT..."
if ! esearch -db sra -query "$BIOPROJECT" \
	| efetch -format runinfo \
	| awk -F',' 'NR > 1 && $1 ~ /^SRR[0-9]+$/ { print $1 }' \
	| sort -u >"$ACCESSION_LIST"; then
	echo "ERROR: Failed to retrieve SRA run information." >&2
	exit 1
fi
if [[ ! -s "$ACCESSION_LIST" ]]; then
	echo "WARNING: No SRR accessions were found for $BIOPROJECT." >&2
	echo "The BioProject may be public but its run records may not be released yet." >&2
	exit 1
fi
echo "  Found $(wc -l <"$ACCESSION_LIST") SRR accessions."

echo "[Step 4] Prefetching SRA accessions (resumable, max 50 GB each)..."
while IFS= read -r accession; do
	echo "  Prefetching $accession..."
	prefetch --max-size 50g --output-directory "$PREFETCH_DIR" "$accession"
done <"$ACCESSION_LIST"

echo "[Step 5] Converting prefetched accessions to FASTQ..."
while IFS= read -r accession; do
	output_pattern="$OUTPUT_DIR/${accession}_1.fastq"
	single_output="$OUTPUT_DIR/${accession}.fastq"
	if [[ -s "$output_pattern" || -s "$single_output" ]]; then
		echo "  FASTQ output already exists for $accession; skipping."
		continue
	fi
	echo "  Converting $accession..."
	fasterq-dump "$PREFETCH_DIR/$accession" \
		--outdir "$OUTPUT_DIR" \
		--threads "$THREADS" \
		--split-files
done <"$ACCESSION_LIST"

echo "[Step 6] Validating FASTQ outputs..."
missing=0
while IFS= read -r accession; do
	found=false
	for fastq in "$OUTPUT_DIR/${accession}.fastq" \
				"$OUTPUT_DIR/${accession}_1.fastq" \
				"$OUTPUT_DIR/${accession}_2.fastq"; do
		if [[ -e "$fastq" ]]; then
			if [[ ! -s "$fastq" ]]; then
				echo "ERROR: FASTQ file is empty: $fastq" >&2
				missing=1
			else
				found=true
			fi
		fi
	done
	if [[ "$found" != true ]]; then
		echo "ERROR: No non-empty FASTQ output found for $accession." >&2
		missing=1
	fi
done <"$ACCESSION_LIST"
if [[ "$missing" -ne 0 ]]; then
	exit 1
fi

if [[ "$CLEANUP" == true ]]; then
	echo "[Step 7] Cleaning up SRA cache files..."
	while IFS= read -r accession; do
		rm -rf "$PREFETCH_DIR/$accession"
	done <"$ACCESSION_LIST"
	rmdir "$PREFETCH_DIR" 2>/dev/null || true
else
	echo "[Step 7] Keeping SRA cache. Re-run with --cleanup to remove it."
fi

echo "[Step 8] Generating sample_list.txt..."
find "$OUTPUT_DIR" -maxdepth 1 -type f \
	\( -name '*.fastq' -o -name '*.fastq.gz' \) \
	-printf '%f\n' \
	| sed -E 's/(_1|_2)?\.fastq(\.gz)?$//' \
	| sort -u >"$SAMPLE_LIST"

echo ""
echo "Download summary"
echo "================"
echo "SRR accessions: $(wc -l <"$ACCESSION_LIST")"
echo "FASTQ files:    $(find "$OUTPUT_DIR" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | wc -l)"
echo "Sample list:    $SAMPLE_LIST"
echo "FASTQ directory: $OUTPUT_DIR"
echo "Done."
