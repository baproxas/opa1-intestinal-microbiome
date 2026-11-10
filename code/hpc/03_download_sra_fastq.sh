#!/bin/bash
# Download FASTQ files for BioProject PRJNA1523974 from NCBI SRA and name them
# by SRA SampleName (OPA01.fastq ... OPA18.fastq) to match
# config/sample_metadata and 02_kraken2_array_core_nt.slurm.
#
# Prerequisites:
#   - SRA Toolkit with prefetch and fasterq-dump on PATH
#   - NCBI Entrez Direct with esearch and efetch on PATH
#
# Safe to rerun: prefetch resumes downloads and finished samples are skipped.
# Set BIOPROJECT=PRJNAxxxx in the environment to test with another project.
#
# Usage:
#   ./03_download_sra_fastq.sh [OUTPUT_DIR]
#   ./03_download_sra_fastq.sh --output-dir OUTPUT_DIR [--threads N] [--cleanup]
#
# Interactive example:
#   ./03_download_sra_fastq.sh data/fastq
#
# SLURM example (create logs/ first):
#   sbatch --cpus-per-task=8 --mem=16G --time=12:00:00 \
#       --wrap='bash code/hpc/03_download_sra_fastq.sh data/fastq'
#
# Optional SLURM header for users who want to save this file as a job script:
#SBATCH --job-name=sra_fastq_download
#SBATCH --account=vkv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

BIOPROJECT="${BIOPROJECT:-PRJNA1523974}"
OUTPUT_DIR="$(pwd)"
OUTPUT_DIR_SET=false
THREADS=8
CLEANUP=false

usage() {
	sed -n '1,35p' "$0"
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
			OUTPUT_DIR_SET=true
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
			if [[ "$OUTPUT_DIR_SET" == true ]]; then
				echo "ERROR: Output directory supplied more than once." >&2
				exit 2
			fi
			OUTPUT_DIR="$1"
			OUTPUT_DIR_SET=true
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
TEMP_DIR="$OUTPUT_DIR/fasterq_tmp"
RUNINFO="$OUTPUT_DIR/sra_runinfo.csv"
RUN_MAP="$OUTPUT_DIR/srr_to_sample.tsv"
ACCESSION_LIST="$OUTPUT_DIR/srr_accessions.txt"
SAMPLE_LIST="$OUTPUT_DIR/sample_list.txt"
mkdir -p "$PREFETCH_DIR" "$TEMP_DIR"

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

# Gate on the SRA runs, not the BioProject record: the record in the
# bioproject database can lag days behind run release. esearch prints an
# ENTREZ_DIRECT summary with <Count>; it never contains <Id> elements.
echo "[Step 2] Checking whether $BIOPROJECT runs are public in SRA..."
SRA_QUERY="${BIOPROJECT}[BioProject]"
if ! sra_count="$(esearch -db sra -query "$SRA_QUERY" \
	| awk -F'[<>]' '/<Count>/ && !seen { print $3; seen = 1 }')"; then
	echo "ERROR: NCBI esearch failed for $SRA_QUERY (network or NCBI outage?)." >&2
	exit 1
fi
if ! [[ "$sra_count" =~ ^[0-9]+$ ]] || [[ "$sra_count" -eq 0 ]]; then
	echo "WARNING: No public SRA records for $BIOPROJECT yet (count='${sra_count}')." >&2
	exit 1
fi
echo "  $sra_count public SRA records found."

echo "[Step 3] Retrieving run info and SRR -> SampleName map..."
if ! esearch -db sra -query "$SRA_QUERY" | efetch -format runinfo >"$RUNINFO"; then
	echo "ERROR: Failed to retrieve SRA run information." >&2
	exit 1
fi
# Locate columns by header name; efetch may repeat the header per batch.
awk -F',' -v OFS='\t' '
	$1 == "Run" {
		r = s = 0
		for (i = 1; i <= NF; i++) { if ($i == "Run") r = i; if ($i == "SampleName") s = i }
		next
	}
	r && s && $r ~ /^[SED]RR[0-9]+$/ { print $r, $s }
' "$RUNINFO" | sort -u >"$RUN_MAP"
if [[ ! -s "$RUN_MAP" ]]; then
	echo "ERROR: No runs with Run/SampleName parsed from $RUNINFO." >&2
	exit 1
fi
while IFS=$'\t' read -r accession sample; do
	if ! [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]]; then
		echo "ERROR: Unsafe or empty SampleName '$sample' for $accession." >&2
		exit 1
	fi
done <"$RUN_MAP"
dup_samples="$(cut -f2 "$RUN_MAP" | sort | uniq -d)"
if [[ -n "$dup_samples" ]]; then
	echo "ERROR: Multiple runs share a SampleName; outputs would collide:" >&2
	echo "$dup_samples" >&2
	exit 1
fi
cut -f1 "$RUN_MAP" >"$ACCESSION_LIST"
echo "  Found $(wc -l <"$RUN_MAP") runs. Map written to $RUN_MAP"

echo "[Step 4] Prefetching SRA accessions (resumable, max 50 GB each)..."
while IFS=$'\t' read -r accession sample; do
	echo "  Prefetching $accession ($sample)..."
	prefetch --max-size 50g --output-directory "$PREFETCH_DIR" "$accession" </dev/null
done <"$RUN_MAP"

echo "[Step 5] Converting to FASTQ and renaming SRR -> SampleName..."
while IFS=$'\t' read -r accession sample; do
	if [[ -s "$OUTPUT_DIR/${sample}.fastq" || -s "$OUTPUT_DIR/${sample}_1.fastq" ]]; then
		echo "  FASTQ output already exists for $sample; skipping."
		continue
	fi
	echo "  Converting $accession -> $sample..."
	# --force overwrites partial SRR-named output from an interrupted run.
	fasterq-dump "$PREFETCH_DIR/$accession" \
		--outdir "$OUTPUT_DIR" \
		--temp "$TEMP_DIR" \
		--threads "$THREADS" \
		--split-files \
		--force </dev/null
	for suffix in "" _1 _2; do
		src="$OUTPUT_DIR/${accession}${suffix}.fastq"
		if [[ -e "$src" ]]; then
			mv -f "$src" "$OUTPUT_DIR/${sample}${suffix}.fastq"
		fi
	done
done <"$RUN_MAP"
rmdir "$TEMP_DIR" 2>/dev/null || true

echo "[Step 6] Validating FASTQ outputs..."
missing=0
while IFS=$'\t' read -r accession sample; do
	found=false
	for fastq in "$OUTPUT_DIR/${sample}.fastq" \
				"$OUTPUT_DIR/${sample}_1.fastq" \
				"$OUTPUT_DIR/${sample}_2.fastq"; do
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
		echo "ERROR: No non-empty FASTQ output found for $sample ($accession)." >&2
		missing=1
	fi
done <"$RUN_MAP"
if [[ "$missing" -ne 0 ]]; then
	exit 1
fi

if [[ "$CLEANUP" == true ]]; then
	echo "[Step 7] Cleaning up SRA cache files..."
	while IFS= read -r accession; do
		rm -rf "${PREFETCH_DIR:?}/$accession"
	done <"$ACCESSION_LIST"
	rmdir "$PREFETCH_DIR" 2>/dev/null || true
else
	echo "[Step 7] Keeping SRA cache. Re-run with --cleanup to remove it."
fi

echo "[Step 8] Generating sample_list.txt..."
cut -f2 "$RUN_MAP" | sort -u >"$SAMPLE_LIST"

echo ""
echo "Download summary"
echo "================"
echo "SRA runs:        $(wc -l <"$RUN_MAP")"
echo "FASTQ files:     $(find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.fastq' | wc -l)"
echo "Sample list:     $SAMPLE_LIST"
echo "SRR map:         $RUN_MAP"
echo "FASTQ directory: $OUTPUT_DIR"
echo "Done."
