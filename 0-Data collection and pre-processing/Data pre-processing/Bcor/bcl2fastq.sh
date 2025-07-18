#!/bin/bash

#SBATCH --job-name=cellranger_mkfastq     # Job name
#SBATCH --output=output.txt               # Stdout log
#SBATCH --error=error.txt                 # Stderr log
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # CPU cores per task
#SBATCH --mem=128G                        # Memory per node
#SBATCH --time=15-00:00:00                # D-HH:MM:SS
#SBATCH --partition=hi                    # Partition name (check that 'hi' is valid)

# Increase file descriptor limit to avoid "too many open files" error
ulimit -n 65536

# Define paths
base_path="/home/pakorns/hp-storage/scRNAseq/Bcor"
SAMPLE_SHEET_PATH="${base_path}/samplesheet_test.csv"
OUTPUT_DIR="${base_path}/outputs"
FLOWCELL_DIR="${base_path}/Bcor_10x_genomics_raw_data"
INTEROP_DIR="${base_path}/stats"

# Create output directories if they don’t exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$INTEROP_DIR"

# Load required modules
module load bcl2fastq/2.20.0

# Run bcl2fastq
bcl2fastq \
  --use-bases-mask=Y28,I10,I10,Y90 \
  --create-fastq-for-index-reads \
  --minimum-trimmed-read-length=8 \
  --mask-short-adapter-reads=8 \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  --sample-sheet="${SAMPLE_SHEET_PATH}" \
  --runfolder-dir="${FLOWCELL_DIR}" \
  --output-dir="${OUTPUT_DIR}" \
  --interop-dir="${INTEROP_DIR}" \
  -r 6 \
  -w 6
