#!/bin/bash
#SBATCH --job-name=sra_to_fastq         # Job name
#SBATCH --output=output.txt             # Stdout log
#SBATCH --error=error.txt               # Stderr log
#SBATCH --time=2-00:00:00               # D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=32G                       # Memory per node
#SBATCH --partition=low                 # Partition name

# Load SRA Toolkit
module load sra-tools/3.0.3

# Define input and output directories
PREFETCH_DIR="/home/pakorns/fastq_files/scRNAseq/phs002458.v2.p1"
OUTER_OUTPUT_DIR="/home/pakorns/fastq_files/scRNAseq/phs002458.v2.p1/fastq"

# Create output directory
mkdir -p "$OUTER_OUTPUT_DIR"

# Loop through SRR directories
for DIR in "$PREFETCH_DIR"/SRR*/
do
    echo "Processing $DIR..."

    # Extract SRR ID and create subdirectory
    SRR_ID=$(basename "$DIR")
    INNER_OUTPUT_DIR="$OUTER_OUTPUT_DIR/$SRR_ID"
    mkdir -p "$INNER_OUTPUT_DIR"

    # Run fasterq-dump
    fasterq-dump --split-files --include-technical "$DIR" -O "$INNER_OUTPUT_DIR" -e 10

    # Check dump status before compressing
    if [ $? -eq 0 ]; then
        echo "✅ Successfully dumped $SRR_ID"

        # Compress resulting FASTQ files
        gzip "$INNER_OUTPUT_DIR"/*.fastq
	
	# Rename compressed FASTQ files for Cell Ranger
	mv "$INNER_OUTPUT_DIR/${SRR_ID}_1.fastq.gz" "$INNER_OUTPUT_DIR/sample_S1_L001_I1_001.fastq.gz"
	mv "$INNER_OUTPUT_DIR/${SRR_ID}_2.fastq.gz" "$INNER_OUTPUT_DIR/sample_S1_L001_R1_001.fastq.gz"
	mv "$INNER_OUTPUT_DIR/${SRR_ID}_3.fastq.gz" "$INNER_OUTPUT_DIR/sample_S1_L001_R2_001.fastq.gz"

        echo "📦 FASTQ files compressed for $SRR_ID"
    else
        echo "❌ Error processing $SRR_ID" >&2
    fi
done

echo "🎉 All jobs completed."
