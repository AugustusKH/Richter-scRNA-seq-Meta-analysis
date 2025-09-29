#!/bin/bash
#SBATCH --job-name=gzip                 # Job name
#SBATCH --output=output.txt             # Stdout log
#SBATCH --error=error.txt               # Stderr log
#SBATCH --time=2-00:00:00               # D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=64G                       # Memory per node
#SBATCH --partition=low                 # Partition name

# Define input and output directories
FASTQ_DIR="/home/pakorns/hc-storage/fastq_files/RNAseq/GSE263238"
OUTPUT_DIR="/home/pakorns/hc-storage/fastq_files/RNAseq/GSE263238/fastq_gzip"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Loop through SRR fastq files
for FASTQ in "$FASTQ_DIR"/SRR*.fastq
do
    echo "Compressing $FASTQ..."

    # Compress and move to output directory
    gzip -c "$FASTQ" > "$OUTPUT_DIR/$(basename "$FASTQ").gz" 

    echo "📦 FASTQ files compressed for $FASTQ"
done

echo "🎉 All jobs completed."
