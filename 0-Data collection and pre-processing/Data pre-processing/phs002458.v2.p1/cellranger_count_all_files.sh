#!/bin/bash

#SBATCH --job-name=cellranger_count   # Job name
#SBATCH --output=output.txt           # Stdout log for array jobs
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=32G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load cellranger
module load cellranger/6.1.2

# Base directory
base_path="/home/pakorns/fastq_files/scRNAseq/phs002458.v2.p1/fastq"

# Loop through all SRR directories
for dir in "$base_path"/SRR*; do
    if [ -d "$dir" ]; then
        srr_id=$(basename "$dir")
        echo "🔧 Processing $srr_id"

        # Create unique output folder per SRR
        output_id="run_count_${srr_id}"

        # Run cellranger count
        cellranger count \
            --id="$output_id" \
            --fastqs="$dir" \
            --sample=sample \
            --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A
    else
        echo "⚠️  Skipped: Not a directory - $dir"
    fi
done
