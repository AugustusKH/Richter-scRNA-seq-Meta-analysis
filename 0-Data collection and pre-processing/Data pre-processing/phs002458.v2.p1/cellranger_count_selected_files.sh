#!/bin/bash

#SBATCH --job-name=cellranger_count   # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=32G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load cellranger
module load cellranger/6.1.2

# Base path to main directory
base_path="/home/pakorns/fastq_files/scRNAseq/phs002458.v2.p1/fastq"

# List of GSM directories
selected_dirs=("SRR22882959" "SRR22882960" "SRR22882961" "SRR22882962" "SRR22882969")

for name in "${selected_dirs[@]}"; do
    dir="$base_path/$name"
    if [ -d "$dir" ]; then
        echo "Processing $dir"
        cd "$dir"

	# Run cell ranger
	cellranger count --id=run_count \
                         --fastqs="$dir" \
                         --sample=sample \
                         --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A
       
    else
        echo "⚠️  Directory not found: $dir"
    fi
done
