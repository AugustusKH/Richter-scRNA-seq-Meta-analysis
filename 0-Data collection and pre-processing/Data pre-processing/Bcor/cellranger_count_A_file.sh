#!/bin/bash

#SBATCH --job-name=cellranger_count   # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=15-00:00:00            # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=128G                    # Memory per node
#SBATCH --partition=hi                # Partition name

# Load cellranger
module load cellranger/6.1.2

# Base path to main directory
base_path="/home/pakorns/hp-storage/scRNAseq/Bcor/outputs/CPT190_Ferrarini_10Xsc"

# List of GSM directories
selected_dirs=("A")

for name in "${selected_dirs[@]}"; do
    dir="$base_path/$name"
    if [ -d "$dir" ]; then
        echo "Processing $dir"
        cd "$dir"

	# Run cell ranger
	cellranger count --id=run_count \
                         --fastqs="$dir" \
                         --sample=Pre \
                         --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A

    else
        echo "⚠️  Directory not found: $dir"
    fi
done
