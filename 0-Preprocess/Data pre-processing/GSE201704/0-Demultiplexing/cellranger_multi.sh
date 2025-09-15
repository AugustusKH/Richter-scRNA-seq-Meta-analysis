#!/bin/bash

#SBATCH --job-name=cellranger_run     # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=15-00:00:00            # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=128G                    # Memory per node
#SBATCH --partition=hi                # Partition name

# Load cellranger (adjust if needed)
# module load cellranger/6.1.2
export PATH=$HOME/local/cellranger-9.0.1:$PATH

# Change to working directory
cd /home/pakorns/fastq_files/scRNAseq/GSE201704/fastq/sample_count

# Loop through lanes 1 to 8
for lane in {1..8}; do
    echo "Running Cell Ranger for lane $lane"
    cellranger multi --id=demultiplexed_samples_lane${lane} --csv=config_lane${lane}.csv
done
