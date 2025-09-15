#!/bin/bash

#SBATCH --job-name=cellranger_run     # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=15-00:00:00            # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=64G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load cellranger
#module load cellranger/6.1.2
export PATH=$HOME/local/cellranger-9.0.1:$PATH

# Base path to main directory
cd /home/pakorns/fastq_files/scRNAseq/GSE201704/fastq/sample_count/data_aggr/CMO305

# Run cellranger multi
cellranger aggr --id=CMO305_aggr --csv=CMO305_aggr.csv
