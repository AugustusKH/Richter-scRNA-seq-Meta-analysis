#!/bin/bash

#SBATCH --job-name=seurat	      # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=15-00:00:00            # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=64G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load R module
module load R/4.4.2

# Run the R script
Rscript merge_object.R
Rscript preprocess_merge.R
