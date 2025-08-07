#!/bin/bash

#SBATCH --job-name=python_job         # Job name
#SBATCH --output=output_%j.txt        # Stdout log
#SBATCH --error=error_%j.txt          # Stderr log
#SBATCH --time=15-00:00:00            # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=64G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load Python module
module load python/3.10.0

# Base path to your expression matrices
BASE_DIR="/home/pakorns/hp-storage/scRNAseq/Nadeu2022_NatMed_scRNAseq_data/expression_matrices"

# List of subprojects
subprojects=("BCLLATLAS_10" "BCLLATLAS_29")

# Loop through each subproject
for subproject in "${subprojects[@]}"; do
    subproject_path="$BASE_DIR/$subproject"

    # Check if subproject path exists
    if [ -d "$subproject_path" ]; then
        # Loop through each GEM ID folder inside the subproject
        for gem_id in "$subproject_path"/*; do
            if [ -d "$gem_id" ]; then
                gem_id_name=$(basename "$gem_id")
                echo "Running doublet_scrublet.py for $subproject and $gem_id_name"
                python3 doublet_scrublet.py "$subproject" "$gem_id_name"
            fi
        done
    else
        echo "Subproject directory not found: $subproject_path"
    fi
done
