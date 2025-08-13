#!/bin/bash

#SBATCH --job-name=cellranger_run 	# Job name
#SBATCH --output=output.txt		# Stdout log
#SBATCH --error=error.txt		# Stderr log
#SBATCH --time=2-00:00:00      		# D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1          		# Number of tasks
#SBATCH --cpus-per-task=4   		# Number of CPU cores per task
#SBATCH --mem=16G            		# Memory per node
#SBATCH --partition=low   		# Partition name

# Load cell ranger
module load cellranger/6.1.2

# Go to the working path
base_path="/home/pakorns/fastq_files/scRNAseq/GSE205704/sample_matrices"

# GSM IDs and matching sample prefixes
declare -A samples
samples=(
  ["GSM6217975"]="L33_EW_GEX"
  ["GSM6217976"]="L34_EW_GEX"
  ["GSM6217978"]="L35_EW_GEX"
  ["GSM6217979"]="L36_EW_GEX"
  ["GSM6217981"]="L37_EW_GEX"
  ["GSM6217983"]="L38_EW_GEX"
)

# Run cellranger count in each GSM datasets
for gsm in "${!samples[@]}"; do
    sample_name="${samples[$gsm]}"
    dir="$base_path/$gsm"
    if [ -d "$dir" ]; then
        echo "Processing $dir with sample $sample_name"
        cd "$dir"
        cellranger count --id=run_count \
            --fastqs="$dir" \
            --sample="$sample_name" \
            --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A
    else
        echo "⚠️  Directory not found: $dir"
    fi
done
