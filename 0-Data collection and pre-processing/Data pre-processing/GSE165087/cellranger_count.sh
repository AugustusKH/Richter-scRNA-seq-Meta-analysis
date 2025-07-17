#!/bin/bash

#SBATCH --job-name=cellranger_run     # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=16G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load cellranger
module load cellranger/6.1.2

# Base path to main directory
base_path="/home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices"

# List of GSM directories
selected_dirs=("GSM5025861" "GSM5025862" "GSM5025863" "GSM5025864" "GSM5025865" "GSM5025866" "GSM5025867" "GSM5025868" "GSM5025869" "GSM5025870" "GSM5025871" "GSM5025871" "GSM5025872" "GSM5025873" "GSM5025874" "GSM5025875" "GSM5025876" "GSM5025877" "GSM5025878" "GSM5025881" "GSM5025882" "GSM5025883" "GSM5025884" "GSM5025885" "GSM5025886" "GSM5025887")

for name in "${selected_dirs[@]}"; do
    dir="$base_path/$name"
    if [ -d "$dir" ]; then
        echo "Processing $dir"
        cd "$dir"

        # Check that the Pool directory exists
        fastq_dir=$(find bamtofastq_output -type d -name "Pool*" | head -n 1)

        if [ -n "$fastq_dir" ]; then
            cellranger count --id=run_count \
                --fastqs="$dir/$fastq_dir" \
                --sample=bamtofastq \
                --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A
        else
            echo "⚠️  No matching Pool* directory found in $dir/bamtofastq_output"
        fi
    else
        echo "⚠️  Directory not found: $dir"
    fi
done
