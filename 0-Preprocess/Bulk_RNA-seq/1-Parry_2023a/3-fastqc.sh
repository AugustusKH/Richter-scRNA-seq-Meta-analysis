#!/bin/bash

#SBATCH --job-name=fastqc             # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=64G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load modules
module load fastqc/0.11.9
module load java/21.0.2+13

# Define paths
working_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/phs002458.v2.p1/fastq"

# List of SRR directories to process
selected_dirs=("SRR22839006" "SRR22839010" "SRR22839013" "SRR22839018" "SRR22956586" "SRR22960566" "SRR22960567" "SRR22964505" "SRR22964506" "SRR22964520")

# Loop for fastqc for each SRRs
for name in "${selected_dirs[@]}"; do
    dir="$working_dir/$name"
    output_dir="$dir/fastqc_results"
    mkdir -p "$output_dir"

    echo "🔄 Processing $dir"

    # Run fastqc
    fastqc -o "$output_dir" "$dir/${name}_1.fastq.gz" "$dir/${name}_2.fastq.gz"

    echo "📦 fastqc for $name is complete!"
done

echo "🎉 All jobs completed."
