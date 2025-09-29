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
fastq_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/GSE263238/fastq_gzip"

# Loop for fastqc for each SRRs
for i in `seq 76 84`; do
    name="SRR285602${i}"
    dir="$fastq_dir/$name"
    output_dir="$dir/fastqc_results"
    mkdir -p "$output_dir"

    echo "🔄 Processing $dir"

    # Run fastqc
    fastqc -o "$output_dir" "$dir/${name}_1.fastq.gz" "$dir/${name}_2.fastq.gz"

    echo "📦 fastqc for $name is complete!"
done

echo "🎉 All jobs completed."
