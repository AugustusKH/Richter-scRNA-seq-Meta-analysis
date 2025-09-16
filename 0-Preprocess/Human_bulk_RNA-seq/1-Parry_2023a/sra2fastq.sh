#!/bin/bash
#SBATCH --job-name=sra2fastq            # Job name
#SBATCH --output=output.txt             # Stdout log
#SBATCH --error=error.txt               # Stderr log
#SBATCH --time=2-00:00:00               # D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=64G                       # Memory per node
#SBATCH --partition=low                 # Partition name

# Load SRA Toolkit
module load sra-tools/3.0.3

# Define input and output directories
PREFETCH_DIR="/home/pakorns/hc-storage/fastq_files/RNAseq/phs002458.v2.p1"
base_path="/home/pakorns/hc-storage/fastq_files/RNAseq/phs002458.v2.p1/fastq"

# List of SRR directories to process
selected_dirs=("SRR22839006" "SRR22839010" "SRR22839013" "SRR22839018" "SRR22956586" "SRR22960566" "SRR22960567" "SRR22964505" "SRR22964506" "SRR22964520")

# Loop through selected SRR IDs
for name in "${selected_dirs[@]}"; do
    dir="$base_path/$name"
    echo "🔄 Processing $dir"

    # Create output directory
    mkdir -p "$dir"

    # Run fasterq-dump
    fasterq-dump --split-files --include-technical "$PREFETCH_DIR/$name" -O "$dir" -e 10

    # Check if the dump was successful
    if [ $? -eq 0 ]; then
        echo "✅ Successfully dumped $name"

        # Compress resulting FASTQ files
        gzip "$dir"/*.fastq
        
        echo "📦 FASTQ files compressed and renamed for $name"
    else
        echo "❌ Error processing $name" >&2
    fi
done

echo "🎉 All jobs completed."
