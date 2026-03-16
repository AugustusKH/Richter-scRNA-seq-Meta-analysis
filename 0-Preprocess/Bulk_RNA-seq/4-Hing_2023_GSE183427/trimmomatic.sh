#!/bin/bash

#SBATCH --job-name=trimmomatic        # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=8             # CPU cores per task
#SBATCH --mem=64G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load modules
module load adoptopenjdk/11.0.12+7
module load trimmomatic/0.39

# Define paths
fastq_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/GSE183427"

# Loop for fastqc for each SRRs
for i in `seq 25 33`; do
    name="GSM55575${i}"
    dir="$fastq_dir/$name"
    output_dir="$dir/trimmed_results"
    mkdir -p "$output_dir"

    echo "🔄 Processing $dir"
    
    # Define input files
    R1="$dir"/SRR*_1.fastq.gz
    R2="$dir"/SRR*_2.fastq.gz

    # Extract sample ID (e.g., SRR15719797)
    sample=$(basename $R1 _1.fastq.gz)

    OUT1P="$output_dir/${sample}_1.trimmed.P.fastq.gz"   # paired output
    OUT1U="$output_dir/${sample}_1.trimmed.U.fastq.gz"   # unpaired output
    OUT2P="$output_dir/${sample}_2.trimmed.P.fastq.gz"
    OUT2U="$output_dir/${sample}_2.trimmed.U.fastq.gz"

    # Run Trimmomatic
    java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 8 -phred33 \
      $R1 $R2 \
      $OUT1P $OUT1U \
      $OUT2P $OUT2U \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

    echo "📦 trimmomatic for $name is complete!"
done

echo "🎉 All jobs completed."
