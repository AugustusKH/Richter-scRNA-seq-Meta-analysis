#!/bin/bash
#SBATCH --job-name=salmon               # Job name
#SBATCH --output=output.txt             # Stdout log
#SBATCH --error=error.txt               # Stderr log
#SBATCH --time=2-00:00:00               # D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=128G                      # Memory per node
#SBATCH --partition=low                 # Partition name

# Load modules
module load salmon/1.5.1
module load gcc/12.1.0

# Define paths
fastq_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/GSE183425"
transcript_ref_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/transcriptome_ref"
index_path="$transcript_ref_dir/mouse_index"

# Make sure transcript_ref_dir exists
mkdir -p "$transcript_ref_dir"

# Download transcriptome if not already present
if [ ! -f "$transcript_ref_dir/Mus_musculus.GRCm39.cdna.all.fa.gz" ]; then
    wget -O "$transcript_ref_dir/Mus_musculus.GRCm39.cdna.all.fa.gz" \
    https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
fi

# Build index if not already present
if [ ! -d "$index_path" ]; then
    salmon index -t "$transcript_ref_dir/Mus_musculus.GRCm39.cdna.all.fa.gz" \
                 -i "$index_path"
fi

# Loop for salmon
for i in `seq 497 504`; do
    name="GSM5557${i}"
    dir="$fastq_dir/$name"
    output_dir="$dir/quant"
    mkdir -p "$output_dir"

    # Get FASTQ files
    R1=$dir/SRR*_1.fastq.gz
    R2=$dir/SRR*_2.fastq.gz

    echo "🔄 Processing $dir"
    
    # Quantifying the samples
    salmon quant -i "$index_path" -l A\
	    -1 $R1 \
	    -2 $R2 \
	    -p 4 --validateMappings --gcBias --seqBias -o "$output_dir" 
    echo "📦 successful quantifying for $name"
done
