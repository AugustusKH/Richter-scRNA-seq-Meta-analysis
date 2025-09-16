#!/bin/bash
#SBATCH --job-name=salmon               # Job name
#SBATCH --output=output.txt             # Stdout log
#SBATCH --error=error.txt               # Stderr log
#SBATCH --time=2-00:00:00               # D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores per task
#SBATCH --mem=64G                       # Memory per node
#SBATCH --partition=low                 # Partition name

# Load modules
module load salmon/1.5.1
module load gcc/12.1.0

# Define paths
fastq_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/phs002458.v2.p1/fastq"
transcript_ref_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/transcriptome_ref"
index_path="$transcript_ref_dir/human_index"

# Make sure transcript_ref_dir exists
mkdir -p "$transcript_ref_dir"

# Download transcriptome if not already present
if [ ! -f "$transcript_ref_dir/Homo_sapiens.GRCh38.cdna.all.fa.gz" ]; then
    wget -O "$transcript_ref_dir/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
    https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

# Build index if not already present
if [ ! -d "$index_path" ]; then
    salmon index -t "$transcript_ref_dir/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
                 -i "$index_path"
fimon index -t "$transcript_ref_dir/Homo_sapiens.GRCh38.cdna.all.fa.gz" -i "$transcript_ref_dir/human_index"

# Define list
selected_dir=("SRR22839006" "SRR22839010" "SRR22839013" "SRR22839018" "SRR22956586" "SRR22960566" "SRR22960567" "SRR22964505" "SRR22964506" "SRR22964520")

# Loop for salmon
for name in "${selected_dir[@]}"; do
    dir="$fastq_dir/$name"
    output_dir="$dir/quant"
    mkdir -p "$output_dir"

    # Check FASTQs exist
    if [[ ! -f "$dir/${name}_1.fastq.gz" || ! -f "$dir/${name}_2.fastq.gz" ]]; then
        echo "❌ Missing FASTQ files for $name, skipping."
        continue
    fi

    echo "🔄 Processing $dir"
    
    # Quantifying the samples
    salmon quant -i "$index_path" -l A\
	    -1 "$dir/${name}_1.fastq.gz" \
	    -2 "$dir/${name}_2.fastq.gz" \
	    -p 8 --validateMappings --gcBias --seqBias -o "$output_dir" 
    echo "📦 successful quantifying for $name"
done


