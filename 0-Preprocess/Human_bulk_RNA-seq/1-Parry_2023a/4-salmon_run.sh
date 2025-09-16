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
working_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/phs002458.v2.p1/fastq"
output_dir=""
transcript_ref_dir="/home/pakorns/hc-storage/fastq_files/RNAseq/transcriptome_ref"

# Download human reference transcriptome and build an index
cd "$transcript_ref_dir"
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i human_index
