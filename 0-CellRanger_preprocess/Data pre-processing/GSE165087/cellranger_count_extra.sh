#!/bin/bash

#SBATCH --job-name=cellranger_run     # Job name
#SBATCH --output=output.txt           # Stdout log
#SBATCH --error=error.txt             # Stderr log
#SBATCH --time=2-00:00:00             # D-HH:MM:SS
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # CPU cores per task
#SBATCH --mem=16G                     # Memory per node
#SBATCH --partition=low               # Partition name

# Load cellranger module
module load cellranger/6.1.2

# Go to working path and run cellranger count
cd /home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025879 

cellranger count --id=run_count \
  --fastqs=/home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025879/bamtofastq_output/SatIB1_0_1_HJC5JDMXX,/home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025879/bamtofastq_output/SatIB1_1_1_HJC5JDMXX \
  --sample=bamtofastq \
  --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A \
  --chemistry=SC5P-R2

cd /home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025880 

cellranger count --id=run_count \
  --fastqs=/home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025880/bamtofastq_output/SatIB2_0_1_HJC5JDMXX,/home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices/GSM5025880/bamtofastq_output/SatIB2_1_1_HJC5JDMXX \
  --sample=bamtofastq \
  --transcriptome=/home/pakorns/scRNAseq_practise/run_cellranger_count/refdata-gex-GRCh38-2020-A \
  --chemistry=SC5P-R2
 
