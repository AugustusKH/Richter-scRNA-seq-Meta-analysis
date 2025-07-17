#!/bin/bash

#SBATCH --job-name=bam_to_fastq 	# Job name
#SBATCH --output=output.txt		# Stdout log
#SBATCH --error=error.txt		# Stderr log
#SBATCH --time=2-00:00:00      		# D-HH:MM:SS (Max run time)
#SBATCH --ntasks=1          		# Number of tasks
#SBATCH --cpus-per-task=4   		# Number of CPU cores per task
#SBATCH --mem=16G            		# Memory per node
#SBATCH --partition=low     		# Partition name

# Enter to working path
cd /home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices

# Load the relevant module
module load cellranger/6.1.2

# Go to each bam files in GSM directories and convert to fastq files
cd GSM5025861
cellranger bamtofastq CLL1_1.1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025862
cellranger bamtofastq CLL1_1.2.possorted_genome.bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025863
cellranger bamtofastq CLL1_3.1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025864
cellranger bamtofastq CLL1_3.2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025865
cellranger bamtofastq CLL1_5.1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025866
cellranger bamtofastq CLL1_5.2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025867
cellranger bamtofastq CLL2_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025868
cellranger bamtofastq CLL2_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025869
cellranger bamtofastq CLL3_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025870
cellranger bamtofastq CLL3_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025871
cellranger bamtofastq CLL4_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025872
cellranger bamtofastq CLL4_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025873
cellranger bamtofastq CLL5_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025874
cellranger bamtofastq CLL5_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025875
cellranger bamtofastq CLL6_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025876
cellranger bamtofastq CLL6_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025877
cellranger bamtofastq CLL6_3.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025878
cellranger bamtofastq CLL6_4.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025879
cellranger bamtofastq CLL7_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025880
cellranger bamtofastq CLL7_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025881
cellranger bamtofastq CLL7_3.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025882
cellranger bamtofastq CLL8_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025883
cellranger bamtofastq CLL8_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025884
cellranger bamtofastq CLL8_3.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025885
cellranger bamtofastq CLL9_1.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025886
cellranger bamtofastq CLL9_2.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..

cd GSM5025887
cellranger bamtofastq CLL9_3.possorted_genome_bam.bam.1 ./bamtofastq_output
cd ..
