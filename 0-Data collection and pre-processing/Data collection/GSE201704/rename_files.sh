#!/bin/bash

# Going to working path
cd /home/pakorns/fastq_files/scRNAseq/GSE201704/fastq/sample_count

# Make directories
mkdir ADT GEX VDJ

# Change the rename of fastq files 
## VDJ files
cp ../SRR18938655_1.fastq ./VDJ/VDJ-Lane-1_S1_L001_I1_001.fastq.gz
cp ../SRR18938655_2.fastq ./VDJ/VDJ-Lane-1_S1_L001_R1_001.fastq.gz
cp ../SRR18938655_3.fastq ./VDJ/VDJ-Lane-1_S1_L001_R2_001.fastq.gz
cp ../SRR18938656_1.fastq ./VDJ/VDJ-Lane-5_S5_L002_I1_001.fastq.gz
cp ../SRR18938656_2.fastq ./VDJ/VDJ-Lane-5_S5_L002_R1_001.fastq.gz
cp ../SRR18938656_3.fastq ./VDJ/VDJ-Lane-5_S5_L002_R2_001.fastq.gz
cp ../SRR18938657_1.fastq ./VDJ/VDJ-Lane-6_S6_L001_I1_001.fastq.gz
cp ../SRR18938657_2.fastq ./VDJ/VDJ-Lane-6_S6_L001_R1_001.fastq.gz
cp ../SRR18938657_3.fastq ./VDJ/VDJ-Lane-6_S6_L001_R2_001.fastq.gz
cp ../SRR18938658_1.fastq ./VDJ/VDJ-Lane-6_S6_L002_I1_001.fastq.gz
cp ../SRR18938658_2.fastq ./VDJ/VDJ-Lane-6_S6_L002_R1_001.fastq.gz
cp ../SRR18938658_3.fastq ./VDJ/VDJ-Lane-6_S6_L002_R2_001.fastq.gz
cp ../SRR18938659_1.fastq ./VDJ/VDJ-Lane-7_S7_L001_I1_001.fastq.gz
cp ../SRR18938659_2.fastq ./VDJ/VDJ-Lane-7_S7_L001_R1_001.fastq.gz
cp ../SRR18938659_3.fastq ./VDJ/VDJ-Lane-7_S7_L001_R2_001.fastq.gz
cp ../SRR18938660_1.fastq ./VDJ/VDJ-Lane-7_S7_L002_I1_001.fastq.gz
cp ../SRR18938660_2.fastq ./VDJ/VDJ-Lane-7_S7_L002_R1_001.fastq.gz
cp ../SRR18938660_3.fastq ./VDJ/VDJ-Lane-7_S7_L002_R2_001.fastq.gz
cp ../SRR18938661_1.fastq ./VDJ/VDJ-Lane-8_S8_L001_I1_001.fastq.gz
cp ../SRR18938661_2.fastq ./VDJ/VDJ-Lane-8_S8_L001_R1_001.fastq.gz
cp ../SRR18938661_3.fastq ./VDJ/VDJ-Lane-8_S8_L001_R2_001.fastq.gz
cp ../SRR18938662_1.fastq ./VDJ/VDJ-Lane-8_S8_L002_I1_001.fastq.gz
cp ../SRR18938662_2.fastq ./VDJ/VDJ-Lane-8_S8_L002_R1_001.fastq.gz
cp ../SRR18938662_3.fastq ./VDJ/VDJ-Lane-8_S8_L002_R2_001.fastq.gz
cp ../SRR18938663_1.fastq ./VDJ/VDJ-Lane-1_S1_L002_I1_001.fastq.gz
cp ../SRR18938663_2.fastq ./VDJ/VDJ-Lane-1_S1_L002_R1_001.fastq.gz
cp ../SRR18938663_3.fastq ./VDJ/VDJ-Lane-1_S1_L002_R2_001.fastq.gz
cp ../SRR18938664_1.fastq ./VDJ/VDJ-Lane-2_S2_L001_I1_001.fastq.gz
cp ../SRR18938664_2.fastq ./VDJ/VDJ-Lane-2_S2_L001_R1_001.fastq.gz
cp ../SRR18938664_3.fastq ./VDJ/VDJ-Lane-2_S2_L001_R2_001.fastq.gz
cp ../SRR18938665_1.fastq ./VDJ/VDJ-Lane-2_S2_L002_I1_001.fastq.gz
cp ../SRR18938665_2.fastq ./VDJ/VDJ-Lane-2_S2_L002_R1_001.fastq.gz
cp ../SRR18938665_3.fastq ./VDJ/VDJ-Lane-2_S2_L002_R2_001.fastq.gz
cp ../SRR18938666_1.fastq ./VDJ/VDJ-Lane-3_S3_L001_I1_001.fastq.gz
cp ../SRR18938666_2.fastq ./VDJ/VDJ-Lane-3_S3_L001_R1_001.fastq.gz
cp ../SRR18938666_3.fastq ./VDJ/VDJ-Lane-3_S3_L001_R2_001.fastq.gz
cp ../SRR18938667_1.fastq ./VDJ/VDJ-Lane-3_S3_L002_I1_001.fastq.gz
cp ../SRR18938667_2.fastq ./VDJ/VDJ-Lane-3_S3_L002_R1_001.fastq.gz
cp ../SRR18938667_3.fastq ./VDJ/VDJ-Lane-3_S3_L002_R2_001.fastq.gz
cp ../SRR18938668_1.fastq ./VDJ/VDJ-Lane-4_S4_L001_I1_001.fastq.gz
cp ../SRR18938668_2.fastq ./VDJ/VDJ-Lane-4_S4_L001_R1_001.fastq.gz
cp ../SRR18938668_3.fastq ./VDJ/VDJ-Lane-4_S4_L001_R2_001.fastq.gz
cp ../SRR18938669_1.fastq ./VDJ/VDJ-Lane-4_S4_L002_I1_001.fastq.gz
cp ../SRR18938669_2.fastq ./VDJ/VDJ-Lane-4_S4_L002_R1_001.fastq.gz
cp ../SRR18938669_3.fastq ./VDJ/VDJ-Lane-4_S4_L002_R2_001.fastq.gz
cp ../SRR18938670_1.fastq ./VDJ/VDJ-Lane-5_S5_L001_I1_001.fastq.gz
cp ../SRR18938670_2.fastq ./VDJ/VDJ-Lane-5_S5_L001_R1_001.fastq.gz
cp ../SRR18938670_3.fastq ./VDJ/VDJ-Lane-5_S5_L001_R2_001.fastq.gz

# GEX files



