#!/bin/bash

# Load module
module load sra-tools/3.0.3

# Enter to working path
cd /home/pakorns/fastq_files/scRNAseq/GSE205704/sample_matrices

# Make directories
mkdir GSM6217975 GSM6217976 GSM6217978 GSM6217979 GSM6217981 GSM6217983

# Download fastq files according to samples
## GSM6217975
cd GSM6217975

fastq-dump --split-files SRR19586157
gzip SRR19586157_1.fastq SRR19586157_2.fastq SRR19586157_3.fastq SRR19586157_4.fastq
mv SRR19586157_1.fastq.gz L33_EW_GEX_S1_L001_I1_001.fastq.gz
mv SRR19586157_2.fastq.gz L33_EW_GEX_S1_L001_I2_001.fastq.gz
mv SRR19586157_3.fastq.gz L33_EW_GEX_S1_L001_R1_001.fastq.gz
mv SRR19586157_4.fastq.gz L33_EW_GEX_S1_L001_R2_001.fastq.gz

fastq-dump --split-files SRR19586158
gzip SRR19586158_1.fastq SRR19586158_2.fastq SRR19586158_3.fastq SRR19586158_4.fastq
mv SRR19586158_1.fastq.gz L33_EW_GEX_S1_L002_I1_001.fastq.gz
mv SRR19586158_2.fastq.gz L33_EW_GEX_S1_L002_I2_001.fastq.gz
mv SRR19586158_3.fastq.gz L33_EW_GEX_S1_L002_R1_001.fastq.gz
mv SRR19586158_4.fastq.gz L33_EW_GEX_S1_L002_R2_001.fastq.gz

cd ..

## GSM6217976
cd GSM6217976

fastq-dump --split-files SRR19586155
gzip SRR19586155_1.fastq SRR19586155_2.fastq SRR19586155_3.fastq SRR19586155_4.fastq
mv SRR19586155_1.fastq.gz L34_EW_GEX_S2_L001_I1_001.fastq.gz
mv SRR19586155_2.fastq.gz L34_EW_GEX_S2_L001_I2_001.fastq.gz
mv SRR19586155_3.fastq.gz L34_EW_GEX_S2_L001_R1_001.fastq.gz
mv SRR19586155_4.fastq.gz L34_EW_GEX_S2_L001_R2_001.fastq.gz

fastq-dump --split-files SRR19586156
gzip SRR19586156_1.fastq SRR19586156_2.fastq SRR19586156_3.fastq SRR19586156_4.fastq
mv SRR19586156_1.fastq.gz L34_EW_GEX_S2_L002_I1_001.fastq.gz
mv SRR19586156_2.fastq.gz L34_EW_GEX_S2_L002_I2_001.fastq.gz
mv SRR19586156_3.fastq.gz L34_EW_GEX_S2_L002_R1_001.fastq.gz
mv SRR19586156_4.fastq.gz L34_EW_GEX_S2_L002_R2_001.fastq.gz

cd ..

## GSM6217978
cd GSM6217978

fastq-dump --split-files SRR19586153
gzip SRR19586153_1.fastq SRR19586153_2.fastq SRR19586153_3.fastq SRR19586153_4.fastq
mv SRR19586153_1.fastq.gz L35_EW_GEX_S3_L001_I1_001.fastq.gz
mv SRR19586153_2.fastq.gz L35_EW_GEX_S3_L001_I2_001.fastq.gz
mv SRR19586153_3.fastq.gz L35_EW_GEX_S3_L001_R1_001.fastq.gz
mv SRR19586153_4.fastq.gz L35_EW_GEX_S3_L001_R2_001.fastq.gz

fastq-dump --split-files SRR19586154
gzip SRR19586154_1.fastq SRR19586154_2.fastq SRR19586154_3.fastq SRR19586154_4.fastq
mv SRR19586154_1.fastq.gz L35_EW_GEX_S3_L002_I1_001.fastq.gz
mv SRR19586154_2.fastq.gz L35_EW_GEX_S3_L002_I2_001.fastq.gz
mv SRR19586154_3.fastq.gz L35_EW_GEX_S3_L002_R1_001.fastq.gz
mv SRR19586154_4.fastq.gz L35_EW_GEX_S3_L002_R2_001.fastq.gz

cd ..

## GSM6217979
cd GSM6217979

fastq-dump --split-files SRR19586151
gzip SRR19586151_1.fastq SRR19586151_2.fastq SRR19586151_3.fastq SRR19586151_4.fastq
mv SRR19586151_1.fastq.gz L36_EW_GEX_S4_L001_I1_001.fastq.gz
mv SRR19586151_2.fastq.gz L36_EW_GEX_S4_L001_I2_001.fastq.gz
mv SRR19586151_3.fastq.gz L36_EW_GEX_S4_L001_R1_001.fastq.gz
mv SRR19586151_4.fastq.gz L36_EW_GEX_S4_L001_R2_001.fastq.gz

fastq-dump --split-files SRR19586152
gzip SRR19586152_1.fastq SRR19586152_2.fastq SRR19586152_3.fastq SRR19586152_4.fastq
mv SRR19586152_1.fastq.gz L36_EW_GEX_S4_L002_I1_001.fastq.gz
mv SRR19586152_2.fastq.gz L36_EW_GEX_S4_L002_I2_001.fastq.gz
mv SRR19586152_3.fastq.gz L36_EW_GEX_S4_L002_R1_001.fastq.gz
mv SRR19586152_4.fastq.gz L36_EW_GEX_S4_L002_R2_001.fastq.gz

cd ..

## GSM6217981
cd GSM6217981

fastq-dump --split-files SRR19586149
gzip SRR19586149_1.fastq SRR19586149_2.fastq SRR19586149_3.fastq SRR19586149_4.fastq
mv SRR19586149_1.fastq.gz 
mv SRR19586149_2.fastq.gz
mv SRR19586149_3.fastq.gz
mv SRR19586149_4.fastq.gz

fastq-dump --split-files SRR19586150
gzip SRR19586150_1.fastq SRR19586150_2.fastq SRR19586150_3.fastq SRR19586150_4.fastq
mv SRR19586150_1.fastq.gz
mv SRR19586150_2.fastq.gz
mv SRR19586150_3.fastq.gz
mv SRR19586150_4.fastq.gz

cd ..

## GSM6217983
cd GSM6217983

fastq-dump --split-files SRR19586147
gzip SRR19586147_1.fastq SRR19586147_2.fastq SRR19586147_3.fastq SRR19586147_4.fastq
mv SRR19586147_1.fastq.gz
mv SRR19586147_2.fastq.gz
mv SRR19586147_3.fastq.gz
mv SRR19586147_4.fastq.gz

fastq-dump --split-files SRR19586148
gzip SRR19586148_1.fastq SRR19586148_2.fastq SRR19586148_3.fastq SRR19586148_4.fastq
mv SRR19586148_1.fastq.gz
mv SRR19586148_2.fastq.gz
mv SRR19586148_3.fastq.gz
mv SRR19586148_4.fastq.gz

cd ..
