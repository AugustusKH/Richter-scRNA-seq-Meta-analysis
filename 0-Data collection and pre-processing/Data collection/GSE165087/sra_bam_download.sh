# This script is used to download bamfiles from ncbi using sra-tools for GSE165087

# Load sra module
module load sra-tools/3.0.3

# Enter the working path
cd /home/pakorns/fastq_files/scRNAseq/GSE165087/sample_matrices

# Make directories and download bam files in individual GSM datasets
mkdir GSM5025861
cd GSM5025861
wget https://sra-pub-src-1.s3.amazonaws.com/SRR13482706/CLL1_1.1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025862
cd GSM5025862
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482707/CLL1_1.2.possorted_genome.bam.bam.1
cd ..

mkdir GSM5025863
cd GSM5025863
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482708/CLL1_3.1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025864
cd GSM5025864
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482709/CLL1_3.2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025865
cd GSM5025865
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482710/CLL1_5.1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025866
cd GSM5025866
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482711/CLL1_5.2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025867
cd GSM5025867
wget https://sra-pub-src-1.s3.amazonaws.com/SRR13482712/CLL2_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025868
cd GSM5025868
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482713/CLL2_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025869
cd GSM5025869
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482714/CLL3_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025870
cd GSM5025870
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482715/CLL3_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025871
cd GSM5025871
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482716/CLL4_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025872
cd GSM5025872
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482717/CLL4_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025873
cd GSM5025873
wget https://sra-pub-src-1.s3.amazonaws.com/SRR13482718/CLL5_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025874
cd GSM5025874
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482719/CLL5_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025875
cd GSM5025875
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482720/CLL6_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025876
cd GSM5025876
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482721/CLL6_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025877
cd GSM5025877
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482722/CLL6_3.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025878
cd GSM5025878
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482723/CLL6_4.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025879
cd GSM5025879
wget https://sra-pub-src-1.s3.amazonaws.com/SRR13482724/CLL7_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025880
cd GSM5025880
wgethttps://sra-pub-src-2.s3.amazonaws.com/SRR13482725/CLL7_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025881
cd GSM5025881
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482726/CLL7_3.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025882
cd GSM5025882
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482727/CLL8_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025883
cd GSM5025883
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482728/CLL8_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025884
cd GSM5025884
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482729/CLL8_3.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025885
cd GSM5025885
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482730/CLL9_1.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025886
cd GSM5025886
wget https://sra-pub-src-2.s3.amazonaws.com/SRR13482731/CLL9_2.possorted_genome_bam.bam.1
cd ..

mkdir GSM5025887
cd GSM5025887
wget https://sra-pub-src-1.s3.amazonaws.com/SRR13482732/CLL9_3.possorted_genome_bam.bam.1
cd ..

