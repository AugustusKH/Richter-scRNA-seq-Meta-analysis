# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
CMO305.data <- Read10X(data.dir = '/home/pakorns/hp-storage/scRNAseq/GSE201704/GSM6069857/CMO305/count_matrices/count/filtered_feature_bc_matrix')

# Create Seurat object
CMO305.seu <- CreateSeuratObject(counts = CMO305.data[['Gene Expression']], min.features = 500, min.cells = 1)

# QC
CMO305.seu[["percent.mt"]] <- PercentageFeatureSet(CMO305.seu, pattern = "^MT-")
CMO305.seu <- subset(CMO305.seu, subset = percent.mt < 20)

# Define other columns in meta-data
CMO305.seu[['HTO']] <- 'CMO305'
CMO305.seu[['ref']] <- 'Rejeski et al. (2022)'
CMO305.seu[['tissue']] <- 'PB'
CMO305.seu[['case']] <- 1
CMO305.seu[['cond']] <- 'RT'
CMO305.seu[['sorted']] <- 'whole'

# Save the object as .rds file
SaveSeuratRds(CMO305.seu, file='/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/CMO305/CMO305.rds')
