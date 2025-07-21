# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
CMO306.data <- Read10X(data.dir = '/home/pakorns/hp-storage/scRNAseq/GSE201704/GSM6069857/CMO306/count_matrices/count/filtered_feature_bc_matrix')

# Create Seurat object
CMO306.seu <- CreateSeuratObject(counts = CMO306.data[['Gene Expression']], min.features = 500, min.cells = 1)

# QC
CMO306.seu[["percent.mt"]] <- PercentageFeatureSet(CMO306.seu, pattern = "^MT-")
CMO306.seu <- subset(CMO306.seu, subset = percent.mt < 20)

# Define other columns in meta-data
CMO306.seu[['HTO']] <- 'CMO306'
CMO306.seu[['ref']] <- 'Rejeski et al. (2022)'
CMO306.seu[['tissue']] <- 'PB'
CMO306.seu[['case']] <- 1
CMO306.seu[['cond']] <- 'RT'
CMO306.seu[['sorted']] <- 'whole'

# Save the object as .rds file
SaveSeuratRds(CMO306.seu, file='/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/CMO306/CMO306.rds')
