# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025883.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025883/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025883.seu <- CreateSeuratObject(counts = GSM5025883.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025883.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025883.seu[['tissue']] <- 'PB'
GSM5025883.seu[['case']] <- 8
GSM5025883.seu[['cond']] <- 'CLL'
GSM5025883.seu[['sorted']] <- 'whole'
GSM5025883.seu[['stage']] <- 'Ibrutinib_timepoint2'
GSM5025883.seu[['time']] <- 'Ibrutinib_7_month'

# Save the object as .rds file
SaveSeuratRds(GSM5025883.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025883/GSM5025883.rds')
