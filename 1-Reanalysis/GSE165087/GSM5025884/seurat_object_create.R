# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025884.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025884/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025884.seu <- CreateSeuratObject(counts = GSM5025884.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025884.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025884.seu[['tissue']] <- 'PB'
GSM5025884.seu[['case']] <- 8
GSM5025884.seu[['cond']] <- 'CLL'
GSM5025884.seu[['sorted']] <- 'whole'
GSM5025884.seu[['stage']] <- 'Ibrutinib_timepoint3'
GSM5025884.seu[['time']] <- 'Ibrutinib_14_month'

# Save the object as .rds file
SaveSeuratRds(GSM5025884.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025884/GSM5025884.rds')
