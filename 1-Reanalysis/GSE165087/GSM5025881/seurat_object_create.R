# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025881.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025881/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025881.seu <- CreateSeuratObject(counts = GSM5025881.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025881.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025881.seu[['tissue']] <- 'PB'
GSM5025881.seu[['case']] <- 7
GSM5025881.seu[['cond']] <- 'CLL'
GSM5025881.seu[['sorted']] <- 'whole'
GSM5025881.seu[['stage']] <- 'CLL Ibrutinib_timepoint3'
GSM5025881.seu[['time']] <- 'Ibrutinib_16_month'

# Save the object as .rds file
SaveSeuratRds(GSM5025881.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025881/GSM5025881.rds')
