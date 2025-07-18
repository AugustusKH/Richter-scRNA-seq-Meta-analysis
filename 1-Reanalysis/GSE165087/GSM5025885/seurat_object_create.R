# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025885.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025885/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025885.seu <- CreateSeuratObject(counts = GSM5025885.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025885.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025885.seu[['tissue']] <- 'PB'
GSM5025885.seu[['case']] <- 9
GSM5025885.seu[['cond']] <- 'CLL'
GSM5025885.seu[['sorted']] <- 'whole'
GSM5025885.seu[['stage']] <- 'Richter_timepoint1'
GSM5025885.seu[['time']] <- 'PB_CLL_phase'

# Save the object as .rds file
SaveSeuratRds(GSM5025885.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025885/GSM5025885.rds')
