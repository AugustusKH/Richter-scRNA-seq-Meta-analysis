# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025887.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025887/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025887.seu <- CreateSeuratObject(counts = GSM5025887.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025887.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025887.seu[['tissue']] <- 'BM'
GSM5025887.seu[['case']] <- 9
GSM5025887.seu[['cond']] <- 'RT'
GSM5025887.seu[['sorted']] <- 'whole'
GSM5025887.seu[['stage']] <- 'Richter_timepoint2'
GSM5025887.seu[['time']] <- 'BM_RT_phase'

# Save the object as .rds file
SaveSeuratRds(GSM5025887.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025887/GSM5025887.rds')
