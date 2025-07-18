# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025886.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025886/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025886.seu <- CreateSeuratObject(counts = GSM5025886.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025886.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025886.seu[['tissue']] <- 'PB'
GSM5025886.seu[['case']] <- 9
GSM5025886.seu[['cond']] <- 'RT'
GSM5025886.seu[['sorted']] <- 'whole'
GSM5025886.seu[['stage']] <- 'Richter_timepoint2 '
GSM5025886.seu[['time']] <- 'PB_RT_phase'

# Save the object as .rds file
SaveSeuratRds(GSM5025886.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025886/GSM5025886.rds')
