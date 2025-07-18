# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025872.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025872/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025872.seu <- CreateSeuratObject(counts = GSM5025872.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025872.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025872.seu[['tissue']] <- 'PB'
GSM5025872.seu[['case']] <- 4
GSM5025872.seu[['cond']] <- 'CLL'
GSM5025872.seu[['sorted']] <- 'whole'
GSM5025872.seu[['stage']] <- 'RIC_relapse_timepoint2'
GSM5025872.seu[['time']] <- 'Relapse_post_RIC'

# Save the object as .rds file
SaveSeuratRds(GSM5025872.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025872/GSM5025872.rds')
