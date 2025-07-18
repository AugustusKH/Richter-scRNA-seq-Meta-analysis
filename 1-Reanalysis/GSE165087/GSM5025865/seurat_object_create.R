# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025865.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025865/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025865.seu <- CreateSeuratObject(counts = GSM5025865.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025865.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025865.seu[['tissue']] <- 'PB'
GSM5025865.seu[['case']] <- 1
GSM5025865.seu[['cond']] <- 'CLL'
GSM5025865.seu[['sorted']] <- 'whole'
GSM5025865.seu[['stage']] <- 'FCR_relapse_timepoint5'
GSM5025865.seu[['time']] <- 'Relapse_post_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025865.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025865/GSM5025865.rds')
