# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025869.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025869/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025869.seu <- CreateSeuratObject(counts = GSM5025869.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025869.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025869.seu[['tissue']] <- 'PB'
GSM5025869.seu[['case']] <- 3
GSM5025869.seu[['cond']] <- 'CLL'
GSM5025869.seu[['sorted']] <- 'whole'
GSM5025869.seu[['stage']] <- 'FCR_relapse_timepoint1'
GSM5025869.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025869.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025869/GSM5025869.rds')
