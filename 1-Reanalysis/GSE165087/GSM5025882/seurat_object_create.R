# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025882.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025882/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025882.seu <- CreateSeuratObject(counts = GSM5025882.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025882.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025882.seu[['tissue']] <- 'PB'
GSM5025882.seu[['case']] <- 8
GSM5025882.seu[['cond']] <- 'CLL'
GSM5025882.seu[['sorted']] <- 'whole'
GSM5025882.seu[['stage']] <- 'Ibrutinib_timepoint1'
GSM5025882.seu[['time']] <- 'Pretreatment_Ibrutinib'

# Save the object as .rds file
SaveSeuratRds(GSM5025882.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025882/GSM5025882.rds')
