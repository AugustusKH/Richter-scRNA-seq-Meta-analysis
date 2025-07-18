# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025875.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025875/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025875.seu <- CreateSeuratObject(counts = GSM5025875.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025875.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025875.seu[['tissue']] <- 'PB'
GSM5025875.seu[['case']] <- 6
GSM5025875.seu[['cond']] <- 'CLL'
GSM5025875.seu[['sorted']] <- 'whole'
GSM5025875.seu[['stage']] <- 'RIC/Ibrutinib_relapse_timepoint1'
GSM5025875.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025875.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025875/GSM5025875.rds')
