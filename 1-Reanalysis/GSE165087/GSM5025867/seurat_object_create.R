# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025867.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025867/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025867.seu <- CreateSeuratObject(counts = GSM5025867.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025867.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025867.seu[['tissue']] <- 'PB'
GSM5025867.seu[['case']] <- 2
GSM5025867.seu[['cond']] <- 'CLL'
GSM5025867.seu[['sorted']] <- 'whole'
GSM5025867.seu[['stage']] <- 'FCR_relapse_timepoint1'
GSM5025867.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025867.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025867/GSM5025867.rds')
