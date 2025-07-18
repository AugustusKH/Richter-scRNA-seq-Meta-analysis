# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025880.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025880/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025880.seu <- CreateSeuratObject(counts = GSM5025880.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025880.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025880.seu[['tissue']] <- 'PB'
GSM5025880.seu[['case']] <- 7
GSM5025880.seu[['cond']] <- 'CLL'
GSM5025880.seu[['sorted']] <- 'whole'
GSM5025880.seu[['stage']] <- 'Ibrutinib_timepoint2'
GSM5025880.seu[['time']] <- 'Ibrutinib_6_month'

# Save the object as .rds file
SaveSeuratRds(GSM5025880.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025880/GSM5025880.rds')
