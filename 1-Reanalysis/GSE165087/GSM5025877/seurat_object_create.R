# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025877.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025877/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025861.seu <- CreateSeuratObject(counts = GSM5025861.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025877.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025877.seu[['tissue']] <- 'PB'
GSM5025877.seu[['case']] <- 6
GSM5025877.seu[['cond']] <- 'CLL'
GSM5025877.seu[['sorted']] <- 'whole'
GSM5025877.seu[['stage']] <- 'RIC/Ibrutinib_relapse_timepoint3'
GSM5025877.seu[['time']] <- 'Ibrutinib_1_month'

# Save the object as .rds file
SaveSeuratRds(GSM5025877.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025877/GSM5025877.rds')
