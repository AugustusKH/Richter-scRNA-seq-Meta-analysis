# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025873.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025873/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025873.seu <- CreateSeuratObject(counts = GSM5025873.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025873.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025873.seu[['tissue']] <- 'PB'
GSM5025873.seu[['case']] <- 5
GSM5025873.seu[['cond']] <- 'CLL'
GSM5025873.seu[['sorted']] <- 'whole'
GSM5025873.seu[['stage']] <- 'RIC_relapse_timepoint1'
GSM5025873.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025873.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025873/GSM5025873.rds')
