# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025870.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025870/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025870.seu <- CreateSeuratObject(counts = GSM5025870.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025870.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025870.seu[['tissue']] <- 'PB'
GSM5025870.seu[['case']] <- 3
GSM5025870.seu[['cond']] <- 'CLL'
GSM5025870.seu[['sorted']] <- 'whole'
GSM5025870.seu[['stage']] <- 'FCR_relapse_timepoint2'
GSM5025870.seu[['time']] <- 'Relapse_post_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025870.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025870/GSM5025870.rds')
