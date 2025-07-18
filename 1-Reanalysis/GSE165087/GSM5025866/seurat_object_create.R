# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025866.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025866/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025866.seu <- CreateSeuratObject(counts = GSM5025866.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025866.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025866.seu[['tissue']] <- 'PB'
GSM5025866.seu[['case']] <- 1
GSM5025866.seu[['cond']] <- 'CLL'
GSM5025866.seu[['sorted']] <- 'whole'
GSM5025866.seu[['stage']] <- 'FCR_relapse_timepoint5'
GSM5025866.seu[['time']] <- 'Relapse_post_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025866.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025866/GSM5025866.rds')
