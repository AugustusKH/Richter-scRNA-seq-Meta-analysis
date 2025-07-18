# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025862.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025862/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025862.seu <- CreateSeuratObject(counts = GSM5025862.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025862.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025862.seu[['tissue']] <- 'PB'
GSM5025862.seu[['case']] <- 1
GSM5025862.seu[['cond']] <- 'CLL'
GSM5025862.seu[['sorted']] <- 'whole'
GSM5025862.seu[['stage']] <- 'FCR_relapse_timepoint1'
GSM5025862.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025862.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM5025862/all/GSM5025862.rds')
