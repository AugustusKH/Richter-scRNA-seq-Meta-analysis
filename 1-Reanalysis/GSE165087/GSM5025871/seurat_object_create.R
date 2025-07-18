# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025871.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025871/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025871.seu <- CreateSeuratObject(counts = GSM5025871.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025871.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025871.seu[['tissue']] <- 'PB'
GSM5025871.seu[['case']] <- 4
GSM5025871.seu[['cond']] <- 'CLL'
GSM5025871.seu[['sorted']] <- 'whole'
GSM5025871.seu[['stage']] <- 'RIC_relapse_timepoint1'
GSM5025871.seu[['time']] <- 'Pretreatment_FCR'

# Save the object as .rds file
SaveSeuratRds(GSM5025871.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025871/GSM5025871.rds')
