# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025879.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025879/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025879.seu <- CreateSeuratObject(counts = GSM5025879.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025879.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025879.seu[['tissue']] <- 'PB'
GSM5025879.seu[['case']] <- 7
GSM5025879.seu[['cond']] <- 'CLL'
GSM5025879.seu[['sorted']] <- 'whole'
GSM5025879.seu[['stage']] <- 'Ibrutinib_timepoint1'
GSM5025879.seu[['time']] <- 'Pretreatment_Ibrutinib'

# Save the object as .rds file
SaveSeuratRds(GSM5025879.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025879/GSM5025879.rds')
