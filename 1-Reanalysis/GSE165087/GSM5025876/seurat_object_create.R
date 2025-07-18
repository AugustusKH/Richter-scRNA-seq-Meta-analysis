# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM5025876.data <- Read10X(data.dir = 'mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025876/count_matrices/outs/filtered_feature_bc_matrix')

# Create Seurat object
GSM5025876.seu <- CreateSeuratObject(counts = GSM5025876.data, min.cells = 3, min.features = 200)

# Define other columns in meta-data
GSM5025876.seu[['ref']] <- 'Penter et al. (2021)' 
GSM5025876.seu[['tissue']] <- 'PB'
GSM5025876.seu[['case']] <- 6
GSM5025876.seu[['cond']] <- 'CLL'
GSM5025876.seu[['sorted']] <- 'whole'
GSM5025876.seu[['stage']] <- 'RIC/Ibrutinib_relapse_timepoint2'
GSM5025876.seu[['time']] <- 'Relapse_post_RIC'

# Save the object as .rds file
SaveSeuratRds(GSM5025876.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087/GSM5025876/GSM5025876.rds')
