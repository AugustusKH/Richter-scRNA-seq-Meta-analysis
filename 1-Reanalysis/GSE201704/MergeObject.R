# Import the package
library(Seurat)

# Read the rds files
CMO305.seu <- readRDS('/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/CMO305/CMO305.rds')
CMO306.seu <- readRDS('/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/CMO306/CMO306.rds')

# Merge the objects
GSM6069857_merged.seu <- merge(CMO305.seu, y = CMO306.seu, add.cell.ids = c('CMO305', 'CMO306'))

# Save merged object
saveRDS(GSM6069857_merged.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/GSM6069857_merged.rds')
