# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
GSM6069857.data <- Read10X(data.dir = "/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/all/count_matrices/count/filtered_feature_bc_matrix")

# Explore the modalities inside the data
names(GSM6069857.data) # return: "Gene Expression"  "Multiplexing Capture"

# Create Seurat object
GSM6069857.seu <- CreateSeuratObject(counts = GSM6069857.data[["Gene Expression"]], min.features = 500, min.cells = 1)

# Explore HTO matrix
hto_counts <- GSM6069857.data[["Multiplexing Capture"]]
class(hto_counts) # Check the class
dim(hto_counts) # Dimensions (genes x barcodes)
## Show row names (feature names) and column names (barcodes)
rownames(hto_counts)[1:2]  # First 5 features (HTO tags) return: CMO305 and CMO306
colnames(hto_counts)[1:5]  # First 5 barcodes

# Subset HTO matrix to only matching barcodes
common_barcodes <- colnames(GSM6069857.seu)  # RNA barcodes
hto_counts_filtered <- hto_counts[, common_barcodes]

# Create and add HTO assay
GSM6069857.seu[["HTO"]] <- CreateAssayObject(counts = hto_counts_filtered)

# Perform HTO demultiplexing
GSM6069857.seu <- NormalizeData(GSM6069857.seu, assay = "HTO", normalization.method = "CLR")
GSM6069857.seu <- HTODemux(GSM6069857.seu, assay = "HTO", positive.quantile = 0.99)

# New metadata columns:
# - classification: GSM6069857.seu$HTO_classification.global
# - assigned tag: GSM6069857.seu$HTO_maxID

# Example: Rename the classification column if desired
GSM6069857.seu$HTO_tag <- GSM6069857.seu$HTO_maxID

# Define other columns in meta-data
GSM6069857.seu[['tissue']] <- 'PB'
GSM6069857.seu[['cond']] <- 'RT'
GSM6069857.seu[['case']] <- 1
GSM6069857.seu[['sorted']] <- 'whole'
GSM6069857.seu[['ref']] <- 'Rejeski et al. (2022)'

# Save the object as .rds file
SaveSeuratRds(GSM6069857.seu, file = '/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/all/GSM6069857.rds')




