# Import the relevant packages
library(dplyr)
library(Seurat)

# Find the working path
getwd()

# Read 10x Genomics data
CMO306.data <- Read10X(data.dir = "/home/pakorns/hp-storage/scRNAseq/GSE201704/GSM6069857/CMO306/count_matrices/count/filtered_feature_bc_matrix")

# Explore the modalities inside the data
names(CMO306.data) # return: "Gene Expression"  "Antibody Capture"

# Create Seurat object
CMO306.seu <- CreateSeuratObject(counts = CMO306.data[["Gene Expression"]], min.features = 500, min.cells = 1)

# Explore ADT matrix
adt_counts <- CMO306.data[["Antibody Capture"]]
class(adt_counts) # Check the class
dim(adt_counts) # Dimensions (genes x barcodes)
## Show row names (feature names) and column names (barcodes)
rownames(adt_counts)  # all features (ADT tags) return: "CMO305" "CMO306" "CD3" "CD4.1" "CD8" "CD19.1" "CD11c" "CD56"
colnames(adt_counts)[1:5]  # First 5 barcodes

# Select ADT matrix containing only surface proteins
surface_rows <- grep("^CD", rownames(adt_count), value = TRUE)
adt_surface <- adt_count[surface_rows, ]

# Subset ADT matrix to only matching barcodes
common_barcodes <- colnames(CMO306.seu)  # RNA barcodes
adt_surface_filtered <- adt_surface[, common_barcodes]

# Create and add ADT assay to the object
adt_assay <- CreateAssayObject(counts = adt_surface_filtered)
CMO306.seu[["ADT"]] <- adt_assay

# QC
CMO306.seu[["percent.mt"]] <- PercentageFeatureSet(CMO306.seu, pattern = "^MT-")
CMO306.seu <- subset(CMO306.seu, subset = percent.mt < 20)

# Define other columns in meta-data
CMO306.seu[["HTO"]] <- "CMO306"
CMO306.seu[["ref"]] <- "Rejeski et al. (2022)"
CMO306.seu[["tissue"]] <- "PB"
CMO306.seu[["case"]] <- 1
CMO306.seu[["cond"]] <- "RT"
CMO306.seu[["sorted"]] <- "whole"

# Save the object as .rds file
SaveSeuratRds(CMO306.seu, file="/mnt/hp-storage/users/pakorns/scRNAseq/GSE201704/GSM6069857/CMO306/CMO306.rds")
