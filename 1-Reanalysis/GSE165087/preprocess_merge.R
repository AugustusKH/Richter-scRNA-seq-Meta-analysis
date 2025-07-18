# Import the relevant package
library(Seurat)
library(SingleR)
library(celldex)

# Set working path
setwd('/mnt/hp-storage/users/pakorns/scRNAseq/GSE165087')

# Read the merge object
merged.seu <- readRDS('GSE165087_QC_merged.rds')

# Join layers in the merged object
merged.seu <- JoinLayers(merged.seu)

# Data pre-processing
merged.seu <- NormalizeData(merged.seu)
merged.seu <- FindVariableFeatures(object = merged.seu)
merged.seu <- ScaleData(object = merged.seu)
merged.seu <- RunPCA(object = merged.seu)
merged.seu <- RunUMAP(object = merged.seu, dims = 1:20)

# Cell type annotation
## Primary annotation
ref <- readRDS('DICE_ref.rds')
seurat_counts <- GetAssayData(merged.seu, layer = 'data')

pred1 <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.main)

merged.seu$singleR.labels_main <- pred1$labels[match(rownames(merged.seu@meta.data), rownames(pred1))]

## Secondary annotation
pred2 <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.fine)

merged.seu$singleR.labels_fine <- pred2$labels[match(rownames(merged.seu@meta.data), rownames(pred2))]

# Save the object
saveRDS(merged.seu, file = 'merged_annotated.rds')
