# Consistency evaluation by Pseudo bulk DEGs analysis in PB from Hing et al. and Nadeu et al.

# Import the relevant packages
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(DESeq2)
library(SingleR)
library(celldex)
library(tidyverse)

# Set working directory
setwd("D:/Bioinfor PhD/Thesis/Progress lab meeting/Datasets/Measure_consistency")

# Import the Seurat datasets from Hing et al.
case1_CLL.seu <- readRDS('Hing_case1_CLL_PB.RDS')
case1_RT.seu <- readRDS('Hing_case1_RS_PB.Rds')
case2_CLL.seu <- readRDS('Hing_case2_CLL_PB.RDS')
case2_RT.seu <- readRDS('Hing_case2_RS_PB.Rds')

# Label origin
case1_CLL.seu[['orig_label']] <- 'CLL1'
case1_RT.seu[['orig_label']] <- 'RT1'
case2_CLL.seu[['orig_label']] <- 'CLL2'
case2_RT.seu[['orig_label']] <- 'RT2'

# QC for case 1 CLL
case1_CLL.seu[["percent.mt"]] <- PercentageFeatureSet(case1_CLL.seu, pattern = "^MT-")
VlnPlot(case1_CLL.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case1_CLL.seu <- subset(case1_CLL.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 1 RT
case1_RT.seu[["percent.mt"]] <- PercentageFeatureSet(case1_RT.seu, pattern = "^MT-")
VlnPlot(case1_RT.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case1_RT.seu <- subset(case1_RT.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 2 CLL
case2_CLL.seu[["percent.mt"]] <- PercentageFeatureSet(case2_CLL.seu, pattern = "^MT-")
VlnPlot(case2_CLL.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case2_CLL.seu <- subset(case2_CLL.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 2 RT
case2_RT.seu[["percent.mt"]] <- PercentageFeatureSet(case2_RT.seu, pattern = "^MT-")
VlnPlot(case2_RT.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case2_RT.seu <- subset(case2_RT.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# Merge case 1 datasets
Hing_PB_merged.seu <- merge(case1_CLL.seu, y=c(case1_RT.seu, case2_CLL.seu, case2_RT.seu))
Hing_PB_merged.seu <- JoinLayers(Hing_PB_merged.seu)

# Merged data pre-processing
Hing_PB_merged.seu <- NormalizeData(Hing_PB_merged.seu)
Hing_PB_merged.seu <- RunALRA(Hing_PB_merged.seu) # Impute the data
imputed_data <- GetAssayData(Hing_PB_merged.seu, layer = "data", assay = "alra") # Calculate imputed raw counts from imputed normalized data and add it to the count layer in the alra assay
count_like_data <- exp(imputed_data) - 1  # Subtract 1 to offset the addition of 1 during log transformation
count_like_data <- round(count_like_data) # Rounding
Hing_PB_merged.seu[['alra']] <- CreateAssayObject(counts = count_like_data)
Hing_PB_merged.seu <- FindVariableFeatures(object = Hing_PB_merged.seu)
Hing_PB_merged.seu <- ScaleData(object = Hing_PB_merged.seu)
Hing_PB_merged.seu <- RunPCA(object = Hing_PB_merged.seu)
Hing_PB_merged.seu <- RunUMAP(object = Hing_PB_merged.seu, dims = 1:20)

# Cell type annotation using SingleR
ref <- celldex::DatabaseImmuneCellExpressionData()
seurat_counts <- GetAssayData(Hing_PB_merged.seu, layer = 'data')

pred <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.main)

Hing_PB_merged.seu$singleR.labels <- pred$labels[match(rownames(Hing_PB_merged.seu@meta.data), rownames(pred))]
DimPlot(Hing_PB_merged.seu, reduction = 'umap', group.by = 'singleR.labels')

# Save the merged dataset
SaveSeuratRds(Hing_PB_merged.seu, file = "Hing_merged.Rds")
Hing_PB_merged.seu <- readRDS("Hing_merged.Rds")

# Pseudobulk analysis using DESeq2
## Data aggregation and manipulation
cts <- AggregateExpression(Hing_PB_merged.seu, 
                           group.by = c("singleR.labels", "orig_label"),
                           assays = "alra",
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$alra
cts.t <- t(cts) # Transpose the matrix
cts.t <- as.data.frame(cts.t)

splitRows <- gsub('_.*', '', rownames(cts.t)) # Substitute the word after _ to empty string for example 'B cells_CLL1' to 'B cells'
## Split the expression matrix to each cell and we will get a list of cell-type expression matrices
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows)) 

## This code is run as a loop to substitute string like 'B cells_CLL1' to 'CLL1' and transpose the split matrix
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub(".*_(.*)", "\\1", rownames(x))
  t(x)
})

cts.split.modified$`B cells`[1:5,1:4] # Explore the first four rows and columns of the first member in the matrix list (B-cell expression matrix)

## Check the percentage of zeros in the count matrix
num_zeros <- sum(cts.split.modified$`B cells` == 0)
matrix_dimensions <- dim(cts.split.modified$`B cells`)
zero_percent <- 100*num_zeros/(matrix_dimensions[1]*matrix_dimensions[2]); zero_percent

## Extract the B-cell expression matrix and metadata table
counts_Bcells <- cts.split.modified$`B cells`
colData <- data.frame(samples = colnames(counts_Bcells))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

## Create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_Bcells,
                              colData = colData,
                              design = ~ condition)

## Filter gene counts summation from the samples more than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")

## Run DESeq
dds <- DESeq(dds)

## Gene expression between RT and non-RT group
Hing_B_res_RT_vs_CLL <- results(dds, contrast = c("condition", "Richter", "CLL"))
Hing_B_res_RT_vs_CLL <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=Hing_B_res_RT_vs_CLL, type = 'normal')
Hing_B_EG <- as.data.frame(Hing_B_res_RT_vs_CLL)

## DEG analysis
Hing_B_DEG <- Hing_B_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
Hing_B_DEG

Hing_up_B_DEG <- Hing_B_DEG %>%
  filter(log2FoldChange > 1)
Hing_up_B_genes <- row.names(Hing_up_B_DEG)

Hing_down_B_DEG <- Hing_B_DEG %>%
  filter(log2FoldChange < -1)
Hing_down_B_genes <- row.names(Hing_down_B_DEG)

# Import the Seurat datasets from Nadeu et al.
case12_CLL.seu <- readRDS('case12_CLL_PB_PT3.RDS')
case12_RT.seu <- readRDS('case12_RT_PB.Rds')
case19_CLL.seu <- readRDS('case19_CLL_PB_PT3.RDS')
case19_RT.seu <- readRDS('case19_RT_PB.Rds')

# Label origin
case12_CLL.seu[['orig_label']] <- 'CLL1'
case12_RT.seu[['orig_label']] <- 'RT1'
case19_CLL.seu[['orig_label']] <- 'CLL2'
case19_RT.seu[['orig_label']] <- 'RT2'

# QC for case 1 CLL
case12_CLL.seu[["percent.mt"]] <- PercentageFeatureSet(case12_CLL.seu, pattern = "^MT-")
VlnPlot(case12_CLL.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case12_CLL.seu <- subset(case1_CLL.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 1 RT
case12_RT.seu[["percent.mt"]] <- PercentageFeatureSet(case12_RT.seu, pattern = "^MT-")
VlnPlot(case12_RT.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case12_RT.seu <- subset(case12_RT.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 2 CLL
case19_CLL.seu[["percent.mt"]] <- PercentageFeatureSet(case19_CLL.seu, pattern = "^MT-")
VlnPlot(case19_CLL.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case19_CLL.seu <- subset(case19_CLL.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# QC for case 2 RT
case19_RT.seu[["percent.mt"]] <- PercentageFeatureSet(case19_RT.seu, pattern = "^MT-")
VlnPlot(case19_RT.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
case19_RT.seu <- subset(case19_RT.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

# Merge case 1 datasets
Nadeu_PB_merged.seu <- merge(case12_CLL.seu, y=c(case12_RT.seu, case19_CLL.seu, case19_RT.seu))
Nadeu_PB_merged.seu <- JoinLayers(Nadeu_PB_merged.seu)

# Merged data pre-processing
Nadeu_PB_merged.seu <- NormalizeData(Nadeu_PB_merged.seu)
Nadeu_PB_merged.seu <- RunALRA(Nadeu_PB_merged.seu) # Impute the data
imputed_data <- GetAssayData(Nadeu_PB_merged.seu, layer = "data", assay = "alra") # Calculate imputed raw counts from imputed normalized data and add it to the count layer in the alra assay
count_like_data <- exp(imputed_data) - 1  # Subtract 1 to offset the addition of 1 during log transformation
count_like_data <- round(count_like_data) # Rounding
Nadeu_PB_merged.seu[['alra']] <- CreateAssayObject(counts = count_like_data)
Nadeu_PB_merged.seu <- FindVariableFeatures(object = Nadeu_PB_merged.seu)
Nadeu_PB_merged.seu <- ScaleData(object = Nadeu_PB_merged.seu)
Nadeu_PB_merged.seu <- RunPCA(object = Nadeu_PB_merged.seu)
Nadeu_PB_merged.seu <- RunUMAP(object = Nadeu_PB_merged.seu, dims = 1:20)

# Cell type annotation using SingleR
ref <- celldex::DatabaseImmuneCellExpressionData()
seurat_counts <- GetAssayData(Nadeu_PB_merged.seu, layer = 'data')

pred <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.main)

Nadeu_PB_merged.seu$singleR.labels <- pred$labels[match(rownames(Nadeu_PB_merged.seu@meta.data), rownames(pred))]
DimPlot(Nadeu_PB_merged.seu, reduction = 'umap', group.by = 'singleR.labels')

# Save the merged dataset
SaveSeuratRds(Nadeu_PB_merged.seu, file = "Nadeu_merged.Rds")
Nadeu_PB_merged.seu <- readRDS("Nadeu_merged.Rds")

# Pseudobulk analysis using DESeq2
## Data aggregation and manipulation
cts <- AggregateExpression(Nadeu_PB_merged.seu, 
                           group.by = c("singleR.labels", "orig_label"),
                           assays = "alra",
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$alra
cts.t <- t(cts) # Transpose the matrix
cts.t <- as.data.frame(cts.t)

splitRows <- gsub('_.*', '', rownames(cts.t)) # Substitute the word after _ to empty string for example 'B cells_CLL1' to 'B cells'
## Split the expression matrix to each cell and we will get a list of cell-type expression matrices
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows)) 

## This code is run as a loop to substitute string like 'B cells_CLL1' to 'CLL1' and transpose the split matrix
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub(".*_(.*)", "\\1", rownames(x))
  t(x)
})

cts.split.modified$`B cells`[1:5,1:4] # Explore the first four rows and columns of the first member in the matrix list (B-cell expression matrix)

## Check the percentage of zeros in the count matrix
num_zeros <- sum(cts.split.modified$`B cells` == 0)
matrix_dimensions <- dim(cts.split.modified$`B cells`)
zero_percent <- 100*num_zeros/(matrix_dimensions[1]*matrix_dimensions[2]); zero_percent

## Extract the B-cell expression matrix and metadata table
counts_Bcells <- cts.split.modified$`B cells`
colData <- data.frame(samples = colnames(counts_Bcells))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

## Create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_Bcells,
                              colData = colData,
                              design = ~ condition)

## Filter gene counts summation from the samples more than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")

## Run DESeq
dds <- DESeq(dds)

## Gene expression between RT and non-RT group
Nadeu_B_res_RT_vs_CLL <- results(dds, contrast = c("condition", "Richter", "CLL"))
Nadeu_B_res_RT_vs_CLL <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=Nadeu_B_res_RT_vs_CLL, type = 'normal')
Nadeu_B_EG <- as.data.frame(Nadeu_B_res_RT_vs_CLL)

## DEG analysis
Nadeu_B_DEG <- Nadeu_B_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
Nadeu_B_DEG

Nadeu_up_B_DEG <- Nadeu_B_DEG %>%
  filter(log2FoldChange > 1)
Nadeu_up_B_genes <- row.names(Nadeu_up_B_DEG)

Nadeu_down_B_DEG <- Nadeu_B_DEG %>%
  filter(log2FoldChange < -1)
Nadeu_down_B_genes <- row.names(Nadeu_down_B_DEG)

# Finding intersection between DEGs from Hing and Nadeu et al.
intersect(Hing_up_B_genes, Nadeu_up_B_genes)
intersect(Hing_down_B_genes, Nadeu_down_B_genes)

# Calculate correlation between log2FC of both datasets
Hing_B_EG$genes <- rownames(Hing_B_EG)
Nadeu_B_EG$genes <- rownames(Nadeu_B_EG)
joined_df <- Hing_B_EG %>%
  inner_join(Nadeu_B_EG, by=c(genes = 'genes')) %>%
  mutate(signFC.x = ifelse(log2FoldChange.x > 0, 1, -1)) %>%
  mutate(signFC.y = ifelse(log2FoldChange.y > 0, 1, -1))
cor(joined_df$log2FoldChange.x, joined_df$log2FoldChange.y, method = 'pearson')
sum(joined_df$signFC.x == joined_df$signFC.y)
nrow(joined_df)

# Spearman correlation of ranks
## Ranked by log2FoldChange
Hing_B_EG$rank <- rank(Hing_B_EG$log2FoldChange)
Nadeu_B_EG$rank <- rank(Nadeu_B_EG$log2FoldChange)
merged <- merge(Hing_B_EG, Nadeu_B_EG, by = "genes")
cor(merged$rank.x, merged$rank.y, method = "spearman")

## Ranked based on upregulated and downregulated DEGs
Hing_B_upEG <- Hing_B_EG[Hing_B_EG$log2FoldChange > 0,]
Hing_B_downEG <- Hing_B_EG[Hing_B_EG$log2FoldChange < 0,]
Nadeu_B_upEG <- Nadeu_B_EG[Nadeu_B_EG$log2FoldChange > 0,]
Nadeu_B_downEG <- Nadeu_B_EG[Nadeu_B_EG$log2FoldChange < 0,]

Hing_B_upEG$rank <- rank(Hing_B_upEG$padj)
Nadeu_B_upEG$rank <- rank(Nadeu_B_upEG$padj)
Hing_B_downEG$rank <- rank(Hing_B_downEG$padj)
Nadeu_B_downEG$rank <- rank(Nadeu_B_downEG$padj)
merged <- merge(Hing_B_upEG, Nadeu_B_upEG, by = "genes")
cor(merged$rank.x, merged$rank.y, method = "spearman")



