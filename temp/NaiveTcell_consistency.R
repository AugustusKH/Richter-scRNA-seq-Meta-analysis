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

# Load the integrated object
integrated.seu <- readRDS("integrated.Rds")

# Level 2 cell type annotation
ref <- celldex::DatabaseImmuneCellExpressionData()
seurat_counts <- GetAssayData(merged.seu, layer = 'data')
pred <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.fine)

integrated.seu$singleR.fine_labels <- pred$labels[match(rownames(integrated.seu@meta.data), rownames(pred))]
DimPlot(integrated.seu, reduction = 'umap', group.by = 'singleR.fine_labels') # Visualised for level 2 cell type annotation

# Subset to Hing and Nadeu datasets
Hing_integrated.seu <- subset(x=integrated.seu, subset = ref=='Hing et al.')
Nadeu_integrated.seu <- subset(x=integrated.seu, subset = ref=='Nadeu et al.')

# Pseudobulk analysis using DESeq2 for Hing dataset in naive T cells
## Data aggregation and manipulation
cts <- AggregateExpression(Hing_integrated.seu, 
                           group.by = c("singleR.fine_labels", "orig_label"),
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

cts.split.modified$`T cells, CD4+, naive`[1:5,1:4] # Explore the first four rows and columns of the first member in the matrix list (B-cell expression matrix)

Hing_naiveCD4_express.mat <- as.data.frame(cts.split.modified$`T cells, CD4+, naive`)
Hing_naiveCD4_express.mat <- Hing_naiveCD4_express.mat %>%
  rename_with(~c('Hing_CLL1', 'Hing_CLL2', 'Hing_RT1', 'Hing_RT2'), c('CLL1', 'CLL2', 'RT1', 'RT2')) %>%
  rownames_to_column('genes')

## Extract the naive T cells expression matrix and metadata table
counts_naiveCD4cells <- cts.split.modified$`T cells, CD4+, naive`
colData <- data.frame(samples = colnames(counts_naiveCD4cells))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

## Create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_naiveCD4cells,
                              colData = colData,
                              design = ~ condition)

## Filter gene counts summation from the samples more than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")

## Run DESeq
dds <- DESeq(dds)

## Gene expression between RT and non-RT group
Hing_naiveCD4_res_RT_vs_CLL <- results(dds, contrast = c("condition", "Richter", "CLL"))
Hing_naiveCD4_res_RT_vs_CLL <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=Hing_naiveCD4_res_RT_vs_CLL, type = 'normal')
Hing_naiveCD4_EG <- as.data.frame(Hing_naiveCD4_res_RT_vs_CLL)

## DEG analysis
Hing_naiveCD4_DEG <- Hing_naiveCD4_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
Hing_naiveCD4_DEG

# Pseudobulk analysis using DESeq2 for Nadeu dataset in naive T cells --> cannot be done due to limitation of samples
## Data aggregation and manipulation
cts <- AggregateExpression(Nadeu_integrated.seu, 
                           group.by = c("singleR.fine_labels", "orig_label"),
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

cts.split.modified$`T cells, CD8+, naive`[1:5,1:2] # Explore the first four rows and columns of the first member in the matrix list (B-cell expression matrix)

Nadeu_naiveCD4_express.mat <- as.data.frame(cts.split.modified$`T cells, CD8+, naive`)
Nadeu_naiveCD4_express.mat <- Nadeu_naiveCD4_express.mat %>%
  rename_with(~c('Nadeu_CLL1', 'Nadeu_RT1'), c('CLL1', 'RT1')) %>%
  rownames_to_column('genes')

# DEGs analysis across Hing and Nadeu et al.
## Join expression matrix between Hing and Nadeu 
naiveCD4_express.mat <- Hing_naiveCD4_express.mat %>%
  inner_join(Nadeu_naiveCD4_express.mat, by=c(genes = 'genes'))

## Extract the B-cell expression matrix and metadata table
counts_naiveCD4cells <- naiveCD4_express.mat %>% 
  tibble::column_to_rownames(var = "genes") %>%
  dplyr::select(-c('Hing_CLL1', 'Hing_RT1'))
colData <- data.frame(samples = colnames(counts_naiveCD4cells))
colData <- colData %>%
  mutate(condition = ifelse(grepl('Hing', samples), 'Hing', 'Nadeu')) %>%
  column_to_rownames(var = "samples")
colData

## Create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_naiveCD4cells,
                              colData = colData,
                              design = ~ condition)

## Filter gene counts summation from the samples more than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Hing")

## Run DESeq
dds <- DESeq(dds)

## Gene expression between RT and non-RT group
naiveCD4_res_Hing_vs_Nadeu <- results(dds, contrast = c("condition", "Nadeu", "Hing"))
naiveCD4_res_Hing_vs_Nadeu <- lfcShrink(dds, contrast = c("condition", "Nadeu", "Hing"), res=naiveCD4_res_Hing_vs_Nadeu, type = 'normal')
naiveCD4_EG <- as.data.frame(naiveCD4_res_Hing_vs_Nadeu)

## DEG analysis
naiveCD4_DEG <- naiveCD4_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
naiveCD4_DEG



