# Pseudo bulk DEGs analysis in PB from Hing et al. and Nadeu et al. (merging and integration before analysis)

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

# Merge all datasets
Hing_merged.seu <- merge(case1_CLL.seu, y=c(case1_RT.seu, case2_CLL.seu, case2_RT.seu),
                         add.cell.ids = c('Hing_CLL1', 'Hing_RT1', 'Hing_CLL2', 'Hing_RT2'))
Nadeu_merged.seu <- merge(case12_CLL.seu, y=c(case12_RT.seu, case19_CLL.seu, case19_RT.seu),
                          add.cell.ids = c('Nadeu_CLL12', 'Nadeu_RT12', 'Nadeu_CLL19', 'Nadeu_RT19'))
merged.seu <- merge(Hing_merged.seu, y=c(Nadeu_merged.seu))
merged.seu <- JoinLayers(merged.seu)

# QC the merged dataset
merged.seu[['orig.ident']] <- merged.seu[['case']] # Set the identity class to orig.ident
Idents(merged.seu) <- merged.seu@meta.data$orig.ident
merged.seu[["percent.mt"]] <- PercentageFeatureSet(merged.seu, pattern = "^MT-")
VlnPlot(merged.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # plot for QC metrices
merged.seu <- subset(merged.seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# Merged data pre-processing
merged.seu <- NormalizeData(merged.seu)
merged.seu <- RunALRA(merged.seu) # Impute the data
imputed_data <- GetAssayData(merged.seu, layer = "data", assay = "alra") # Calculate imputed raw counts from imputed normalized data and add it to the count layer in the alra assay
count_like_data <- exp(imputed_data) - 1  # Subtract 1 to offset the addition of 1 during log transformation
count_like_data <- round(count_like_data) # Rounding
merged.seu[['alra']] <- CreateAssayObject(counts = count_like_data)
merged.seu <- FindVariableFeatures(object = merged.seu)
merged.seu <- ScaleData(object = merged.seu)
merged.seu <- RunPCA(object = merged.seu)
merged.seu <- RunUMAP(object = merged.seu, dims = 1:20)

# Cell type annotation using SingleR
ref <- celldex::DatabaseImmuneCellExpressionData()
seurat_counts <- GetAssayData(merged.seu, layer = 'data')

pred <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.main)

merged.seu$singleR.labels <- pred$labels[match(rownames(merged.seu@meta.data), rownames(pred))]
DimPlot(merged.seu, reduction = 'umap', group.by = 'singleR.labels')

# Save the merged file as .RDS
SaveSeuratRds(merged.seu, file = "merged.Rds")
merged.seu <- readRDS("merged.Rds")

# Integrate the merged dataset
obj.list <- SplitObject(merged.seu, split.by = "orig.ident")

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated.seu <- IntegrateData(anchorset = anchors)

## Downstreaming analysis of the integrated dataset
integrated.seu <- ScaleData(object = integrated.seu)
integrated.seu <- RunPCA(object = integrated.seu)
integrated.seu <- RunUMAP(object = integrated.seu, dims = 1:20)
integrated.seu <- FindNeighbors(integrated.seu, dims = 1:20)
integrated.seu <- FindClusters(integrated.seu, resolution = 0.05)

seurat_counts <- GetAssayData(integrated.seu, layer = 'data')

pred <- SingleR(test = seurat_counts,
                ref = ref,
                labels = ref$label.main)

integrated.seu$singleR.labels <- pred$labels[match(rownames(integrated.seu@meta.data), rownames(pred))]
DimPlot(integrated.seu, reduction = 'umap', group.by = 'singleR.labels')

# Save the integrated file as .RDS
SaveSeuratRds(integrated.seu, file = "integrated.Rds")
integrated.seu <- readRDS("integrated.Rds")

# Subset to Hing and Nadeu datasets
Hing_integrated.seu <- subset(x=integrated.seu, subset = ref=='Hing et al.')
Nadeu_integrated.seu <- subset(x=integrated.seu, subset = ref=='Nadeu et al.')

# Pseudobulk analysis using DESeq2 for Hing dataset
## Data aggregation and manipulation
cts <- AggregateExpression(Hing_integrated.seu, 
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

Hing_B_express.mat <- as.data.frame(cts.split.modified$`B cells`)
Hing_B_express.mat <- Hing_B_express.mat %>%
  rename_with(~c('Hing_CLL1', 'Hing_CLL2', 'Hing_RT1', 'Hing_RT2'), c('CLL1', 'CLL2', 'RT1', 'RT2')) %>%
  rownames_to_column('genes')

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

# Pseudobulk analysis using DESeq2 for Nadeu dataset
## Data aggregation and manipulation
cts <- AggregateExpression(Nadeu_integrated.seu, 
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

Nadeu_B_express.mat <- as.data.frame(cts.split.modified$`B cells`)
Nadeu_B_express.mat <- Nadeu_B_express.mat %>%
  rename_with(~c('Nadeu_CLL1', 'Nadeu_CLL2', 'Nadeu_RT1', 'Nadeu_RT2'), c('CLL1', 'CLL2', 'RT1', 'RT2')) %>%
  rownames_to_column('genes')

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
#Hing_B_EG$rank <- rank(Hing_B_EG$log2FoldChange)
#Nadeu_B_EG$rank <- rank(Nadeu_B_EG$log2FoldChange)
Hing_B_EG$rank <- rank(Hing_B_EG$padj)
Nadeu_B_EG$rank <- rank(Nadeu_B_EG$padj)
merged <- merge(Hing_B_EG, Nadeu_B_EG, by = "genes")
cor(merged$rank.x, merged$rank.y, method = "spearman")

## Ranked by up- or downregulated genes  
Hing_B_upEG <- Hing_B_EG[Hing_B_EG$log2FoldChange > 0,]
Hing_B_downEG <- Hing_B_EG[Hing_B_EG$log2FoldChange < 0,]
Nadeu_B_upEG <- Nadeu_B_EG[Nadeu_B_EG$log2FoldChange > 0,]
Nadeu_B_downEG <- Nadeu_B_EG[Nadeu_B_EG$log2FoldChange < 0,]

Hing_B_downEG$rank <- rank(Hing_B_downEG$padj)
Nadeu_B_downEG$rank <- rank(Nadeu_B_downEG$padj)
merged <- merge(Hing_B_downEG, Nadeu_B_downEG, by = "genes")
cor(merged$rank.x, merged$rank.y, method = "spearman")

# Join expression matrix between Hing and Nadeu 
B_express.mat <- Hing_B_express.mat %>%
  inner_join(Nadeu_B_express.mat, by=c(genes = 'genes'))
cor(B_express.mat$Hing_CLL1, B_express.mat$Nadeu_CLL2, method='spearman')

# PCA to identify batch effects
library(factoextra)
B_express_row.mat <- B_express.mat %>% 
  remove_rownames %>% 
  column_to_rownames(var="genes")
t_B_express.mat <- t(B_express_row.mat)
zero_var_genes <- apply(t_B_express.mat, 2, function(x) var(x) == 0)
#sum(zero_var_genes)
t_B_express.mat <- t_B_express.mat[, !zero_var_genes]
res.pca <- prcomp(t_B_express.mat, scale = TRUE)

fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Batch correction using Combat-Seq
library(sva)
batch <- c(0, 0, 0, 0, 1, 1, 1, 1)
group <- c(0, 0, 1, 1, 0, 0, 1, 1)
adjusted <- ComBat_seq(as.matrix(B_express_row.mat), batch=batch, group=group)

t_adjusted <- t(adjusted)
zero_var_genes <- apply(t_adjusted, 2, function(x) var(x) == 0)
#sum(zero_var_genes)
t_adjusted <- t_adjusted[, !zero_var_genes]
res.pca <- prcomp(t_adjusted, scale = TRUE)

fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
