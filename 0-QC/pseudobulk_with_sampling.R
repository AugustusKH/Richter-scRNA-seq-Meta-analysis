# Pseudobulk analysis in the same patients but different conditions using random sampling

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

# Import the integrated object
integrated.seu <- readRDS("integrated.Rds")
DefaultAssay(integrated.seu)

# Subset B cells in each study
Hing_B_cell.seu <- subset(integrated.seu, subset = ref=='Hing et al.' & singleR.labels=='B cells')
Nadeu_B_cell.seu <- subset(integrated.seu, subset = ref=='Nadeu et al.' & singleR.labels=='B cells')

# Subset patients and conditions in each study
case1_CLL.seu <- subset(Hing_B_cell.seu, subset = case=='1' & cond=='CLL')
case1_RT.seu <- subset(Hing_B_cell.seu, subset = case=='1' & cond=='RS')
case2_CLL.seu <- subset(Hing_B_cell.seu, subset = case=='2' & cond=='CLL')
case2_RT.seu <- subset(Hing_B_cell.seu, subset = case=='2' & cond=='RS')
case12_CLL.seu <- subset(Nadeu_B_cell.seu, subset = case=='12' & cond=='CLL')
case12_RT.seu <- subset(Nadeu_B_cell.seu, subset = case=='12' & cond=='RT')
case19_CLL.seu <- subset(Nadeu_B_cell.seu, subset = case=='19' & cond=='CLL')
case19_RT.seu <- subset(Nadeu_B_cell.seu, subset = case=='19' & cond=='RT')

# Count the number of cells in each object
ncol(case1_CLL.seu) #250
ncol(case1_RT.seu) #7542
ncol(case2_CLL.seu) #98
ncol(case2_RT.seu) #1378
ncol(case12_CLL.seu) #2536
ncol(case12_RT.seu) #457
ncol(case19_CLL.seu) #983
ncol(case19_RT.seu) #1833

# Sampling cells from the patient objects
set.seed(111)

## Case 1 
### CLL
case1_CLL_all.cell <- colnames(case1_CLL.seu)
n_half <- floor(length(case1_CLL_all.cell) / 2)
cells_A <- sample(case1_CLL_all.cell, n_half)
cells_B <- setdiff(case1_CLL_all.cell, cells_A)
case1_CLL_A <- subset(case1_CLL.seu, cells = cells_A)
case1_CLL_B <- subset(case1_CLL.seu, cells = cells_B)

case1_CLL_A[['sampling']] <- 'case1_CLL_A'
case1_CLL_B[['sampling']] <- 'case1_CLL_B'

### RT
case1_RT_all.cell <- colnames(case1_RT.seu)
n_half <- floor(length(case1_RT_all.cell) / 2)
cells_A <- sample(case1_RT_all.cell, n_half)
cells_B <- setdiff(case1_RT_all.cell, cells_A)
case1_RT_A <- subset(case1_RT.seu, cells = cells_A)
case1_RT_B <- subset(case1_RT.seu, cells = cells_B)

case1_RT_A[['sampling']] <- 'case1_RT_A'
case1_RT_B[['sampling']] <- 'case1_RT_B'

## Case 2 
### CLL
case2_CLL_all.cell <- colnames(case2_CLL.seu)
n_half <- floor(length(case2_CLL_all.cell) / 2)
cells_A <- sample(case2_CLL_all.cell, n_half)
cells_B <- setdiff(case2_CLL_all.cell, cells_A)
case2_CLL_A <- subset(case2_CLL.seu, cells = cells_A)
case2_CLL_B <- subset(case2_CLL.seu, cells = cells_B)

case2_CLL_A[['sampling']] <- 'case2_CLL_A'
case2_CLL_B[['sampling']] <- 'case2_CLL_B'

### RT
case2_RT_all.cell <- colnames(case2_RT.seu)
n_half <- floor(length(case2_RT_all.cell) / 2)
cells_A <- sample(case2_RT_all.cell, n_half)
cells_B <- setdiff(case2_RT_all.cell, cells_A)
case2_RT_A <- subset(case2_RT.seu, cells = cells_A)
case2_RT_B <- subset(case2_RT.seu, cells = cells_B)

case2_RT_A[['sampling']] <- 'case2_RT_A'
case2_RT_B[['sampling']] <- 'case2_RT_B'

## Case 12 
### CLL
case12_CLL_all.cell <- colnames(case12_CLL.seu)
n_half <- floor(length(case12_CLL_all.cell) / 2)
cells_A <- sample(case12_CLL_all.cell, n_half)
cells_B <- setdiff(case12_CLL_all.cell, cells_A)
case12_CLL_A <- subset(case12_CLL.seu, cells = cells_A)
case12_CLL_B <- subset(case12_CLL.seu, cells = cells_B)

case12_CLL_A[['sampling']] <- 'case12_CLL_A'
case12_CLL_B[['sampling']] <- 'case12_CLL_B'

### RT
case12_RT_all.cell <- colnames(case12_RT.seu)
n_half <- floor(length(case12_RT_all.cell) / 2)
cells_A <- sample(case12_RT_all.cell, n_half)
cells_B <- setdiff(case12_RT_all.cell, cells_A)
case12_RT_A <- subset(case12_RT.seu, cells = cells_A)
case12_RT_B <- subset(case12_RT.seu, cells = cells_B)

case12_RT_A[['sampling']] <- 'case12_RT_A'
case12_RT_B[['sampling']] <- 'case12_RT_B'

## Case 19 
### CLL
case19_CLL_all.cell <- colnames(case19_CLL.seu)
n_half <- floor(length(case19_CLL_all.cell) / 2)
cells_A <- sample(case19_CLL_all.cell, n_half)
cells_B <- setdiff(case19_CLL_all.cell, cells_A)
case19_CLL_A <- subset(case19_CLL.seu, cells = cells_A)
case19_CLL_B <- subset(case19_CLL.seu, cells = cells_B)

case19_CLL_A[['sampling']] <- 'case19_CLL_A'
case19_CLL_B[['sampling']] <- 'case19_CLL_B'

### RT
case19_RT_all.cell <- colnames(case19_RT.seu)
n_half <- floor(length(case19_RT_all.cell) / 2)
cells_A <- sample(case19_RT_all.cell, n_half)
cells_B <- setdiff(case19_RT_all.cell, cells_A)
case19_RT_A <- subset(case19_RT.seu, cells = cells_A)
case19_RT_B <- subset(case19_RT.seu, cells = cells_B)

case19_RT_A[['sampling']] <- 'case19_RT_A'
case19_RT_B[['sampling']] <- 'case19_RT_B'

## Merge all sampling datasets
sampling_merged.seu <- merge(case1_CLL_A, y=c(case1_CLL_B, case1_RT_A, case1_RT_B, case2_CLL_A, case2_CLL_B, case2_RT_A, case2_RT_B,
                                              case12_CLL_A, case12_CLL_B, case12_RT_A, case12_RT_B, case19_CLL_A, case19_CLL_B, 
                                              case19_RT_A, case19_RT_B),
                             add.cell.ids = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'))

## Construct an expression matrix
cts <- AggregateExpression(sampling_merged.seu, 
                           group.by = c("sampling"),
                           assays = "alra",
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$alra
express_mat <- as.data.frame(cts)

## Perform DEG analysis
### Case 1
case1_express_mat <- express_mat %>%
  select(contains('case1-'))
colData <- data.frame(samples = colnames(case1_express_mat))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(countData = case1_express_mat,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")
dds <- DESeq(dds)

case1_res <- results(dds, contrast = c("condition", "Richter", "CLL"))
case1_res <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=case1_res, type = 'normal')
case1_EG <- as.data.frame(case1_res)

case1_DEG <- case1_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
case1_DEG

case1_upDEG <- case1_DEG %>%
  filter(log2FoldChange > 1)
case1_upgenes <- row.names(case1_upDEG)

case1_downDEG <- case1_DEG %>%
  filter(log2FoldChange < -1)
case1_downgenes <- row.names(case1_downDEG)

### Case 2
case2_express_mat <- express_mat %>%
  select(contains('case2-'))
colData <- data.frame(samples = colnames(case2_express_mat))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(countData = case2_express_mat,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")
dds <- DESeq(dds)

case2_res <- results(dds, contrast = c("condition", "Richter", "CLL"))
case2_res <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=case2_res, type = 'normal')
case2_EG <- as.data.frame(case2_res)

case2_DEG <- case2_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
case2_DEG

case2_upDEG <- case2_DEG %>%
  filter(log2FoldChange > 1)
case2_upgenes <- row.names(case2_upDEG)

case2_downDEG <- case2_DEG %>%
  filter(log2FoldChange < -1)
case2_downgenes <- row.names(case2_downDEG)

### Case 12
case12_express_mat <- express_mat %>%
  select(contains('case12-'))
colData <- data.frame(samples = colnames(case12_express_mat))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(countData = case12_express_mat,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")
dds <- DESeq(dds)

case12_res <- results(dds, contrast = c("condition", "Richter", "CLL"))
case12_res <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=case12_res, type = 'normal')
case12_EG <- as.data.frame(case12_res)

case12_DEG <- case12_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
case12_DEG

case12_upDEG <- case12_DEG %>%
  filter(log2FoldChange > 1)
case12_upgenes <- row.names(case12_upDEG)

case12_downDEG <- case12_DEG %>%
  filter(log2FoldChange < -1)
case12_downgenes <- row.names(case12_downDEG)

### Case 19
case19_express_mat <- express_mat %>%
  select(contains('case19-'))
colData <- data.frame(samples = colnames(case19_express_mat))
colData <- colData %>%
  mutate(condition = ifelse(grepl('RT', samples), 'Richter', 'CLL')) %>%
  column_to_rownames(var = "samples")

dds <- DESeqDataSetFromMatrix(countData = case19_express_mat,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "CLL")
dds <- DESeq(dds)

case19_res <- results(dds, contrast = c("condition", "Richter", "CLL"))
case19_res <- lfcShrink(dds, contrast = c("condition", "Richter", "CLL"), res=case19_res, type = 'normal')
case19_EG <- as.data.frame(case19_res)

case19_DEG <- case19_EG %>%
  filter(abs(log2FoldChange)>1 & padj<0.05)
case19_DEG

case19_upDEG <- case19_DEG %>%
  filter(log2FoldChange > 1)
case19_upgenes <- row.names(case19_upDEG)

case19_downDEG <- case19_DEG %>%
  filter(log2FoldChange < -1)
case19_downgenes <- row.names(case19_downDEG)

# Intersection analysis of upregualted DEGs
## Combine into a named list
deg_list <- list(
  Case1 = case1_upgenes,
  Case2 = case2_upgenes,
  Case12 = case12_upgenes,
  Case19 = case19_upgenes
)

# Intersection of all sets
upDEG_common_elements <- Reduce(intersect, deg_list)

## Compute intersection counts
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) length(intersect(deg_list[[x]], deg_list[[y]])))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

## Plot the heatmap
library(pheatmap)

pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         main = "Intersection of upregulated DEGs")

## Plot the Venn diagram
library(ggvenn)
ggvenn(
  deg_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + 
  ggtitle("Venn diagram of upregulated DEGs")

## Define Jaccard similarity
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

## Jaccard similarity
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) jaccard(deg_list[[x]], deg_list[[y]]))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "deeppink2"))(100),
         main = "Jaccard similarity of upregulated DEGs")

# Intersection analysis of downregualted DEGs
## Combine into a named list
deg_list <- list(
  Case1 = case1_downgenes,
  Case2 = case2_downgenes,
  Case12 = case12_downgenes,
  Case19 = case19_downgenes
)

# Intersection of all sets
downDEG_common_elements <- Reduce(intersect, deg_list)

## Compute intersection counts
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) length(intersect(deg_list[[x]], deg_list[[y]])))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

## Plot the heatmap
pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         main = "Intersection of downgulated DEGs")

## Plot the Venn diagram
ggvenn(
  deg_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + 
  ggtitle("Venn diagram of downregulated DEGs")

## Jaccard similarity
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) jaccard(deg_list[[x]], deg_list[[y]]))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "deeppink2"))(100),
         main = "Jaccard similarity of downregulated DEGs")

# Correlation analysis
library(purrr)
case1_EG <- tibble::rownames_to_column(case1_EG, var = "gene")
case2_EG <- tibble::rownames_to_column(case2_EG, var = "gene")
case12_EG <- tibble::rownames_to_column(case12_EG, var = "gene")
case19_EG <- tibble::rownames_to_column(case19_EG, var = "gene")

## List of data frames
df_list <- list(case1_EG, case2_EG, case12_EG, case19_EG)

## Perform inner joins across all
allcase_EG <- reduce(df_list, inner_join, by = "gene")
log2FC_df <- allcase_EG %>%
  select(log2FoldChange.x, log2FoldChange.y, log2FoldChange.x.x, log2FoldChange.y.y)

deg_list <- list(
  Case1 = log2FC_df$log2FoldChange.x,
  Case2 = log2FC_df$log2FoldChange.y,
  Case12 = log2FC_df$log2FoldChange.x.x,
  Case19 = log2FC_df$log2FoldChange.y.y
)

## Compute Pearson correlation
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) cor(deg_list[[x]], deg_list[[y]], method = 'pearson'))))

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

## Plot the heatmap
pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Pearson correlation of log2FC")

## Compute Spearman correlation
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) cor(deg_list[[x]], deg_list[[y]], method = 'spearman'))))

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

## Plot the heatmap
pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Spearman correlation of log2FC")

# Gene Ranking analysis
## Upregulated DEGs
case1_upDEG <- tibble::rownames_to_column(case1_upDEG, var = "gene")
case2_upDEG <- tibble::rownames_to_column(case2_upDEG, var = "gene")
case12_upDEG <- tibble::rownames_to_column(case12_upDEG, var = "gene")
case19_upDEG <- tibble::rownames_to_column(case19_upDEG, var = "gene")

case1_upDEG$rank <- rank(case1_upDEG$padj)
case2_upDEG$rank <- rank(case2_upDEG$padj)
case12_upDEG$rank <- rank(case12_upDEG$padj)
case19_upDEG$rank <- rank(case19_upDEG$padj)

df_list <- list(case1_upDEG, case2_upDEG, case12_upDEG, case19_upDEG)
upDEG_merged <- reduce(df_list, inner_join, by = "gene")

### Compute Spearman correlation
deg_list <- list(
  Case1 = upDEG_merged$rank.x,
  Case2 = upDEG_merged$rank.y,
  Case12 = upDEG_merged$rank.x.x,
  Case19 = upDEG_merged$rank.y.y
)

inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) cor(deg_list[[x]], deg_list[[y]], method = 'spearman')))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

### Plot the heatmap
pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "green"))(100),
         main = "Spearman correlation of rank padj upregulated DEGs")

## Downregulated DEGs
case1_downDEG <- tibble::rownames_to_column(case1_downDEG, var = "gene")
case2_downDEG <- tibble::rownames_to_column(case2_downDEG, var = "gene")
case12_downDEG <- tibble::rownames_to_column(case12_downDEG, var = "gene")
case19_downDEG <- tibble::rownames_to_column(case19_downDEG, var = "gene")

case1_downDEG$rank <- rank(case1_downDEG$padj)
case2_downDEG$rank <- rank(case2_downDEG$padj)
case12_downDEG$rank <- rank(case12_downDEG$padj)
case19_downDEG$rank <- rank(case19_downDEG$padj)

df_list <- list(case1_downDEG, case2_downDEG, case12_downDEG, case19_downDEG)
downDEG_merged <- reduce(df_list, inner_join, by = "gene")

### Compute Spearman correlation
deg_list <- list(
  Case1 = downDEG_merged$rank.x,
  Case2 = downDEG_merged$rank.y,
  Case12 = downDEG_merged$rank.x.x,
  Case19 = downDEG_merged$rank.y.y
)

inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) cor(deg_list[[x]], deg_list[[y]], method = 'spearman')))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

### Plot the heatmap
pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "green"))(100),
         main = "Spearman correlation of rank padj downregulated DEGs")

# GSEA comparison
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)

## Hallmark set
set.seed(2020)
hs_hm_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  collection = "H"
)

## Case 1 
lfc_vector <- case1_EG$log2FoldChange # Create a named vector ranked based on the log2 fold change
names(lfc_vector) <- case1_EG$gene
lfc_vector <- na.omit(lfc_vector) # Omit any NA values
lfc_vector <- sort(lfc_vector, decreasing = TRUE) # Sort the log2 fold change in descending order

HM_gse <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hm_sets,
    gs_name,
    gene_symbol
  )
)

case1_GSEA.result <- HM_gse@result

## Case 2 
lfc_vector <- case2_EG$log2FoldChange # Create a named vector ranked based on the log2 fold change
names(lfc_vector) <- case2_EG$gene
lfc_vector <- na.omit(lfc_vector) # Omit any NA values
lfc_vector <- sort(lfc_vector, decreasing = TRUE) # Sort the log2 fold change in descending order

HM_gse <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hm_sets,
    gs_name,
    gene_symbol
  )
)

case2_GSEA.result <- HM_gse@result

## Case 12 
lfc_vector <- case12_EG$log2FoldChange # Create a named vector ranked based on the log2 fold change
names(lfc_vector) <- case12_EG$gene
lfc_vector <- na.omit(lfc_vector) # Omit any NA values
lfc_vector <- sort(lfc_vector, decreasing = TRUE) # Sort the log2 fold change in descending order

HM_gse <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hm_sets,
    gs_name,
    gene_symbol
  )
)

case12_GSEA.result <- HM_gse@result

## Case 19 
lfc_vector <- case19_EG$log2FoldChange # Create a named vector ranked based on the log2 fold change
names(lfc_vector) <- case19_EG$gene
lfc_vector <- na.omit(lfc_vector) # Omit any NA values
lfc_vector <- sort(lfc_vector, decreasing = TRUE) # Sort the log2 fold change in descending order

HM_gse <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    hs_hm_sets,
    gs_name,
    gene_symbol
  )
)

case19_GSEA.result <- HM_gse@result

## Intersection analysis of significant functional enrichment
deg_list <- list(
  Case1 = case1_GSEA.result$ID,
  Case2 = case2_GSEA.result$ID,
  Case12 = case12_GSEA.result$ID,
  Case19 = case19_GSEA.result$ID
)

# Intersection of all sets
GSEA_common_elements <- Reduce(intersect, deg_list)

inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) length(intersect(deg_list[[x]], deg_list[[y]])))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "deepskyblue2"))(100),
         main = "Intersection of significant enrichment terms")

## Plot the Venn diagram
ggvenn(
  deg_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + 
  ggtitle("Venn diagram of GSEA")

## Jaccard similarity
inter_matrix <- outer(
  names(deg_list),
  names(deg_list),
  Vectorize(function(x, y) jaccard(deg_list[[x]], deg_list[[y]]))
)

dimnames(inter_matrix) <- list(names(deg_list), names(deg_list))

pheatmap(inter_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c("white", "deepskyblue2"))(100),
         main = "Jaccard similarity of enrichment terms")
