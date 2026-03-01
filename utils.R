# Functions used in this study

## Import the relevant library
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
library(DT)
library(glue)
library(matrixStats)
library(ggsci)

## Load the data for enrichment analysis
hs_hm_sets <- read.csv(here::here("list_enrichment/hs_hm_sets.csv"))
hs_go_sets <- read.csv(here::here("list_enrichment/hs_go_sets.csv"))
hs_kegg_sets <- read.csv(here::here("list_enrichment/hs_kegg_sets.csv"))
hs_reactome_sets <- read.csv(here::here("list_enrichment/hs_reactome_sets.csv"))
hs_wikipath_sets <- read.csv(here::here("list_enrichment/hs_wikipath_sets.csv"))

## Define a function to filter out doublets based DoubletFinder package
doublet_filter <- function(seurat.obj, dims=1:30, sct=TRUE){
  # Set the default layer 
  if(sct){
    assay_to_use <- "SCT"
  } else {
    assay_to_use <- "RNA"
  }
  DefaultAssay(seurat.obj) <- assay_to_use
  
  # Create 10X multiple rate table (https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled)
  multiple_rate.tbl <- data.frame(
    'Multiplet_rate'= c(0.004, 0.008, 0.016, 0.024, 0.032, 0.04, 0.048, 0.056, 0.064, 0.072, 0.08),
    'Loaded_cells' = c(825, 1650, 3300, 4950, 6600, 8250, 9900, 11550, 13200, 14850, 16500),
    'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  
  # Calculate rate for recovery
  multiplet_rate <- multiple_rate.tbl %>% 
    dplyr::filter(Recovered_cells < nrow(seurat.obj@meta.data)) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
    dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) 
  
  # pK Identification (no ground-truth)
  sweep_list <- paramSweep(seurat.obj, PCs = dims, sct = sct) 
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- seurat.obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi <- round(multiplet_rate * nrow(seurat.obj@meta.data))
  nExp.poi.adj <- round(nExp.poi*(1-homotypic.prop))
  
  # Run DoubletFinder
  seurat.obj <- doubletFinder(seu = seurat.obj, PCs = dims, pK = optimal.pk, nExp = nExp.poi.adj, sct = sct) 
  
  # Change name of metadata column with Singlet/Doublet information
  colnames(seurat.obj@meta.data)[grepl('DF.classifications.*', colnames(seurat.obj@meta.data))] <- "doublet_finder"
  
  # Filter out doublets
  singlet <- with(seurat.obj@meta.data, doublet_finder == "Singlet")
  
  if (any(singlet)) {
    # Subset only if at least one cell passes
    seurat.obj <- subset(seurat.obj, cells = colnames(seurat.obj)[singlet])
  } 
  return(seurat.obj)
} 

## Define a function for cell type annotation using SingleR

### Define a function for primary cell type annotation using SingleR
singleR_primary_annotation <- function(seurat.obj, ref){
  seurat_counts <- GetAssayData(seurat.obj, layer = 'data')
  pred <- SingleR(test = seurat_counts,
                  ref = ref,
                  labels = ref$label.main)
  seurat.obj$singleR.labels_main <- pred$labels[match(rownames(seurat.obj@meta.data), rownames(pred))]
  return(seurat.obj)
}

### Define a function for secondary cell type annotation using SingleR
singleR_secondary_annotation <- function(seurat.obj, ref){
  seurat_counts <- GetAssayData(seurat.obj, layer = 'data')
  pred <- SingleR(test = seurat_counts,
                  ref = ref,
                  labels = ref$label.fine)
  seurat.obj$singleR.labels_fine <- pred$labels[match(rownames(seurat.obj@meta.data), rownames(pred))]
  return(seurat.obj)
}

### Define a main function for annotation
singleR_cell_annotation <- function(seurat.obj, annotation="primary", ref = NULL){
  if (is.null(ref)) {
    ref <- celldex::DatabaseImmuneCellExpressionData()
  }
  
  if(annotation == "primary"){
    seurat.obj <- singleR_primary_annotation(seurat.obj, ref)
  } else if (annotation == "secondary"){
    seurat.obj <- singleR_secondary_annotation(seurat.obj, ref)
  } else if (annotation == "both"){
    seurat.obj <- singleR_primary_annotation(seurat.obj, ref)
    seurat.obj <- singleR_secondary_annotation(seurat.obj, ref)
  } else {
    stop("annotation must be one of 'primary', 'secondary', or 'both'")
  }
  
  return(seurat.obj)  
}

## Define a function to calculate MAD cutoff for cell quality
get_mad_cutoffs <- function(x, k = 3) {
  med <- median(x, na.rm = TRUE)
  madv <- mad(x, constant = 1, na.rm = TRUE)  # default constant 1 gives median absolute deviation
  lower <- med - k * madv
  upper <- med + k * madv
  return(c(lower = lower, upper = upper))
}

## Define a function for seurat cell quality control
seurat_quality_control <- function(seurat.obj, 
                                   species="human", 
                                   min_nCount = 600, max_nCount = 15000,
                                   min_nFeature = 200, max_nFeature = 4000,
                                   max_percent_mt = 10, 
                                   min_percent_ribo = 10, max_percent_ribo = 60,
                                   mad_cutoff=FALSE){
  # Add QC metrics
  if (species == "human"){
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern="^MT-.")
    seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj, pattern="^RP[SL]")
  } else if (species == "mouse") {
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
    seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj, pattern = "^Rp[sl]")
  } else {
    message("Add species between human or mouse!")
  }
  
  # Create logical vector of cells passing QC
  if(mad_cutoff){
    # Compute thresholds
    nFeature <- FetchData(seurat.obj, "nFeature_RNA")$nFeature_RNA
    nCount <- FetchData(seurat.obj, "nCount_RNA")$nCount_RNA
    
    # Median ± 3*MAD for upper bound
    max_nFeature <- ceiling(get_mad_cutoffs(nFeature, k = 3)["upper"])
    max_nCount <- ceiling(get_mad_cutoffs(nCount, k = 3)["upper"])
  }
  
  pass_qc <- with(seurat.obj@meta.data,
                  nCount_RNA > min_nCount & nCount_RNA < max_nCount &
                    nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature &
                    percent.ribo > min_percent_ribo & percent.ribo < max_percent_ribo &
                    percent.mt < max_percent_mt)
  
  n_pass <- sum(pass_qc)
  
  if (n_pass > 200) {
    seurat.obj <- subset(seurat.obj, cells = colnames(seurat.obj)[pass_qc])
  } else {
    warning("Few or no cells passed QC — returning empty object")
    return(NULL) # return empty object instead:
  }
  return(seurat.obj)
}

## Define a function for seurat pre-processing 
seurat_preprocess <- function(seurat.obj, dims = 1:30){
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)
  seurat.obj <- FindNeighbors(seurat.obj, dims = dims)
  seurat.obj <- FindClusters(seurat.obj, resolution=0.25)
  seurat.obj <- RunUMAP(seurat.obj, dims = dims)
  return(seurat.obj)
}

## Define a function for seurat pre-processing with SCTransform
seurat_SCTransform <- function(seurat.obj, dims = 1:30){
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- SCTransform(seurat.obj, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
  seurat.obj <- RunPCA(seurat.obj)
  seurat.obj <- FindNeighbors(seurat.obj, dims = dims)
  seurat.obj <- FindClusters(seurat.obj, resolution=0.25)
  seurat.obj <- RunUMAP(seurat.obj, dims = dims)
  return(seurat.obj)
}

## Define a function for cell clustering by varying resolutions
plot_varying_cluster_resolution <- function(seurat.obj, assay_to_use = "SCT") {
  # Pick the right prefix depending on assay
  prefix <- ifelse(tolower(assay_to_use) == "sct", "SCT_snn_res.", "RNA_snn_res.")
  
  # Check which clustering columns exist
  res.cols <- grep(paste0("^", prefix), colnames(seurat.obj@meta.data), value = TRUE)
  
  if (length(res.cols) == 0) {
    stop("No clustering results found for assay: ", assay_to_use)
  }
  
  # Generate one UMAP per resolution, with labels and no legend
  plots <- lapply(seq(0.1, 0.9, by = 0.1), function(res) {
    DimPlot(
      seurat.obj,
      reduction = "umap",
      group.by = str_c(prefix, res),
      label = TRUE,
      repel = TRUE,
      label.size = 4
    ) + NoLegend() 
  })
  
  # Combine all plots into a grid
  clusters.plt <- wrap_plots(plots, ncol = 3)
  
  return(clusters.plt)
}

## Define a function for pseudobulk analysis using DESeq2 package
pseudobulk_deseq2 <- function(seurat.obj, aggregate_group.by = "case", RT_subjects, padj_NA = FALSE){
  # Data aggregation
  cts <- AggregateExpression(seurat.obj,
                             group.by = aggregate_group.by,
                             assays = "RNA",
                             slot = "counts",
                             resturn.seurat = FALSE) 
  cts <- cts$RNA
  
  # Create pseudocount matrix for tumour cells
  pseudo_count.mat <- as.data.frame(cts)
  
  # Check RT_subjects are present
  if (!all(RT_subjects %in% colnames(pseudo_count.mat))){
    message("Some RT_subject(s) not found in column names. Please check input.")
  } 
  
  # Create metadata for DESeq2 object
  colData <- data.frame(samples = colnames(pseudo_count.mat))
  colData$condition <- ifelse(colData$samples %in% RT_subjects, "RT", "CLL")
  rownames(colData) <- colData$samples
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = pseudo_count.mat,
                                colData = colData,
                                design = ~ condition)
  
  # Keep genes with at least 10 counts across all samples
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Run DESeq2
  dds$condition <- relevel(dds$condition, ref = "CLL")
  dds <- DESeq(dds)
  
  # Generate a table for gene expression
  res <- results(dds,
                 contrast = c("condition", "RT", "CLL"),
                 independentFiltering = padj_NA)
  res <- lfcShrink(dds, 
                   contrast = c("condition", "RT", "CLL"), 
                   res = res,
                   type = "normal")
  
  return(as.data.frame(res))
}

## Define a function for volcano plot
volcano_plot <- function(dea.df, 
                         FC=1, 
                         p=0.05, 
                         x_col="log2FoldChange",
                         y_col="padj",
                         xlim=c(-5,5), 
                         ylim=c(-0.1,12), 
                         boxedLabels = FALSE,
                         genes_of_interest=NULL){
  
  # Create custom colours
  keyvals <- rep('grey75', nrow(dea.df))
  names(keyvals) <- rep('NS', nrow(dea.df))
  
  keyvals[which(abs(dea.df[[x_col]]) > FC & dea.df[[y_col]] > p)] <- 'grey50'
  names(keyvals)[which(abs(dea.df[[x_col]]) > FC & dea.df[[y_col]] > p)] <- 'log2FoldChange'
  
  keyvals[which(abs(dea.df[[x_col]]) < FC & dea.df[[y_col]] < p)] <- 'grey25'
  names(keyvals)[which(abs(dea.df[[x_col]])  < FC & dea.df[[y_col]] < p)] <- '-Log10Q'
  
  keyvals[which(dea.df[[x_col]] < -FC & dea.df[[y_col]] < p)] <- '#0000CD'
  names(keyvals)[which(dea.df[[x_col]]  < -FC & dea.df[[y_col]] < p)] <- 'Signif. down-regulated'
  
  keyvals[which(dea.df[[x_col]] > FC & dea.df[[y_col]] < p)] <- '#DC143C'
  names(keyvals)[which(dea.df[[x_col]] > FC & dea.df[[y_col]] < p)] <- 'Signif. up-regulated'
  
  # Plot
  volcano.plt <- EnhancedVolcano(
    dea.df,
    lab = rownames(dea.df),
    x = x_col,
    y = y_col,
    selectLab = genes_of_interest,   # only label genes of interest
    xlim = xlim,
    ylim = ylim,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10] ~ italic(P)),
    title = 'Custom colour over-ride',
    pCutoff = p,
    FCcutoff = FC,
    max.overlaps = Inf, # show all overlapping labels
    boxedLabels = boxedLabels,  # box the labels so they are legible
    parseLabels = FALSE,
    pointSize = 1,
    labSize = 4,
    colCustom = keyvals,
    colAlpha = 0.3,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    colConnectors = 'grey50',
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black'
  )
  
  return(volcano.plt)
}
  
## Define a function for GSEA in seurat
HM_GSEA <- function(lfc_vector, gene_ident="gene_symbol", minGSSize = 25, maxGSSize = 500){
  set.seed(2020)
  hm_gsea <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = minGSSize, # Minimum gene set size
    maxGSSize = maxGSSize, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      hs_hm_sets,
      gs_name,
      gene_ident
    )
  )
  return(hm_gsea)
} 

GO_GSEA <- function(lfc_vector, gene_ident="gene_symbol", minGSSize = 25, maxGSSize = 500){
  set.seed(2020)
  gobp_gsea <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = minGSSize, # Minimum gene set size
    maxGSSize = maxGSSize, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      hs_go_sets,
      gs_name,
      gene_ident
    )
  )
  return(gobp_gsea)
}

KEGG_GSEA <- function(lfc_vector, gene_ident="gene_symbol", minGSSize = 25, maxGSSize = 500){
  set.seed(2020)
  kegg_gsea <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = minGSSize, # Minimum gene set size
    maxGSSize = maxGSSize, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      hs_kegg_sets,
      gs_name,
      gene_ident
    )
  )
  return(kegg_gsea)
}

REACTOME_GSEA <- function(lfc_vector, gene_ident="gene_symbol", minGSSize = 25, maxGSSize = 500){
  set.seed(2020)
  reactome_gsea <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = minGSSize, # Minimum gene set size
    maxGSSize = maxGSSize, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      hs_reactome_sets,
      gs_name,
      gene_ident
    )
  )
  return(reactome_gsea)
}

WIKIPATH_GSEA <- function(lfc_vector, gene_ident="gene_symbol", minGSSize = 25, maxGSSize = 500){
  set.seed(2020)
  wikipath_gsea <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = minGSSize, # Minimum gene set size
    maxGSSize = maxGSSize, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      hs_wikipath_sets,
      gs_name,
      gene_ident
    )
  )
  return(wikipath_gsea)
}

# Adapt the bubble plot code from the CellChat package
netVisual_bubble_custom <- function(object, 
                                    focus.cell = NULL, # <--- NEW ARGUMENT HERE
                                    sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, sort.by.source.priority = TRUE, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                                    comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                                    min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL, dot.size.min = NULL, dot.size.max = NULL,
                                    title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                                    grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                                    return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  
  # --- LOGIC FOR SINGLE OBJECT ---
  if (!is.list(object@net[[1]])) {
    message("Comparing communications on a single object \n")
    # (Standard parameter setup omitted for brevity, logic follows below)
    if (is.null(vjust.x) | is.null(hjust.x)) {
      angle=c(0, 45, 90); hjust=c(0, 1, 1); vjust=c(0, 1, 0.5)
      vjust.x = vjust[angle == angle.x]; hjust.x = hjust[angle == angle.x]
    }
    if (length(color.heatmap) == 1) {
      color.use <- tryCatch({ RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap) }, 
                            error = function(e) { scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors) })
    } else { color.use <- color.heatmap }
    if (direction == -1) { color.use <- rev(color.use) }
    
    if (is.null(comparison)) {
      cells.level <- levels(object@idents)
      if (is.null(sources.use)) sources.use <- cells.level
      if (is.null(targets.use)) targets.use <- cells.level
      
      df.net <- subsetCommunication(object, slot.name = "net", sources.use = sources.use, targets.use = targets.use, signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
      
      # --- CUSTOM FILTER INSERTION START ---
      if (!is.null(focus.cell)) {
        df.net <- df.net[df.net$source %in% focus.cell | df.net$target %in% focus.cell, ]
      }
      # --- CUSTOM FILTER INSERTION END ---
      
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      
      df.net$pval[df.net$pval > 0.05] = 1
      df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
      df.net$pval[df.net$pval <= 0.01] = 3
      df.net$prob[df.net$prob == 0] <- NA
      df.net$prob.original <- df.net$prob
      df.net$prob <- -1/log(df.net$prob)
      
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net <- with(df.net, df.net[order(interaction_name_2),])
      df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
      cells.order <- group.names
      df.net$source.target <- factor(df.net$source.target, levels = cells.order)
      df <- df.net
    }
  } else {
    # --- LOGIC FOR MERGED OBJECT (YOUR CASE) ---
    message("Comparing communications on a merged object \n")
    if (is.null(vjust.x) | is.null(hjust.x)) {
      angle=c(0, 45, 90); hjust=c(0, 1, 1); vjust=c(0, 1, 0.5)
      vjust.x = vjust[angle == angle.x]; hjust.x = hjust[angle == angle.x]
    }
    if (length(color.heatmap) == 1) {
      color.use <- tryCatch({ RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap) }, 
                            error = function(e) { scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors) })
    } else { color.use <- color.heatmap }
    if (direction == -1) { color.use <- rev(color.use) }
    
    cells.level <- levels(object@idents$joint)
    if (is.null(sources.use)) sources.use <- cells.level
    if (is.null(targets.use)) targets.use <- cells.level
    
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", sources.use = sources.use, targets.use = targets.use, signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.all <- data.frame()
    
    for (ii in 1:length(comparison)) {
      df.net <- df.net.all[[comparison[ii]]]
      
      # --- CUSTOM FILTER INSERTION START ---
      if (!is.null(focus.cell)) {
        # Keep rows where Source OR Target is the focus cell
        df.net <- df.net[df.net$source %in% focus.cell | df.net$target %in% focus.cell, ]
      }
      # --- CUSTOM FILTER INSERTION END ---
      
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      
      # Logic to handle isolates based ONLY on the filtered data
      # Note: We rebuild source.target list based on what is actually present if filtered
      if (!is.null(focus.cell)) {
        relevant.sources <- intersect(sources.use, unique(c(df.net$source, df.net$target))) # Simplified logic
        relevant.targets <- intersect(targets.use, unique(c(df.net$source, df.net$target)))
        # Actually, best to just let the code process what is there.
        # We skip the "forcing isolate" step for the filtered-out cells to prevent gaps.
        source.target.isolate <- character(0) 
      } else {
        source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
        source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      }
      
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], ")")
      
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      } else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), ncol = 5))
        colnames(df.net) <- c("interaction_name_2","source.target","prob","pval","prob.original")
        df.net$source.target <- group.names0
      }
      
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    
    # Process Probabilities
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T)*1.1, max(df.all$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2),])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
    
    # Clean Axis ordering
    cells.order <- c()
    dataset.name.order <- c()
    # Logic: Only iterate through group names that actually exist in the data now
    present.groups <- unique(gsub(" \\(.*\\)", "", df$source.target)) 
    
    for (i in 1:length(present.groups)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(present.groups[i], " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  
  # --- REMAINDER OF FUNCTION (PLOTTING) ---
  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  
  # Generate Plot
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  
  # Legends and Scales
  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  if (is.null(dot.size.max)) dot.size.max = max(df$pval)
  if (is.null(dot.size.min)) dot.size.min = min(df$pval)
  g <- g + scale_radius(range = c(dot.size.min, dot.size.max), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  
  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  # Grid lines
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  return(g)
}

# Plot gene expression based on pseudotime
# Import color scales
npg <- pal_npg("nrc")(10)
aaas <- pal_aaas()(10)
nature_scale <- c(npg, "darkgoldenrod1")
npg_scale_bars <- scale_fill_manual(values=nature_scale)
npg_scale <- scale_colour_manual(values = nature_scale)
aaas_scale <- c(aaas, "darkgoldenrod1")

plot_expression_on_pseudotime <- function(cds_ordered, cell_annotations, genes,
                                         facet_wrap_plot=FALSE,
                                         pseudotime_boundary = 30) { 
  
  # Subset cds containing selected cell type annotations and genes  
  cds_subset <- cds_ordered[rowData(cds_ordered)$gene_short_name %in% genes,
                           colData(cds_ordered)$level2_annotate %in% cell_annotations]
  
  # Create gene expression dataframe
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  colnames(cds_exprs) <- c("Feature_id", "Cell", "Norm_expression")
  
  # Add pseudotime values for each cell
  pseudotime <- pseudotime(cds_subset)
  cds_exprs$pseudotime <- pseudotime[match(cds_exprs$Cell, 
                                          names(pseudotime))]
  
  # Add annotations 
  annotations <- colData(cds_subset)$level2_annotate
  names(annotations) <- rownames(colData(cds_subset))
  cds_exprs$annotations <- annotations[match(cds_exprs$Cell, 
                                            names(annotations))]
  
  # Order panels according to order of genes provided
  cds_exprs$Feature_id = factor(cds_exprs$Feature_id,
                                levels = genes)
  
  # Fit a cubic spline to normalized expression values across pseudotime
  p <-  ggplot(cds_exprs, aes(pseudotime, Norm_expression, fill=Feature_id)) + 
    geom_smooth(method = "lm", formula = y ~ splines::ns(x, 3), 
                color="black") 
  p <- p + 
    xlim(0, pseudotime_boundary) + 
    scale_y_log10() + 
    ylab("Normalised expression") + 
    xlab("Pseudotime") +
    theme_bw() + 
    theme(aspect.ratio = 1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    npg_scale_bars
  
  # If we want to plot genes on different panels
  if (facet_wrap_plot) { 
    p <- p + 
      facet_wrap(~Feature_id, scales = "free_y") + 
      theme(legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "white", colour = "black"))
    return(p)
  }
  return(p)
}

miller_discrete_scale <- function(style="points", option=1) { 
  npg = pal_npg("nrc")(10)
  nature_scale = c(npg, "darkgoldenrod1")
  nature_scale2 = c("#6A3D9A", nature_scale[-1])
  
  if (style=="bars") {
    if (option==1) {
      discr_scale = scale_fill_manual(values = nature_scale)
    } else if (option==2) {
      discr_scale = scale_fill_manual(values = nature_scale2)
    }
  } else {
    if (option==1) {
      discr_scale = scale_colour_manual(values = nature_scale)
    } else if (option==2) {
      discr_scale = scale_colour_manual(values = nature_scale2)
    }
  }
  
  return(discr_scale)
}





