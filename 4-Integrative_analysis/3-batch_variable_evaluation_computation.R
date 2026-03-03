# Batch removal and biological preservation evaluation

# Import the relevant package
library(kBET)
library(FNN)
library(lisi)
library(cluster)
library(Seurat)
library(stats)
library(DescTools)
library(data.table)
library(tidyverse)
library(ggsci)

# Source functions
source(here::here("bin/utils.R"))

# Define paths to load and save
path_to_obj <- here::here("scRNAseq/meta_analysis/R_objects")
path_to_integrated_lst <- paste0(path_to_obj, "/integrated_seu_lst.rds")
path_to_save_lisi_df <- paste0(path_to_obj, "/lisi_df.rds")
path_to_save_kbet_df <- paste0(path_to_obj, "/kbet_df.rds")
path_to_save_sil_df <- paste0(path_to_obj, "/sil_df.rds")
path_to_save_pcr_df <- paste0(path_to_obj, "/pcr_df.rds")

# Load the integrated list
harmony_seu.lst <- readRDS(path_to_integrated_lst) 

# Named list of datasets + their batch variables
harmony_inputs <- list(
  sample       = list(obj = harmony_seu.lst$sample,       batch = "sample_id"),
  ref          = list(obj = harmony_seu.lst$ref,          batch = "ref"),
  case         = list(obj = harmony_seu.lst$case,         batch = "patient"),
  sample_ref   = list(obj = harmony_seu.lst$sample_ref,   batch = "ref"),
  sample_case  = list(obj = harmony_seu.lst$sample_case,  batch = "patient"),
  ref_case     = list(obj = harmony_seu.lst$ref_case,     batch = "ref"),
  all          = list(obj = harmony_seu.lst$all,          batch = "ref")
)

# Define a function for subsampling
subsample_seu <- function(seu, prop = 0.5) {
  set.seed(NULL)
  n <- ncol(seu)
  keep <- sample(colnames(seu), size = floor(prop * n))
  return(seu[, keep])
}

# Set sampling parameters
n_repeats <- 5
subsample_prop <- 0.04

##########################
# 1. LISI scores         #
##########################

# Define function to compute batch LISI (bLISI), cell LISIS (cLISI), and cluster LISI (clLISI)
compute_lisi_scores <- function(seu, batch_var, label_var = c("singleR.labels_main", "SCT_snn_res.0.3")) {
  pca_list <- seu@reductions$pca
  emb <- pca_list@cell.embeddings[, 1:30]
  meta <- seu@meta.data[, c(batch_var, label_var), drop = FALSE]
  lisi_res <- compute_lisi(emb, meta, c(batch_var, label_var))
  df <- as.data.frame(lisi_res) # convert to data frame with clear names
  colnames(df) <- c("bLISI", "cLISI", "clLISI")
  return(df)
}

# Compute all LISI results
lisi_results <- list()

for (rep in 1:n_repeats) {
  message("LISI subsampling repeat: ", rep)
  
  tmp <- lapply(names(harmony_inputs), function(name) {
    item <- harmony_inputs[[name]]
    seu_sub <- subsample_seu(item$obj, prop = subsample_prop) # subsample
    df <- compute_lisi_scores(seu_sub, batch_var = item$batch) # compute
    df$dataset <- name
    df$repeat_n  <- rep
    return(df)
  })
  lisi_results[[rep]] <- rbindlist(tmp)
}

# Merge all into one tidy data frame
lisi_df <- rbindlist(lisi_results, use.names = TRUE, fill = TRUE)
saveRDS(lisi_df, path_to_save_lisi_df)

##########################
# 2. kBET                #
##########################

# Define function to compute kBET
compute_kbet <- function(seu, batch_var, k_range) {
  pca_list <- seu@reductions$pca
  emb <- pca_list@cell.embeddings[, 1:30]
  batch <- seu@meta.data[[batch_var]]
  kbet_list <- vector("list", length(k_range))
  
  for (i in seq_along(k_range)) {
    
    batch_est <- kBET(
      df = emb, 
      batch = batch, 
      k0 = k_range[i], 
      do.pca = FALSE, 
      plot = FALSE, 
      heuristic = FALSE,
      adapt = FALSE,
      verbose = TRUE
    )
    
    df <- data.frame(
      class = rep(c("observed", "expected"),
                  each = length(batch_est$stats$kBET.observed)),
      data  = c(batch_est$stats$kBET.observed,
                batch_est$stats$kBET.expected),
      k     = paste0("k_", k_range[i])
    )
    
    kbet_list[[i]] <- df
  }
  return(rbindlist(kbet_list, use.names = TRUE, fill = TRUE))
}

# Define a base for k range
base_k <- c(10, 25, 50, 100, 500, 1000, 2500)

# Compute all kBET results
kbet_results <- list()

for (rep in 1:n_repeats) {
  message("kBET subsampling repeat: ", rep)
  
  tmp <- lapply(names(harmony_inputs), function(name) {
    item <- harmony_inputs[[name]]
    seu_sub <- subsample_seu(item$obj, prop = subsample_prop)
    
    # Set k_range
    batch_vec <- seu_sub@meta.data[[item$batch]]
    batch_labels <- as.factor(batch_vec)
    k0 = floor(mean(table(batch_labels))) / 4 # set k0 as the original article suggested
    k_range <- unique(sort(c(base_k, k0)))
    n_cells <- ncol(seu_sub)
    k_range <- k_range[k_range < n_cells] # remove impossible k
    
    # Run kBET
    df <- compute_kbet(seu_sub, batch_var = item$batch, k_range)
    df$dataset <- name
    df$repeat_n  <- rep
    return(df)
  })
  kbet_results[[rep]] <- rbindlist(tmp)
}

# Merge all into one tidy data frame
kbet_df <- rbindlist(kbet_results, use.names = TRUE, fill = TRUE)
saveRDS(kbet_df, path_to_save_kbet_df)

##########################
# 3. Silhouette scores   #
##########################

# Define a function to compute cluster average sihouette width (clASW)
compute_sil <- function(seu, res){
  pca_list <- seu@reductions$pca
  emb <- pca_list@cell.embeddings[, 1:30]
  dist_mat <- dist(emb)
  sil_list <- vector("list", length(res))
  
  for (i in seq_along(res)){
    clusters <- seu[[paste0("SCT_snn_res.", res[i])]][,1]
    sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist_mat)
    df <- data.frame(cluster=sil[, 1], neighbor=sil[, 2], sil_width=sil[, 3], res=res[i])
    sil_list[[i]] <- df
  }
  return(rbindlist(sil_list, use.names = TRUE, fill = TRUE))
}

# Define a cluster resolution range
res <- seq(0.1, 0.9, by = 0.1)

# Compute all silhouette results
sil_results <- list()

for (rep in 1:n_repeats) {
  message("Silhouette subsampling repeat: ", rep)
  
  tmp <- lapply(names(harmony_inputs), function(name) {
    item <- harmony_inputs[[name]]
    seu_sub <- subsample_seu(item$obj, prop = subsample_prop)
    df <- compute_sil(seu_sub, res)
    df$dataset <- name
    df$repeat_n  <- rep
    return(df)
  })
  sil_results[[rep]] <- rbindlist(tmp)
}

# Merge all into one tidy data frame
sil_df <- rbindlist(sil_results, use.names = TRUE, fill = TRUE)
saveRDS(sil_df, path_to_save_sil_df)

##########################
# 4. PCRegression        #
##########################

# Define a function to compute PCR
compute_pcr <- function(seu, batch_var){
  pca_list <- seu@reductions$pca
  pca_data <- list(x = pca_list@cell.embeddings, sdev = pca_list@stdev)
  class(pca_data) <- "prcomp"  # pcRegression expects a prcomp object
  batch <- as.factor(seu[[batch_var, drop = TRUE]])  # use column dynamically
  pcr <- kBET::pcRegression(pca.data=pca_data,
                            batch=batch,
                            n_top=30)
  return(pcr$R2Var)
}

# Compute all PCR results
pcr_results <- list()

for (rep in 1:n_repeats) {
  message("PCR subsampling repeat: ", rep)
  
  tmp <- lapply(names(harmony_inputs), function(name) {
    item <- harmony_inputs[[name]]
    seu_sub <- subsample_seu(item$obj, prop = subsample_prop)
    R2Var <- compute_pcr(seu_sub, batch_var = item$batch)
    df <- data.frame(R2Var = R2Var, dataset = name, repeat_n = rep)
    return(df)
  })
  pcr_results[[rep]] <- rbindlist(tmp)
}

# Merge all into one tidy data frame
pcr_df <- rbindlist(pcr_results, use.names = TRUE, fill = TRUE)
saveRDS(pcr_df, path_to_save_pcr_df)




