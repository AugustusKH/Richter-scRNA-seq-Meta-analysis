# Data merging and integration

# Import the relevant package
library(Seurat)
library(harmony)
library(dplyr)
library(SingleR)
library(celldex)
library(sctransform)

# Source functions
source(here::here("bin/utils.R"))

# Define paths to load and save
path_to_penter_lst <- here::here("scRNAseq/systematic_review/1-Penter_2021_GSE165087/R_objects/Penter_seu_lst.rds")
path_to_nadeu_lst <- here::here("scRNAseq/systematic_review/2-Nadeu_2022/R_objects/Nadeu_qc_seu_lst.rds")
path_to_rejeski_lst <- here::here("scRNAseq/GSE201704/GSM6069857/R_objects/Rejeski_2022_GSE201704_seu_lst.rds")
path_to_hing_lst <- here::here("scRNAseq/systematic_review/3-Hing_2023_GSE205704/R_objects/Hing_seu_lst.rds")
path_to_parry_a_lst <- here::here("scRNAseq/Parry_2023a/R_objects/Parry_2023a_seu_lst.rds")
path_to_parry_b_lst <- here::here("scRNAseq/systematic_review/4-Parry_2023b/R_objects/Parry_seu_lst.rds")
path_to_rigo_lst <- here::here("scRNAseq/systematic_review/5-Rigo_2024/R_objects/Rigo_seu_lst.rds")
path_to_object <- here::here("scRNAseq/meta_analysis/R_objects")
path_to_save_all_lst <- paste0(path_to_object, "/all_seu_lst.rds")
path_to_save_merged_obj <- paste0(path_to_object, "/merged_seu.rds")
path_to_save_downstream_merged_obj <- paste0(path_to_object, "/downstream_merged_seu.rds")
path_to_save_integrated_lst <- paste0(path_to_object, "/integrated_seu_lst.rds")

# Load each data sets and collect into one list
path_lst <- c(path_to_penter_lst, path_to_nadeu_lst, path_to_rejeski_lst, path_to_hing_lst, path_to_parry_a_lst, path_to_parry_b_lst, path_to_rigo_lst)
seurat_lst <- purrr::map(path_lst, function(path){
  lst <- readRDS(path)
  seurat.obj <- purrr::map(lst, function(obj){
    return(obj) 
  })
  return(seurat.obj)
})

# Flatten (combine) all Seurat objects from all studies
seurat_all <- purrr::flatten(seurat_lst)

# Save the list
saveRDS(seurat_all, path_to_save_all_lst)

# Merge all data sets
seurat_all <- readRDS(path_to_save_all_lst)
seurat <- seurat_all[[1]]
for (i in seq(2, length(seurat_all))) {
  message("Merging object ", i, " of ", length(seurat_all))
  seurat <- merge(x = seurat, y = seurat_all[[i]])
  gc()
  if (i %% 2 == 0) {
    saveRDS(seurat, file = paste0(path_to_object, "/merged_partial_", i, ".rds"))
  }
}
seurat$patient <- paste0(seurat$case, "_", seurat$ref) # Some cases have the same id but different across those studies
seurat$sample_id <- seurat$case_id # Convert to the right name as case_id refer truly to samples not cases

# Free up memory afterwards
rm(seurat_all)
gc()

# Save the merged object
saveRDS(seurat, path_to_save_merged_obj)

# Downstream analysis for the merged object
DefaultAssay(seurat) <- "RNA"
seurat <- JoinLayers(seurat)
seurat <- SCTransform(seurat, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
seurat <- RunPCA(seurat)
saveRDS(seurat, path_to_save_downstream_merged_obj) # save the merge object after downstream analysis

# Integration using harmony for various metadata variables
vars_list <- list(
  sample = "sample_id",
  ref  = "ref",
  case = "patient",
  sample_ref = c("sample_id", "ref"),
  sample_case = c("sample_id", "patient"),
  ref_case = c("ref", "patient"),
  all = c("sample_id", "ref", "patient")
)

harmony_results <- list()

for (name in names(vars_list)) {
  message("Running Harmony for: ", name)
  
  so <- seurat
  
  so <- RunHarmony(
    so,
    group.by.vars = vars_list[[name]],
    dims = 1:30,
    verbose = FALSE
  )
  
  so <- RunUMAP(so, dims = 1:30, reduction = "harmony")
  so <- FindNeighbors(so, dims = 1:30, reduction = "harmony")
  so <- FindClusters(
    so,
    random.seed = 12345,
    verbose = FALSE,
    resolution = seq(0.1, 0.9, by = 0.1)
  )
  
  so <- singleR_cell_annotation(so, annotation="both")
  so$harmony_model <- name
  
  harmony_results[[name]] <- so
  
  saveRDS(so, paste0(path_to_object, "/harmony_", name, ".rds"))
}

# Save the integrated list
saveRDS(harmony_results, path_to_save_integrated_lst)
