library(Seurat)

# Base path where all folders are
base_path <- "/home/pakorns/hp-storage/scRNAseq/GSE165087"

# List subdirectories and select only contain "GSM"
subdirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
subdirs <- subdirs[grepl("/GSM", subdirs)]

# Function to load, QC, and return filtered Seurat object
load_and_qc <- function(dir) {
  rds_file <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_file) == 1) {
    seurat_obj <- readRDS(rds_file)

    # Add percent.mt (adjust pattern depending on species)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    # Optional: print summary before filtering
    print(paste0("QC summary for: ", basename(dir)))
    print(summary(seurat_obj[["nFeature_RNA"]]))
    print(summary(seurat_obj[["nCount_RNA"]]))
    print(summary(seurat_obj[["percent.mt"]]))

    # Filter cells (adjust thresholds to suit your data)
    seurat_obj <- subset(seurat_obj, 
                         subset = nFeature_RNA > 200 & 
                                  nFeature_RNA < 2500 &
				  nCount_RNA < 10000 & 
                                  percent.mt < 20)

    # Track source
    seurat_obj$orig.ident <- basename(dir)

    return(seurat_obj)
  } else {
    warning(paste("No or multiple RDS files in", dir))
    return(NULL)
  }
}

# Apply QC to all datasets
seurat_list <- lapply(subdirs, load_and_qc)
seurat_list <- Filter(Negate(is.null), seurat_list)  # Remove NULLs

# Merge cleaned datasets
combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = basename(subdirs))

# Save merged object
saveRDS(combined_seurat, file = file.path(base_path, "GSE165087_QC_merged.rds"))

cat("✅ Merged object saved as GSE165087_QC_merged.rds\n")
