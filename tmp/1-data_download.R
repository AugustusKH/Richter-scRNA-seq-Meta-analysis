# This script is used for downloading GSE183426 which is mouse scRNA-seq count matrices for Hing et al. (2023)

## Import the relevant packages
library(GEOquery)
library(glue)

## Load the data
gse <- getGEO("GSE183426", GSEMatrix = FALSE)
saveRDS(gse, "gse_obj.rds")
gsm_list <- GSMList(gse)

for (gsm_name in names(gsm_list)) {
  download_dir <- here::here(gsm_name)
  dir.create(download_dir, recursive = TRUE)
  gsm <- gsm_list[[gsm_name]]
  for (i in 1:3) {
    url <- str_c("https", substring(gsm@header[[glue("supplementary_file_{i}")]], 4))
    command <- glue::glue("wget -P {download_dir} {url}")
    system(
      command,
      intern = FALSE,
      ignore.stdout = TRUE,
      ignore.stderr = TRUE,
      wait = TRUE
    )
  }
}




