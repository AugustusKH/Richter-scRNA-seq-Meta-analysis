## Attach the packages you will need for the analysis.
library(GEOquery)
library(glue)

## Check root path
here::here()

## Define paths to load and save 
path_to_study <- here::here("scRNAseq/Hing_2023_GSE183426")

## Download GEO files
# Collect GSM_IDs in a list
gse <- getGEO("GSE183426", GSEMatrix = FALSE)
gsm_list <- GSMList(gse)

for (gsm in names(gsm_list)) {
  print(gsm)
  download_dir <- str_c(path_to_study, "/", gsm)
  dir.create(download_dir, recursive = TRUE)
  gsm_data <- gsm_list[[gsm]]
  for (i in 1:4) {
    url <- gsm_data@header[[glue("supplementary_file_{i}")]]
    url <- sub("^ftp", "http", url)
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











