# Separate the data for Nadeu et al.

# Import the relevant packages
library(dplyr)
library(Seurat)

# Set working directory
setwd("D:/Bioinfor PhD/Thesis/Progress lab meeting/Datasets/Measure_consistency")

# Load the dataset
Nadeu.seu <- readRDS('3.seurat_filtered.rds')
str(Nadeu.seu)

# Classidy individual datasets to CLL or RT
Nadeu.seu@meta.data$cond[Nadeu.seu@meta.data$sample_description_FN == 'RS'] <- 'RT'
Nadeu.seu@meta.data$cond[Nadeu.seu@meta.data$sample_description_FN != 'RS'] <- 'CLL'

# Classidy patients to their id
unique(Nadeu.seu[['library_name']])

Nadeu.seu@meta.data <- Nadeu.seu@meta.data %>%
  mutate(case = ifelse(grepl('012', library_name), 12, 
                       ifelse(grepl('019', library_name), 19, 
                              ifelse(grepl('3299', library_name), 3299,
                                     ifelse(grepl('365', library_name), 365,
                                            ifelse(grepl('15-194', library_name), 365,
                                                   ifelse(grepl('05_662', library_name), 12,
                                                          ifelse(grepl('00_30', library_name), 12, 19))))))))
unique(Nadeu.seu[['case']])

Nadeu.seu[['ref']] <- 'Nadeu et al.'

# Select patients with PB tissue and CLL post-treatment3
case12.CLL_PB_PT3 <- subset(x=Nadeu.seu, subset = case=='12' & tissue=='PB' & sample_description_FN=='posttreatment_3')
case19.CLL_PB_PT3 <- subset(x=Nadeu.seu, subset = case=='19' & tissue=='PB' & sample_description_FN=='posttreatment_3')

# Select patients with PB tissue and RT
case12.RT_PB <- subset(x=Nadeu.seu, subset = case=='12' & tissue=='PB' & sample_description_FN=='RS')
case19.RT_PB <- subset(x=Nadeu.seu, subset = case=='19' & tissue=='PB' & sample_description_FN=='RS')

# Save .RDS files
SaveSeuratRds(case12.CLL_PB_PT3, file = "case12_CLL_PB_PT3.Rds")
SaveSeuratRds(case19.CLL_PB_PT3, file = "case19_CLL_PB_PT3.Rds")
SaveSeuratRds(case12.RT_PB, file = "case12_RT_PB.Rds")
SaveSeuratRds(case19.RT_PB, file = "case19_RT_PB.Rds")
