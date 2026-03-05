# Integrative single-cell analysis in Richter's Transformation
Richter’s transformation (RT), characterised by progression from chronic lymphocytic leukaemia (CLL) to an aggressive lymphoma phenotype, remains associated with poor clinical outcomes. This progression is strongly influenced by the tumour microenvironment (TME), which promotes tumour survival and facilitates immune evasion. Single-cell RNA sequencing (scRNA-seq) analyses have been conducted to investigate tumour biology in RT; however, studies focusing on the TME remain limited, largely due to the insufficient representation of microenvironmental cells, which constrains statistical power. 

Herein, we performed an integrative single-cell data analysis with batch correction across seven human RT studies identified through PubMed, Scopus, and Web of Science. The integrated dataset, comprising 593,424 cells, enabled the construction of a simulated biological system of RT progression and facilitated the investigation of tumour-microenvironment interplay through gene expression profiling, cell–cell communication, and trajectory analyses. Consistent with findings from previous individual tumour studies, tumour heterogeneity remained evident following batch correction. Nevertheless, our integrated analysis demonstrated that RT tumour cells exhibited reduced reliance on extrinsic activation pathways, suggesting an adaptive mechanism associated with immune evasion. This shift was accompanied by enhanced intrinsic activation and metabolic reprogramming towards aerobic glycolysis (the Warburg effect), collectively promoting tumour cell survival and proliferation. We also identified novel candidate RT-related genes, including GINS2 and SLC25A4, which were uniquely upregulated differentially expressed genes (DEGs) detected in the integrated dataset but not showed in individual human bulk or single-cell studies. Notably, both genes demonstrated consistent upregulation across all five murine bulk datasets. 

Microenvironmental and cell-cell communication analyses indicated that monocytes represent a key driver of RT progression, characterised by enhanced BAFF signalling between tumour and monocytic cells. This interaction was supported by the upregulation of TNFSF13B within the nurse-like cell (NLC)-like monocyte subtype. An imbalance in tissue sampling between RT and CLL conditions was also observed, which may impose limitations on TME analyses. Nevertheless, taken together, our integrated analysis extends the biological understanding of RT by enabling comprehensive characterisation of both tumour and microenvironmental cell populations, highlighting potential implications for future therapeutic strategies.

All code and scripts used throughout this study are included in this repository and are publicly available to ensure reproducibility.

## Data accessibility

All analysed datasets can be downloaded from [https://doi.org/10.5281/zenodo.18874149]. There are 6 main available files. The first file is raw count matrices including metadata for individual bulk and single-cell transcriptomic studies. The file can be downloaded directly from the Zenodo link or via the command line below:

```bash
  wget https://zenodo.org/records/18874149/files/raw_count_matrices.zip
  unzip raw_count_matrices.zip
  cd raw_count_matrices
```
This command can also be used to download other ZIP files, such as `mouse_objects.zip` and `microenvironment.zip`. These files contain analysed Seurat objects derived from mouse scRNA-seq data and integrated human non-tumour cell datasets, respectively. 

We also provide additional objects obtained from human single-cell data integration as `.rds` files (e.g. `all_seu_lst.rds`, `merged_seu.rds`, and `integrated_seu.rds`). These files can be downloaded either from the link above or directly via the command line, as shown below, without requiring extraction.

```bash
  wget https://zenodo.org/records/18874149/files/integrated_seu.rds
```
## Important packages

We list the essential packages used in this study, along with the versions used, below. Additional packages used for the individual method are shown on `Session Information` at the end of the Quarto HTML reports. 

* [Seurat 5.3.1 ](https://satijalab.org/seurat/articles/get_started_v5_new)
* sctransform_0.4.2
* [harmony 1.2.4](https://github.com/immunogenomics/harmony)
* [SingleR 2.8.0](https://github.com/dviraran/SingleR)
* [DESeq2 1.46.0](https://github.com/thelovelab/DESeq2)
* [monocle3 1.4.26](https://cole-trapnell-lab.github.io/monocle3/)
* [CellChat 2.2.0](https://github.com/jinworks/CellChat)
* clusterProfiler 4.14.6
* msigdbr 25.1.1
* org.Hs.eg.db 3.20.0
* EnhancedVolcano 1.24.0
* UpSetR 1.4.0
* stringr 1.6.0
* dplyr 1.1.4
* tibble 3.3.0
* ggplot2 4.0.1
* pheatmap 1.0.13
* ComplexHeatmap 2.22.0
* DT 0.34.0

## Publication article

Further details in both methodology and results are described in our article via this [link]()












