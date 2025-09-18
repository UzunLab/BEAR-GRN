library(Seurat)
library(dplyr)
library(patchwork)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(tidyverse)
#renv::install("Seurat@5.1.0")

data_dir <- 'RAW_DATA'

# Load scRNA-seq and scATAC-seq data using correct path
scRNA_data <- readRDS(paste0(data_dir, '/scRNA/scRNA.Seurat.rds'))
scATAC_data <- readRDS(paste0(data_dir, '/scATAC/scATAC.PeakMatrix_summarized_experiment.rds'))
scRNA_data
#Extract only E7.5_rep1 RNA data
E7.5_rna_cells <- grep("E7.5_rep1", colnames(scRNA_data), value = TRUE)
E7.5_rna_cells
scRNA_E7.5 <- subset(scRNA_data, cells = E7.5_rna_cells)
scRNA_E7.5
# E7.5_rep1 ATAC data
E7.5_atac_cells <- grep("E7.5_rep1", colnames(scATAC_data), value = TRUE)
E7.5_atac_cells
scATAC_E7.5 <- scATAC_data[, E7.5_atac_cells]
scATAC_E7.5
#find common cells
common_cells <- intersect(E7.5_rna_cells, E7.5_atac_cells)
length(common_cells)
# Subset RNA and ATAC data to keep only common cells
scRNA_common <- subset(scRNA_data, cells = common_cells)

scATAC_common <- scATAC_data[, common_cells]
scRNA_common
scATAC_common

# Create a Seurat object with both RNA and ATAC data
multiomic_data <- CreateSeuratObject(counts = scRNA_common@assays$RNA@counts)
multiomic_data
multiomic_data[["ATAC"]] <- CreateAssayObject(counts = assay(scATAC_common, "PeakMatrix"))

# Add metadata if needed
atac_metadata <- as.data.frame(colData(scATAC_common))
atac_metadata
rownames(atac_metadata) <- colnames(scATAC_common)
atac_metadata
multiomic_data <- AddMetaData(multiomic_data, atac_metadata)
multiomic_data@assays$RNA
multiomic_data@assays$ATAC
VlnPlot(multiomic_data,
        features = c("nFeature_RNA", # Number of genes per cell
                     "nCount_RNA", # Number of RNA per cell
                     "mitochondrial_percent_RNA", # Percent mitochondrial
                     "nFeature_ATAC",
                     "TSSEnrichment_atac",
                     "FRIP"# Number of peaks
        ),
        ncol = 3,
        pt.size = 0.1)

multiomic_data_filter <- subset(multiomic_data,
                              subset = nFeature_RNA > 1000 &
                                nFeature_RNA < 7500 &
                                mitochondrial_percent_RNA < 20 &
                                nFeature_ATAC > 1000 &
                                nFeature_ATAC < 30000 &
                                TSSEnrichment_atac > 1 &
                                FRIP > 0.3)


multiomic_data_filter
saveRDS(multiomic_data_filter,"multiomic_data_filtered_L1_E7.5_rep1.rds")
write.csv(multiomic_data_filter@assays$RNA$counts,"multiomic_data_filtered_L1_E7.5_rep1_RNA.csv")
write.csv(multiomic_data_filter@assays$ATAC$counts,"multiomic_data_filtered_L1_E7.5_rep1_ATAC.csv")
multiomic_data_filter@assays$ATAC$counts
multiomic_data@meta.data
