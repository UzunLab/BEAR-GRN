#renv_lockfile_load(project = New, strict = TRUE)
setwd("~/Multi_omics_GRN/Pando/New/K562")
#renv_lockfile_load(project = New, strict = TRUE)
#renv::status()
renv::activate()
renv::restore()

# Load libraries
library(Seurat)
#install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#install.packages("Signac")
library(Signac)
library(Matrix)
#devtools::install_github('quadbio/Pando')
#renv::snapshot()
#setwd("~/Multi_omics_GRN/Pando/New/K562/")
# Load scRNA and scATAC data from text files
scRNA_data <- read.table("subsampled_RNA_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
scATAC_data <- read.table("subsampled_ATAC_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
# Convert scRNA data to sparse matrix
#scRNA_data_sparse <- as(as.matrix(scRNA_data), "dgCMatrix")
# Create Seurat object with scRNA-seq data
seurat_rna <- CreateSeuratObject(counts = scRNA_data, assay = "RNA")
# Convert scATAC data to sparse matrix
scATAC_data_sparse <- as(as.matrix(scATAC_data), "sparseMatrix")
# Create ChromatinAssay object with scATAC data
#chrom_assay <- CreateChromatinAssay(counts = scATAC_data, sep = c(":", "-"), assay = "peaks")

# Add the chromatin assay to the Seurat object
#seurat_rna[["peaks"]] <- chrom_assay
#Assays(seurat_rna)
#head(seurat_rna[["peaks"]]@counts)

# Load libraries
library(GenomicRanges)
#library(Signac)
#library(rtracklayer)
###########################################################
# Load the GTF file for hg38
#gtf_file <- "gencode.v38.annotation.gtf"  # Replace with your GTF file path
#gene_annotations <- rtracklayer::import(gtf_file)

# Filter the GTF file to include only gene annotations
#gene_annotations <- gene_annotations[gene_annotations$type == "gene"]
# Add the ChromatinAssay to the Seurat object (if not already done)
# Let's assume 'chrom_assay' is your ATAC data
#seurat_rna[["peaks"]] <- chrom_assay  # Replace 'chrom_assay' with your actual object

# Add gene annotations to the ChromatinAssay
#Annotation(seurat_rna[["peaks"]]) <- gene_annotations
#################################
# add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(scATAC_data), sep = c(":", "-"))
grange.counts
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac.counts <- scATAC_data[as.vector(grange.use), ]
dim(atac.counts)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
#saveRDS(annotations,"EnsDb.Mmusculus.v79_UCSC_mm10.rds")
#file <- "e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"
#path <- file.path(dir.in, file)
#annotations <- readRDS("../EnsDb.Mmusculus.v79_UCSC_mm10.rds")
chrom.assay <- CreateChromatinAssay(
  counts = atac.counts,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  annotation = annotations
)
chrom.assay
seurat_rna[["peaks"]] <- chrom.assay
#e18[["ATAC"]] <- chrom.assay
# Preprocess RNA data
seurat_rna <- NormalizeData(seurat_rna)
#seurat_rna <- FindVariableFeatures(seurat_rna)

# Preprocess ATAC data
seurat_rna <- RunTFIDF(seurat_rna, assay = "peaks")
#seurat_rna <- FindTopFeatures(seurat_rna, assay = "peaks", min.cutoff = "q0")
#seurat_rna <- RunSVD(seurat_rna, assay = "peaks")

# Optionally, integrate both assays
#seurat_rna <- FindMultiModalNeighbors(seurat_rna, reduction.list = list("pca", "lsi"), 
#                                     dims.list = list(1:30, 1:30), modality.weight.name = "RNA_ATAC.weight")
#seurat_rna <- RunUMAP(seurat_rna, nn.name = "weighted.nn", reduction.name = "wnn.umap")


library(Pando)
seurat_object <- initiate_grn(seurat_rna)

library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs)
motifs
seurat_object <- find_motifs(
  seurat_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
genes <- row.names(scRNA_data)
seurat_object <- infer_grn(
  seurat_object,
  genes = genes,
  peak_to_gene_method = 'Signac',
  method = "glm"
  )
	     
#"glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge",             "bayesian_ridge"))
coef(seurat_object)
#seurat_object <- find_modules(seurat_object)
#modules <- NetworkModules(seurat_object)
#modules@meta
GRN <- GetGRN(seurat_object)
#saveRDS(GRN, file = "regulatory_network.rds")
#save network as csv
grn <- GRN@networks$glm_network@coefs
write.csv(grn,"K562_raw_network.csv",row.names = FALSE)
library(dplyr)
#network <- readRDS("regulatory_network.rds")

# Filter based on p-value threshold (e.g., pval < 0.05)
filtered_data <- GRN@networks$glm_network@coefs %>%
  filter(pval < 0.05)

# Save the filtered data into a CSV file
write.csv(filtered_data, "K562_filtered_network.csv", row.names = FALSE)
