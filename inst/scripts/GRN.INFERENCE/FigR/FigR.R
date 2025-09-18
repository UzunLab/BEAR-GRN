library(FigR)
setwd("/gpfs/Home/kmk7420/Multi_omics_GRN/FigR/K562")
rna_counts <- read.table("subsampled_RNA_K562_human_filtered.csv", sep = ",", row.names = 1, header = TRUE)
atac_counts <- read.table("subsampled_ATAC_K562_human_filtered.csv", sep = ",", row.names = 1, header = TRUE)
# Filter out rows where row names do not start with 'chr' or are 'chrM'
atac_counts <- atac_counts[grepl("^chr[0-9XY]", rownames(atac_counts)), ]
rna_counts
atac_counts
#BiocManager::install("SummarizedExperiment")

# Load necessary library
library(SummarizedExperiment)
library(GenomicRanges)
# Assuming you have peak counts in a matrix (peaks x cells)
peak_counts_matrix <- as.matrix(atac_counts)  # Your ATAC count matrix

# Check the first few entries
#head(peak_granges)
#peak_granges# Extract chromosome, start, and end positions from row names
peak_info <- do.call(rbind, strsplit(rownames(peak_counts_matrix), "[:|-]"))
#peak_info
chromosomes <- peak_info[, 1]
chromosomes
start_positions <- as.numeric(peak_info[, 2])
end_positions <- as.numeric(peak_info[, 3])

# Create the GRanges object
peak_granges <- GRanges(
  seqnames = chromosomes, 
  ranges = IRanges(start = start_positions, end = end_positions)
)

# Check the first few entries
head(peak_granges)
# Create a RangedSummarizedExperiment
atac_se <- SummarizedExperiment(assays = list(counts = peak_counts_matrix),
                                      rowRanges = peak_granges)
assay(atac_se, "counts") <- as(assay(atac_se, "counts"), "dgCMatrix")
# Check the class of the created object
class(atac_se)# Assuming atac_counts is the scATAC-seq matrix (peaks x cells)
# Create a SummarizedExperiment object
#atac_se <- SummarizedExperiment(
#  assays = list(counts = atac_counts),  # Add the count matrix as an assay
#  rowData = DataFrame(peak = rownames(atac_counts)),  # Assign peak names to rowData
#  colData = DataFrame(cell = colnames(atac_counts))  # Assign cell names to colData
#)

# Check the SummarizedExperiment object
#atac_se

# Load necessary library
library(Matrix)

# Assuming rna_counts is your scRNA-seq matrix (genes x cells)
# Convert to a sparseMatrix
rna_counts_matrix <- as.matrix(rna_counts)

rna_sparse <- as(rna_counts_matrix, "CsparseMatrix")
# Check the sparseMatrix object
rna_sparse
# Run using multiple cores if parallel support
cisCor <- runGenePeakcorr(ATAC.se = atac_se,
                          RNAmat = rna_sparse,
                          genome = "hg38", # Also supports mm10 and hg38
                          nCores = 12, 
                          p.cut=NULL)
saveRDS(cisCor,"cisCor.rds")
#cisCor <- readRDS("cisCor.rds")
# Filter peak-gene correlations by p-value                    
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)

# Determine DORC genes
dorcGenes <- cisCor.filt %>% dorcJPlot(cutoff=1, # Default
                                       returnGeneList = TRUE)
library(FNN)
# Get DORC scores
dorcMat <- getDORCScores(ATAC.se = atac_se,dorcTab=cisCor.filt,geneList=dorcGenes,nCores=12)
#dorcMat
######################test#################################
#cellKNN.mat <- knn.index(atac_counts, k=4)
library(Seurat)
library(Signac)
common_cells <- intersect(colnames(rna_sparse), colnames(atac_se))
#common_cells
rna_subset <- rna_sparse[, common_cells]
atac_subset <- atac_se[, common_cells]
#rna_subset
#atac_subset

# RNA Seurat Object
seurat_rna <- CreateSeuratObject(counts = rna_subset, assay = "RNA")
num_genes <- dim(seurat_rna)[1]
# ATAC Seurat Object
seurat_atac <- CreateSeuratObject(counts = assays(atac_subset)$counts, assay = "ATAC")

seurat_rna[["ATAC"]] <- CreateAssayObject(counts = assays(atac_subset)$counts)

seu.multi <- seurat_rna
seu.multi
# Normalize RNA data
seu.multi <- NormalizeData(seu.multi)
seu.multi <- FindVariableFeatures(seu.multi, selection.method = "vst", nfeatures = num_genes)
seu.multi <- ScaleData(seu.multi)
seu.multi <- RunPCA(seu.multi, npcs = 50)

# Normalize and process ATAC data
seu.multi <- RunTFIDF(seu.multi, assay = "ATAC")
seu.multi <- FindTopFeatures(seu.multi, min.cutoff = 1, assay = "ATAC")
seu.multi <- RunSVD(seu.multi, assay = "ATAC", n = 50)

# Find multi-modal neighbors (WNN)
seu.multi <- FindMultiModalNeighbors(
  seu.multi, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight",
  k.nn = 30
)

# Extract WNN kNN matrix
cellKNN.mat <- seu.multi@neighbors$weighted.nn@nn.idx
head(cellKNN.mat)

rownames(cellKNN.mat) <- colnames(dorcMat)
############################################################
# Smooth DORC scores (using cell KNNs)
library(doParallel)
library(foreach)
dorcMat.smooth <- smoothScoresNN(NNmat=cellKNN.mat,mat=dorcMat,nCores=12)
dorcMat.smooth

rnaMat.smooth <- smoothScoresNN(NNmat = cellKNN.mat, mat = rna_sparse, nCores = 12)
#dorcMat.smooth[is.na(dorcMat.smooth)] <- 0
#dorcMat.smooth <- dorcMat.smooth[rowSums(is.na(dorcMat.smooth)) == 0, ]
#dorcMat.smooth
dorcMat.smooth <- as.matrix(dorcMat.smooth)
# Compute mean and standard deviation
col_means <- colMeans(dorcMat.smooth, na.rm = TRUE)
col_sds <- apply(dorcMat.smooth, 2, sd, na.rm = TRUE)

# Replace 0 standard deviations with 1 to avoid division by zero
col_sds[col_sds == 0] <- 1  

# Standardize the matrix
dorcMat.smooth <- sweep(dorcMat.smooth, 2, col_means, FUN = "-") 
dorcMat.smooth <- sweep(dorcMat.smooth, 2, col_sds, FUN = "/")

# Confirm that no NAs exist
sum(is.na(dorcMat.smooth))  # Should return 0

# Run FigR
library(parallel)
nCores <- 12  # Adjust this based on your system
cl <- makeCluster(nCores, type = "PSOCK")  # Use PSOCK for cross-platform compatibility
registerDoParallel(cl)
foreach::getDoParWorkers()  # Should return 12

fig.d <- runFigRGRN(ATAC.se= atac_se,
                    rnaMat=rnaMat.smooth, # Smoothed RNA matrix using paired cell kNNs
                    dorcMat=dorcMat.smooth,
                    dorcTab=cisCor.filt,
                    genome="hg38",
                    dorcGenes=dorcGenes,
                    nCores=nCores)
saveRDS(fig.d,"fig.d.rds")
df_filtered <- fig.d[fig.d$Score != 0, ]
write.csv(df_filtered,"K562_filtered_network.csv",row.names=FALSE,quote=FALSE)
#fig.d
#rankDrivers(fig.d,rankBy = "meanScore")
#rankDrivers(fig.d,score.cut = 1.5,rankBy = "nTargets",interactive = TRUE)
#library(ComplexHeatmap)

#plotfigRHeatmap(figR.d = fig.d,
#                score.cut = 1.5,
#                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
#                show_row_dend = FALSE # from ComplexHeatmap
#)
#library(caTools)
#library(networkD3)
#plotfigRNetwork(fig.d,
#                score.cut = 1.5,
#                weight.edges = TRUE)


