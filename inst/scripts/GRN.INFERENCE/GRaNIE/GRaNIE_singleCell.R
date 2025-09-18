library(Seurat)
library(Signac)
library(GRaNIE)
library(biomaRt)
library(readr)

# Set working directory
setwd("/gpfs/Home/kmk7420/Multi_omics_GRN/GRaNIE/K562")

rna_data <- read.table("subsampled_RNA_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
atac_data <- read.table("subsampled_ATAC_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")

# Merge RNA and ATAC
########################
pbmc <- CreateSeuratObject(counts = rna_data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
library(BSgenome.Hsapiens.UCSC.hg38)

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_data), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_data <- atac_data[as.vector(grange.use), ]
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

#frag.file <- "/Users/lihuazhang/Documents/DIRECT-NET/PBMC/Data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_data,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

# Preprocess RNA
DefaultAssay(pbmc) <- "RNA"
pbmc
pbmc <- SCTransform(pbmc, verbose = FALSE, return.only.var.genes = FALSE)
#seurat_rna <- SCTransform(seurat_rna, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50)

# Preprocess ATAC
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:50)

############################
# Weighted Nearest Neighbor integration
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
#pbmc@neighbors
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", dims = 1:50)
#pbmc@graphs
# Clustering
seurat_obj <- FindClusters(pbmc, graph.name = "wknn", resolution = 10)

# Generate pseudobulk samples
pseudobulk <- AggregateExpression(seurat_obj, assays = c("RNA", "ATAC"), slot = "counts", fun = mean)
#pseudobulk
countsRNA.df <- as.data.frame(pseudobulk$RNA)
countsPeaks.df <- as.data.frame(pseudobulk$ATAC)
#View(countsRNA.df)
# Convert gene symbols to Ensembl IDs
# Convert row names into a column
countsRNA.df$ENSEMBL <- rownames(countsRNA.df)
# Optionally, remove the old row names
countsRNA.df <- countsRNA.df[, c("ENSEMBL", colnames(countsRNA.df)[-ncol(countsRNA.df)])]

#rownames(countsRNA.df_with_header) <- NULL
countsPeaks.df$peakID <- rownames(countsPeaks.df)
# Optionally, remove the old row names
countsPeaks.df <- countsPeaks.df[, c("peakID", colnames(countsPeaks.df)[-ncol(countsPeaks.df)])]

###################################
library(biomaRt)

# Connect to the Ensembl database for human genes
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve mapping between gene symbols (external_gene_name) and Ensembl gene IDs
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = countsRNA.df$ENSEMBL,
  mart = ensembl_mart
)
# Create a named vector mapping gene symbols to Ensembl IDs
mapping_vector <- setNames(gene_mapping$ensembl_gene_id, gene_mapping$external_gene_name)

# Replace the ENSEMBL column by mapping each gene symbol to its Ensembl ID
countsRNA.df$ENSEMBL <- mapping_vector[countsRNA.df$ENSEMBL]
rownames(countsRNA.df) <- NULL
rownames(countsPeaks.df) <- NULL

genomeAssembly = "hg38"
objectMetadata.l = list(name = "K562", genomeAssembly = genomeAssembly)
GRN = initializeGRN(objectMetadata = objectMetadata.l, outputFolder = "./", genomeAssembly = genomeAssembly)

# Add data to GRaNIE
GRN = addData(GRN, counts_peaks = countsPeaks.df, normalization_peaks = "none", idColumn_peaks = "peakID",
              counts_rna = countsRNA.df, normalization_rna = "none", idColumn_RNA = "ENSEMBL", forceRerun = TRUE)

# Add TFBS using HOCOMOCO v12
tfbs_folder = "/gpfs/Home/kmk7420/Multi_omics_GRN/GRaNIE/H12INVIVO/"
GRN = addTFBS(GRN, motifFolder = tfbs_folder, TFs = "all", filesTFBSPattern = "_TFBS", fileEnding = ".bed.gz", forceRerun = TRUE)

# Overlap peaks with TFBS
GRN = overlapPeaksAndTFBS(GRN, nCores = 12, forceRerun = TRUE)

# Add TF-Peak and Peak-Gene connections using Spearman correlation
GRN = addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, connectionTypes = c("expression"), corMethod = "spearman", maxFDRToStore = 0.3, forceRerun = TRUE)
GRN = addConnections_peak_gene(GRN, corMethod = "spearman", promoterRange = 250000, TADs = NULL, nCores = 12, plotDiagnosticPlots = FALSE, forceRerun = TRUE)
#?filterGRNAndConnectGenes
# Filter and finalize GRN
GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.3, peak_gene.fdr.threshold = 0.3, peak_gene.fdr.method = "BH", gene.types = c("all"), forceRerun = TRUE)
GRN = add_TF_gene_correlation(GRN, corMethod = "spearman", nCores = 12, forceRerun = TRUE)
saveRDS(GRN, "final_GRN_singlecell.rds")
GRN_connections.all = getGRNConnections(GRN, type = "all.filtered", include_TF_gene_correlations = TRUE,
                                        include_geneMetadata = TRUE)
GRN_connections.all
write.csv(GRN_connections.all,"GRN_connections_all_unfiltered_scK562_all_column.csv",row.names=FALSE,quote=FALSE)
selected_col_GRN <- GRN_connections.all[,c(2,16,25,26)]
selected_col_GRN
#write.csv(GRN_connections.all,"GRN_connections_all_unfiltered_scmESC_all_column.csv", row.names = FALSE, quote=FALSE)
write.csv(selected_col_GRN,"GRN_connections_all_unfiltered_scK562_selected_column.csv",row.names = FALSE, quote=FALSE)

library(dplyr)

filtered_sorted_unique_GRN <- selected_col_GRN %>%
  filter(TF_gene.p_raw < 0.05) %>%   # Filter for significant interactions
  distinct(TF.name, gene.name, .keep_all = TRUE) %>%  # Remove duplicate TF-gene pairs
  arrange(desc(TF_gene.r))           # Sort by TF_gene.r in descending order

# View top results
filtered_sorted_unique_GRN
write.csv(filtered_sorted_unique_GRN,"GRN_connections_filtered_sorted_scK562_uniq.csv", row.names = FALSE, quote=FALSE)
