#suppressWarnings(invisible(library(stream)))
library(stream2)
#packageVersion("Seurat")
setwd("~/Multi_omics_GRN/STREAM/")
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
########################
rna_counts <- read.table("subsampled_RNA_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
atac_counts <- read.table("subsampled_ATAC_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
#rna_counts
#################### Create Seurat object
library(Seurat)
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
library(EnsDb.Hsapiens.v86)
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg38"
annotations <- readRDS("EnsDb.Hsapiens.v86_annot.rds")
#frag.file <- "/Users/lihuazhang/Documents/DIRECT-NET/PBMC/Data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 41,
  annotation = annotations
)
chrom_assay
pbmc[["ATAC"]] <- chrom_assay


###############RNA analysis########################
DefaultAssay(pbmc) <- "RNA"
#library(Signac)
library(dplyr)
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

###################Integration analysis######################################
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
#> Performing TF-IDF normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunSVD(pbmc)
#> Running SVD
#> Scaling cell embeddings
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
#pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "SeuratProject"
pbmc@meta.data

#############

#pbmc1 <- qs::qread("tutorial_pbmc.qsave")
#pbmc1
#pbmc@meta.data
#library(SeuratObject)
#library(Seurat)
###########################################
#pbmc1 <- RegionStats(pbmc, peak.assay="peaks", genome = BSgenome.Hsapiens.UCSC.hg38)
#pbmc1@assays$peaks@meta.features
head(pbmc@meta.data)
#??run_stream
var_genes <- length(row.names(pbmc@assays$RNA))
var_peaks <- length(row.names(pbmc@assays$ATAC))
options(future.globals.maxSize = 1000 * 1024^2)
en.regs <- run_stream(obj = pbmc,
                      qubic.path = "/gpfs/Home/kmk7420/Multi_omics_GRN/STREAM/qubic2/qubic",
                      peak.assay = "ATAC",
                      var.genes = 9000,
                      top.peaks = 50000,
                      min.cells = 41,
                      org = "hg38",
                      c.cutoff = 1.0,
                      distance = 5e+05,
                      BlockOverlap = 0.30,
                      Extension = 1.0
)
en.regs
saveRDS(pbmc,"pbmc.rds")
saveRDS(en.regs,"en.regs.rds")
#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
library(Seurat)
qs::qsave(en.regs, "eRegulons1.qsave")
names(en.regs[[1]])
head(en.regs[[1]]$genes)
head(en.regs[[1]]$peaks)
head(en.regs[[1]]$links)
#message ("The running time of eRegulon prediction is: ", time1["elapsed"])
#pbmc@meta.data["predicted.id"] <- pbmc@meta.data$orig.ident
#head(pbmc@meta.data["predicted.id"])
#head(colnames(pbmc@meta.data))
pbmc$celltype <- pbmc$seurat_clusters
pbmc$celltype

###########################
get_en_GRNs_singletype <- function(obj, en.regs, celltype.name = "my_celltype") {
  # Extract all enhancer-GRN links
  links <- unlist(as(lapply(en.regs, function(y) {
    GenomicRanges::mcols(y$links)$TF <- y$TF
    y$links
  }), "GRangesList"))
  
  grn <- list(
    links = links,
    cell.type = celltype.name,
    cells = colnames(obj)
  )
  
  message("Constructed enhancer GRN for: ", celltype.name)
  return(grn)
}

time2 <- system.time( en.grns <- get_en_GRNs_singletype(obj = pbmc,
                                                 en.regs = en.regs) )
en.grns
gr <- en.grns$links
df <- as.data.frame(gr)
write.csv(df, file = "K562_eGRN_links_without_filter_ct.csv", row.names = FALSE)
##########################

time2 <- system.time( en.grns <- get_cts_en_GRNs(obj = pbmc, celltype = "celltype",
                                                 en.regs = en.regs, peak.assay = "ATAC",
                                                 rna.dims = 30, atac.dims = 30,
                                                 padj.cutoff = 0.05,
                                                 out.dir = "./") )
#deafault padj 0.05
qs::qsave(en.grns, "eGRNs.qsave")
names(en.grns[[1]])
en.grns
head(en.grns[[1]]$links)
en.grns[[1]]$links
message ("The running time of eGRN construction is: ", time2["elapsed"])
time3 <- system.time( cts.en.regs <- get_cts_en_regs(obj = pbmc, peak.assay = "ATAC", de.genes = NULL,
                                                     cts.en.grns = en.grns, out.dir = "./", celltype = "celltype",
                                                     min.pct = 0.25, logfc.threshold = 0.25, padj.cutoff = 0.05) )
qs::qsave(cts.en.regs, "cell-type-specific-eRegulons.qsave")
names(cts.en.regs[[1]])
head(cts.en.regs[[1]]$genes)
head(cts.en.regs[[1]]$enhancers)
head(cts.en.regs[[1]]$links)
write.csv(cts.en.regs[[1]]$links,"cell-type-specific-eRegulons.csv", row.names = FALSE, quote = FALSE)
message ("The running time of cell-type-specific eRegulon discovery is: ", time3["elapsed"])
#ZZcts.en.regs
#TF_target <- qs::qread("TF_target_pairs.qsave")
#eRegulon <-  qs::qread("eRegulons.qsave")
#eRegulon[[1]]
#TF_target1 <- qs::qread("TF_target_pairs1.qsave")
#TF_target1
