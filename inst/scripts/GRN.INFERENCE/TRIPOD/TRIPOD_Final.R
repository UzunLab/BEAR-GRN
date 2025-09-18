
# Install the remotes package
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}

# Insall Seurat and Seurat Object versions required for running TRIPOD
#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

# Install required dependency for Signac
#BiocManager::install("ggbio")
#BiocManager::install("GenomeInfoDb")
renv::snapshot()
# Install Signac ( this version works with Seurat v4)
#devtools::install_version(package = 'Signac', version = package_version('1.1.1'))
#install.packages("devtools")

# Install Seurat Disk
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
remotes::install_github("mojaveazure/seurat-disk")

# Intstall Other Packages needed
#BiocManager::install("packageName")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


setwd("~/Multi_omics_GRN/TRIPOD/mESC_data/test_TRIPOD/filtered/Env")
renv::activate()
renv::restore()
setwd("/gpfs/Home/kmk7420/Multi_omics_GRN/TRIPOD/mESC_data/test_TRIPOD/filtered/Env/test_tripod/k562/all_genes")


# Load Libraries

library(BiocParallel)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DescTools)


# Read in the data. The 10x hdf5 file contains both data types. 
#inputdata.10x <- Read10X_h5("../input_data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

# View data and modalities present
#str(inputdata.10x)

# extract RNA and ATAC data
#rna_counts <- inputdata.10x$`Gene Expression`
#atac_counts <- inputdata.10x$Peaks
rna_counts <- read.csv("subsampled_RNA_K562_human_filtered.csv", row.names = 1)
atac_counts <- read.csv("subsampled_ATAC_K562_human_filtered.csv", row.names = 1)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)

# QC by removing cells with higher mitochondiral counts
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# view the new variable added to the metadata
View(pbmc@meta.data)

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes

# Filter and Prepare ATAC Data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))

# Change to UCSC since the data was mapped to hg19
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get Gene Annotations for ATAC Peaks
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# Read in the fragment file
#frag.file <- "../input_data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
#View(frag.file)
# Create a chromatin assay (Add gene annonations to ATAC data)
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  annotation = annotations
)

View(atac_counts)

# Add ATAC to the Seurat Object created

pbmc[["ATAC"]] <- chrom_assay

# view the Updated Seurat Object
#View(pbmc@meta.data)


#VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
#        log = TRUE, pt.size = 0) + NoLegend()

# Now the pbmc object contains both RNA and ATAC data

# Perform QC on both assays extracting high quality cells for downstream analysis
#pbmc <- subset(
#  x = pbmc,
#  subset = nCount_ATAC < 7e4 &
#    nCount_ATAC > 5e3 &
#    nCount_RNA < 25000 &
#    nCount_RNA > 1000 &
#    percent.mt < 20
#)

# View the update Seurat object after subsetting
#View(pbmc)

# Randomly select 1000 cells
#set.seed(123)
#selected_cells <- sample(Cells(pbmc), size = 1000, replace = FALSE)
#pbmc <- subset(pbmc, cells = selected_cells)

#View(pbmc@meta.data)


# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',
          verbose = FALSE)



# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#View(pbmc@meta.data)


# # Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(pbmc) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
pwm_set
motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = TRUE)
motif.matrix
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)

# Know the number of motifs present 
motif.object # 633 motifs in 106056 regions

pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Set multiprocessing options before running chromvar (this is a must to do by the chromvar package developers )
# I used SerialParam because I use windows
# please check here for other options https://greenleaflab.github.io/chromVAR/articles/Introduction.html
register(SerialParam())
register(MulticoreParam())
# Run chromVAR - for sparse chromatin accessibility analysis
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

# Save chromVAR results
write.csv(pbmc@assays$chromvar@data, file = "chromVAR_Motif_Results.csv")

# Save seurat object after chromvar
saveRDS(pbmc, file = "pbmc_after_chromvar.rds")


# Cell Annotations - Mapping to a Multimodal reference (from CiteSeq)

# load reference data
#reference <- LoadH5Seurat("../input_data/pbmc_multimodal.h5seurat")


# Visualize the reference data
#p1=DimPlot(object = reference, 
#           reduction = "wnn.umap",
#           group.by = "celltype.l1", 
#           label = TRUE, 
#           label.size = 3, 
#           repel = TRUE) + NoLegend() + ggtitle("CITE-seq celltype l1")
#p2=DimPlot(object = reference, 
#           reduction = "wnn.umap",
#           group.by = "celltype.l2", 
#           label = TRUE, 
#           label.size = 3, 
#           repel = TRUE) + NoLegend() + ggtitle("CITE-seq celltype l2")
#p1+p2

# Find anchors between reference and query(in this case our seurat object, pbmc)
# using a precomputed supervised PCA (spca) transformation for this example. 
#DefaultAssay(pbmc) <- "SCT"
#anchors <- FindTransferAnchors(
#  reference = reference,
#  query = pbmc,
#  normalization.method = "SCT",
#  reference.reduction = "spca",
#  dims = 1:50,
#  recompute.residuals = FALSE
#)

# Transfer cell type labels and protein data from the reference to the query
# project the query data onto the umap structure of the reference
#pbmc <- MapQuery(
#  anchorset = anchors,
#  query = pbmc,
#  reference = reference,
#  refdata = list(
#    celltype.l1 = "celltype.l1",
#    celltype.l2 = "celltype.l2",
#    predicted_ADT = "ADT"
#  ),
#  reference.reduction = "spca", 
#  reduction.model = "wnn.umap"
#)

#pbmc$celltype=as.factor(pbmc$predicted.celltype.l2)
#pbmc$celltype="Macrophage"

# Remove some minor cell types (with less than 100 cells)
#pbmc$celltype[pbmc$celltype=='cDC2']='cDC'
#pbmc=pbmc[,!pbmc$celltype%in%names(which(table(pbmc$celltype)<100))]
#pbmc$celltype=droplevels(pbmc$celltype)

# Renormalize after QC
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', verbose = FALSE)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# For each cluster, we only keep the main cell type
cell.keep=rep(FALSE, ncol(pbmc))
pbmc$celltype = "Macrophage"
temp=table(pbmc$seurat_clusters, pbmc$celltype)
temp
print('test1')
for(i in 1:nrow(temp)){
  clusteri=as.numeric(rownames(temp)[i])
  celltypei=colnames(temp)[which.max(temp[i,])]
  cell.keep[which(pbmc$seurat_clusters==clusteri & pbmc$celltype==celltypei)]=TRUE
}
pbmc=pbmc[,cell.keep]
pbmc
# Renormalize after QC
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', 
          verbose = FALSE)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
#p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


# Save seurat object after Preprocessing
saveRDS(pbmc, file = "pbmc.rds")



###################################################################
#-------------Preparing the Seurat Object fo TRIPOD -----------------------#
###################################################################

# Install TRIPOD and load library
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#devtools::install_github("yuchaojiang/TRIPOD/package")

library(TRIPOD)
library(dendextend)

# Creating objects required for model fitting
# this will contain a list of object with informattion necessary for the fitting
# this only keeps motifs when genes encoding corresponding TFs are present in the RNA-seq data.
# this ensures that the same set of genes is included in all objects.
tripod.pbmc <- getObjectsForModelFit(object = pbmc, chr = paste0("chr", 1:22))

# Extract  the inidvidual object from the model
transcripts.gr <- tripod.pbmc$transcripts.gr
peaks.gr <- tripod.pbmc$peaks.gr
motifxTF <- tripod.pbmc$motifxTF
peakxmotif <- tripod.pbmc$peakxmotif


# Filter the seurat object 
# This exclude cases where a TF corresponds to multiple motifs or vice versa
# Well the TRIPOD authors hope to include such in a future version of the package
pbmc <- filterSeuratObject(object = pbmc, tripod.object = tripod.pbmc)

# perform normalization, dimension reduction, nearest neighbor graph construction, 
# and Uniform Manifold Approximation and Projection (UMAP) embedding for the filtered data
pbmc <- processSeuratObject(object = pbmc, dim.rna = 1:50, dim.atac = 2:50,
                            verbose = FALSE)


# Visualize the data on recalculated graphs
p1 <- DimPlot(pbmc, reduction = "umap.rna",  group.by = "celltype", 
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac",  group.by = "celltype", 
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", 
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("WNN")
#p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


# Partition the single cells into clusters to account for the issue of sparsity

# First generete a data frame containing the resolutions, number of resultant clusters
# and the number of clusters with fewer than threshold number of single cells
# The threshold is set to 20 cells
set.seed(123)
num.clusters <- optimizeResolution(
  object = pbmc,
  graph.name = "wsnn",
  assay.name = "WNN",
  resolutions = seq(10, 35, 5),
  min.num = 20
)

# Check the number of clusters
num.clusters

# Resolution 15 was selected to get the clusters
res <- 15
pbmc <- getClusters(object = pbmc, graph.name = "wsnn", algorithm = 3,
                    resolution = res, verbose = FALSE)


# obtain matrices containing normalized RNA expression and chromatin accessibility per metacell (referred to as cluster)
metacell.pbmc <- getMetacellMatrices(object = pbmc,
                                     cluster.name = "seurat_clusters"
)

# Extract Individual Matrices
metacell.rna <- metacell.pbmc$rna
metacell.peak <- metacell.pbmc$peak
View(metacell.rna)
View(metacell.peak)
# Save the Seurat Object for TRIPOD Analysis
saveRDS(pbmc, file = 'pbmc_metacells.rds')


# Assign colors to cell types and create objects storing mapping between them
pbmc$celltype=factor(pbmc$celltype)
#color.pbmc <- getColors(object = pbmc)

# Extracts the mapping between the metacells and cell types
metacell.celltype <- pbmc$celltype
#metacell.celltype.col <- color.pbmc$metacell$color

# obtain highly variable genes based on SCT
DefaultAssay(pbmc) <- "SCT"
hvg.pbmc <- VariableFeatures(pbmc)
#rownames(pbmc)
# Save Intermediate files for TRIPOD
pbmc.transcripts.gr=transcripts.gr
save(pbmc.transcripts.gr, file = 'pbmc.transcripts.gr.rda')
pbmc.peaks.gr=peaks.gr
save(pbmc.peaks.gr, file='pbmc.peaks.gr.rda')
pbmc.motifxTF=motifxTF
save(pbmc.motifxTF, file='pbmc.motifxTF.rda')
pbmc.peakxmotif=peakxmotif
save(pbmc.peakxmotif, file='pbmc.peakxmotif.rda')
pbmc.metacell.rna=metacell.rna
save(pbmc.metacell.rna, file='pbmc.metacell.rna.rda')
pbmc.metacell.peak=metacell.peak
save(pbmc.metacell.peak, file='pbmc.metacell.peak.rda')
pbmc.metacell.celltype=metacell.celltype
save(pbmc.metacell.celltype, file='pbmc.metacell.celltype.rda')
#pbmc.metacell.celltype.col=metacell.celltype.col
#save(pbmc.metacell.celltype.col, file='../input_data/pbmc.metacell.celltype.col.rda')


##############################################################################
#--------- Detecting Candidate Trio Regulatory Relationships---------------------------
###############################################################################

# Fit model
hvg.pbmc
# Usually over 1000 target genes are used for the model fitting. 
# I used just two because I was getting issues with parallelization
#genes <- c("CCR7", "GNLY")
#genes <- c("MYO16","HBZ","SLC25A37")
#genes <- hvg.pbmc
genes <- rownames(pbmc)
genes
# set the size of the window around transcription start sites (TSS) to 100 kb.
# For a given gene, ATAC peaks within the windows are considered as possible cis-regulatory regions.
ext.upstream <- ext.downstream <- 1e5

# NB// For a given trio consisting of the gene g the peak t and the TF j,
# TRIPOD characterizes the target gene expression (denoted as Yg) as a function of the chromatin accessibility 
# in the peak region (denoted as Xt) and the TF (denoted as Yj) expression. The Xt variables are extracted from 
# the ATAC-seq data, whereas the Yg and Yj variables are extracted from the RNA-seq data.

#We first invoke the getXYMatrices() function. 
# For a given gene, the function returns a list containing a set of Yg, Xt and Yj which they call a xymats list.
# Here, we obtain a list of xymats objects for the two genes.
bpparam <- BiocParallel::bpparam()
#bpparam$verbose <- TRUE

xymats.list <- bplapply(
  genes,
  getXYMatrices,
  ext.upstream = ext.upstream,
  transcripts.gr = pbmc.transcripts.gr,
  peaks.gr = pbmc.peaks.gr,
  metacell.rna = pbmc.metacell.rna,
  metacell.peak = pbmc.metacell.peak,
  peakxmotif = pbmc.peakxmotif,
  motifxTF = pbmc.motifxTF,
  metacell.celltype = pbmc.metacell.celltype
)
names(xymats.list) <- genes

#xymatlist_new$Xt
#xymats.list$MYO16$Xt
# Filter xymats.list to retain genes where Xt is a 2D matrix
filtered_xymats <- xymats.list[sapply(xymats.list, function(gene_data) {
  # Check if "Xt" exists and is a 2D matrix
  "Xt" %in% names(gene_data) && is.matrix(gene_data$Xt) && ncol(gene_data$Xt) > 1
})]
genes <- names(filtered_xymats)
length(genes)
# run TRIPOD matching Xt (Running TRIPOD by characterizing target gene expression as a function of chromatin accessibility)
xymats.tripod.Xt.list <- bplapply(
  filtered_xymats,
  fitModel,
  model.name = "TRIPOD",
  match.by = "Xt",
  BPPARAM = bpparam
)

names(xymats.tripod.Xt.list) <- genes
xymats.tripod.Xt.list
# run TRIPOD matching Yj (Running TRIPOD by characterizing target gene expression as a function of TF expression)
xymats.tripod.Yj.list <- bplapply(
  filtered_xymats,
  fitModel,
  model.name = "TRIPOD",
  match.by = "Yj"
)
names(xymats.tripod.Yj.list) <- genes


## Exploring the 2 Conditional levels for each of the 2 models above susggested in the TRIPOD paper
# Get list of hits by Benjamini-Hochberg (BH) method with a certain FDR threshold. We typically use a threshold of 0.01 
# and focus on hits with significantly positive coefficients.

# set FDR < 0.01
fdr.thresh <- 1

# focus on positive sign
sign <- "positive"


# NB// The getTrios() function returns a data frame with columns containing the target gene name, peak number, 
# TF number, peak name, TF name, estimates of the coefficient, nominal p-value, and BH-adjusted p-value. 

# TRIPOD level 1 matching Xt
xymats.tX1.pos.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 1
)
head(xymats.tX1.pos.df[order(xymats.tX1.pos.df$adj), ])
write.csv(xymats.tX1.pos.df,"xymats.tX1.pos.df.csv")
# TRIPOD level 2 matching Xt
xymats.tX2.pos.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 2
)
head(xymats.tX2.pos.df[order(xymats.tX2.pos.df$adj), ])
write.csv(xymats.tX2.pos.df,"xymats.tX2.pos.df.csv")

# TRIPOD level 1 matching Yj
xymats.tY1.pos.df <- getTrios(
  xymats.list = xymats.tripod.Yj.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 1
)
head(xymats.tY1.pos.df[order(xymats.tY1.pos.df$adj), ])
write.csv(xymats.tY1.pos.df,"xymats.tY1.pos.df.csv")
# TRIPOD level 2 matching Xt
xymats.tY2.pos.df <- getTrios(
  xymats.list = xymats.tripod.Yj.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 2
)
head(xymats.tY2.pos.df[order(xymats.tY2.pos.df$adj), ])
write.csv(xymats.tY2.pos.df,"xymats.tY2.pos.df.csv")



# focus on negative sign
sign <- "negative"


# NB// The getTrios() function returns a data frame with columns containing the target gene name, peak number,
# TF number, peak name, TF name, estimates of the coefficient, nominal p-value, and BH-adjusted p-value.

# TRIPOD level 1 matching Xt
xymats.tX1.neg.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 1
)
head(xymats.tX1.neg.df[order(xymats.tX1.neg.df$adj), ])
write.csv(xymats.tX1.neg.df,"xymats.tX1.neg.df.csv")
# TRIPOD level 2 matching Xt
xymats.tX2.neg.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 2
)
head(xymats.tX2.neg.df[order(xymats.tX2.neg.df$adj), ])
write.csv(xymats.tX2.neg.df,"xymats.tX2.neg.df.csv")

# TRIPOD level 1 matching Yj
xymats.tY1.neg.df <- getTrios(
  xymats.list = xymats.tripod.Yj.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 1
)
head(xymats.tY1.neg.df[order(xymats.tY1.neg.df$adj), ])
write.csv(xymats.tY1.neg.df,"xymats.tY1.neg.df.csv")
# TRIPOD level 2 matching Xt
xymats.tY2.neg.df <- getTrios(
  xymats.list = xymats.tripod.Yj.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 2
)
head(xymats.tY2.neg.df[order(xymats.tY2.neg.df$adj), ])
write.csv(xymats.tY2.neg.df,"xymats.tY2.neg.df.csv")

