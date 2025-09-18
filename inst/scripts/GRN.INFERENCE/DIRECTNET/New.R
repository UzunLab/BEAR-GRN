setwd("~/Multi_omics_GRN/DIRECTNET/K562")
library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

##########load data#############################
#rna_counts <- inputdata.10x$`Gene Expression`
#atac_counts <- inputdata.10x$Peaks
rna_counts <- read.table("subsampled_RNA_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")
atac_counts <- read.table("subsampled_ATAC_K562_human_filtered.csv", header = TRUE, row.names = 1, sep = ",")

genome.info <- read.table(file = "hg38.TSS.regions.txt", header = TRUE)
#genome.info
#names(genome.info) <- c("Chrom","Starts","Ends","genes")
#genes <- lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
#genes <- lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
#genes <- unlist(genes)
#genome.info$genes <- genes
unik <- !duplicated(genome.info$genes)# filter out different transcript
#unik
genome.info <- genome.info[unik,]
#genome.info
#################### Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
print("test1")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
print("test2")
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
print("test3")
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

#frag.file <- "/Users/lihuazhang/Documents/DIRECT-NET/PBMC/Data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

#VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
#        log = TRUE, pt.size = 0) + NoLegend()

# pbmc <- subset(
#   x = pbmc,
#   subset = nCount_ATAC < 7e4 &
#     nCount_ATAC > 5e3 &
#     nCount_RNA < 25000 &
#     nCount_RNA > 1000 &
#     percent.mt < 20
# )

###############RNA analysis########################
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

###################Integration analysis######################################
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
#> Performing TF-IDF normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 1)
pbmc <- RunSVD(pbmc)
#> Running SVD
#> Scaling cell embeddings
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
#pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "Buffer1"
#Idents(pbmc) <- "sub.cluster"
Idents(pbmc)
# add annotations
# pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
# pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
# pbmc <- RenameIdents(pbmc, '17' = 'Naive B', '11' = 'Intermediate B', '10' = 'Memory B', '21' = 'Plasma')
# pbmc <- RenameIdents(pbmc, '7' = 'NK')
# pbmc <- RenameIdents(pbmc, '4' = 'CD4 TEM', '13'= "CD4 TCM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
# pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2')
# pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
# pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')

pbmc$celltype <- Idents(pbmc) 
# p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
# p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
DefaultAssay(pbmc) <- "RNA"

#pbmc1 <- FindVariableFeatures(pbmc, assay = "RNA", selection.method = "vst", nfeatures = 2000)
#pbmc1
#markers <- pbmc1@assays$RNA@var.features
#markers
markers <- row.names(pbmc)
print("Markers")
markers
pbmc <- Run_DIRECT_NET(pbmc, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, 
                       max_overlap=0.5, size_factor_normalize = FALSE, 
                       genome.info = genome.info, focus_markers = markers)
direct.net_result <- Misc(pbmc, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind,direct.net_result)) # links for markers
direct.net_result

direct.net_result$function_type <- gsub("HF","HC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("Rest","MC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("LF","LC",direct.net_result$function_type)
##############Visualize the links####################################
temp <- tempfile()
download.file("https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#marker <- "ENO1"
#Plot_connections(direct.net_result, gene_anno, marker, cutoff = 0.5, upstream = 100000, downstream = 10000)

###############Construct regulatory networks##################################
### identify differential accessible peaks (DA)
DefaultAssay(pbmc) <- 'ATAC'
# #focused_markers <- markers[which(markers$group %in% c("CD4 TEM", "CD4 TCM")), , drop = FALSE]
# groups <- unique(focused_markers$group)
# groups  
# da_peaks <- FindMarkers(
#   object = pbmc,
#   min.pct = 0.2,
#   logfc.threshold = 0.6,
#   ident.1 = groups,
#   group.by = "celltype",
#   test.use = 'LR',
#   only.pos = TRUE
# )
# pbmc@meta.data$celltype
pbmc
#pbmc2 <- FindTopFeatures(pbmc, min.cutoff = 100)  # Adjust cutoff as needed
#pbmc2
# Extract the top variable peaks
variable_peaks <- VariableFeatures(pbmc)
variable_peaks
#######################Detect CRE-TF connections###############################
direct.net_result
# Assuming 'markers' is your list/vector of genes
focused_markers <- data.frame(
  gene = markers,  # Column for gene names
  group = rep("Buffer1", length(markers))  # Assign "Buffer1" to all genes
)
focused_markers
# CRE-gene connections
CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = focused_markers)

CREs_Gene$promoter
# Find focused CREs which is overlapped with DA
variable_peaks1 <- list()
variable_peaks1$Buffer1 <- variable_peaks
variable_peaks1
Focused_CREs <- generate_CRE(L_G_record = CREs_Gene$distal, P_L_G_record = CREs_Gene$promoter, variable_peaks1)
Focused_CREs
# detect TFs for distal CREs
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install("JASPAR2016")
L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$distal, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)
L_TF_record
# detect TFs for Promoters
P_L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$promoter, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = focused_markers)
P_L_TF_record
groups <- pbmc$celltype
groups
network_links <- generate_links_for_Cytoscape(L_G_record = Focused_CREs$L_G_record, L_TF_record, P_L_G_record = Focused_CREs$P_L_G_record, P_L_TF_record,groups)
#??generate_links_for_Cytoscape
network_links
write.csv(network_links,"Network_links.csv", row.names = FALSE, quote = FALSE)
Node_attribute <- generate_node_for_Cytoscape(network_links,markers = focused_markers)
Node_attribute
