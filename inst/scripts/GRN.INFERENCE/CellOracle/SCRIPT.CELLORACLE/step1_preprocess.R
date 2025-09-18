library(cicero)
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
library(monocle3)
library(Matrix)
#setwd("~/Multi_omics_GRN/CellOracle/mESC_DS014/filtered_data//subsamples/cells_1000/scATAC")
setwd("~/Multi_omics_GRN/CellOracle/VERSION1.Macrophase_DS026/Macrophase_data")
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the input and output filenames
input_file <- args[1]  # First argument: input file
output_file <- args[2] # Second argument: output file
#input_file <- "data/E7.5_rep1/multiomic_data_1000_cells_E7.5_rep1.rds"  # First argument: input file
#output_file <- "out/1000_cells_E7.5_rep1" # Second argument: output file

# Print the input and output for confirmation
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
#indata <- readRDS(input_file)
#DefaultAssay(indata) <- "ATAC"
indata <- read.csv(input_file, header = TRUE, row.names = 1, sep = ",")
#indata <- read.csv("E7.5_rep1_1000_cells_ATAC.csv", header = TRUE, row.names = 1, sep = ",")
dim(indata)
head(indata)
#Extract cell info
cellinfo <- data.frame(cells = colnames(indata))
row.names(cellinfo) <- cellinfo$cells
head(cellinfo)


# Extract peak information from row names
peak_names <- row.names(indata)
peak_names
#peak_names <- indata[['genes']]
peakinfo <- data.frame(site_name = peak_names)
peakinfo <- cbind(do.call(rbind, strsplit(peak_names, "[:-]")), peakinfo)
names(peakinfo)[1:3] <- c("chr", "bp1", "bp2")
row.names(peakinfo) <- peakinfo$site_name

head(peakinfo)
#sparse_matrix <- indata@assays$ATAC$data
# Convert to a sparse matrix (dgCMatrix)
sparse_matrix <- as(as.matrix(indata), "dgCMatrix")

# Create the input_cds object
input_cds <- new_cell_data_set(
  expression_data = sparse_matrix,
  cell_metadata = cellinfo,
  gene_metadata = peakinfo
)


# Detect genes (or peaks in this case)
input_cds <- detect_genes(input_cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# View the summary of the input_cds object
input_cds

# Visualize peak_count_per_cell
#hist(Matrix::colSums(exprs(input_cds)))
# Filter cells by peak_count
# Please set an appropriate threshold values according to your data
max_count <- 50000
min_count <- 2000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count]
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count]
# Data preprocessing
set.seed(2017)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP


cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
class(input_cds)
umap_coords
# Save Cds object (Optional)
#output_folder = "ATAC_DS14/"
cicero_cds_file <- paste0(output_file, "_cicero_cds.Rds")
saveRDS(cicero_cds, cicero_cds_file)
#data("human.hg19.genome")
#chromosome_length <- human.hg19.genome

# If your scATAC-seq was aligned to the mm10 reference genome, you can read the chromosome length file using the following command.
#download.file(url = "https://raw.githubusercontent.com/morris-lab/CellOracle/master/docs/demo_data/mm10_chromosome_length.txt",
#              destfile = "./mm10_chromosome_length.txt")
chromosome_length <- read.table("hg38.fa.sizes")

# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

# Save results (Optional)
connfile <- paste0(output_file,"_cicero_connections.Rds")
saveRDS(conns, connfile)

# Check results
head(conns)
peakfile <- paste0(output_file,"_all_peaks.csv")
all_peaks <- row.names(exprs(input_cds))
all_peaks_cleaned <- gsub("[:\\-]", "_", all_peaks)
write.csv(x = all_peaks_cleaned, file = peakfile, row.names = FALSE)
cicero_conn <- paste0(output_file,"_cicero_connections.csv")
write.csv(x = conns, file = cicero_conn, row.names = FALSE)
