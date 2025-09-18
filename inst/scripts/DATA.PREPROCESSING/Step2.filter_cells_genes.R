# Load the gene expression matrix from a TSV file
#gene_expression <- read.csv("RAW_DATA/scRNA/filtered_data_RNA_E7.5_rep1.csv", row.names = 1)
gene_expression <- read.csv("multiomic_data_filtered_L1_E7.5_rep1_RNA.csv", row.names = 1)
dim(gene_expression)
# Identify mitochondrial genes (assuming they start with "mt-" or "MT-")
mito_genes <- grep("^mt-|^MT-", rownames(gene_expression), value = TRUE)

# Remove mitochondrial genes from the dataset
gene_expression_filtered <- gene_expression[!rownames(gene_expression) %in% mito_genes, ]

# Check the dimensions of the filtered dataset
dim(gene_expression_filtered)


##################filter cells##############################

# Filter cells with fewer than 1000 genes expressed
# Calculate the number of expressed genes per cell
expressed_genes_per_cell <- colSums(gene_expression_filtered > 0)

# Select cells with at least 1000 genes expressed
cells_above_threshold <- names(expressed_genes_per_cell[expressed_genes_per_cell >= 1000])

# Filter the gene expression matrix to include only these cells
filtered_gene_expression <- gene_expression_filtered[, cells_above_threshold]
dim(filtered_gene_expression)
###################filter genes#################################
# Calculate the total number of cells
total_cells <- ncol(filtered_gene_expression)
total_cells
# Calculate the threshold for 10% of cells
threshold <- 0.1 * total_cells
threshold
# Calculate the number of non-zero counts for each gene
gene_counts <- rowSums(filtered_gene_expression > 0)
gene_counts
# Filter for genes expressed in at least 10% of cells
genes_above_threshold <- names(gene_counts[gene_counts >= threshold])
matrix_filtered <- filtered_gene_expression[genes_above_threshold,]
dim(matrix_filtered)
genes <- row.names(matrix_filtered)
write.csv(genes, file = "RNA.SC.DATA.PROCESSED/filtered_genes_multiomics_E7.5_rep1_L2_RNA.csv", row.names = FALSE)
write.csv(matrix_filtered, file = "RNA.SC.DATA.PROCESSED/filtered_multiomics_E7.5_rep1_L2_RNA.csv")

