# Load the RNA and ATAC datasets
#rna_data <- read.csv("RNA.SC.DATA.PROCESSED/Expression_Matrix_cells_1000g_and_genes_10prcnt_cells_E7.5_rep1.csv", row.names = 1)
#atac_data <- read.csv("RAW_DATA/scATAC/filtered_data_ATAC_E7.5_rep1.csv", row.names = 1)
rna_data <- read.csv("RNA.SC.DATA.PROCESSED/filtered_multiomics_E7.5_rep1_L2_RNA.csv", row.names = 1)
atac_data <- read.csv("multiomic_data_filtered_L1_E7.5_rep1_ATAC.csv", row.names = 1)
# Extract cell names (column names) from both datasets
rna_cells <- colnames(rna_data)
atac_cells <- colnames(atac_data)

# Identify the common cells
common_cells <- intersect(rna_cells, atac_cells)

# Filter both datasets to include only the common cells
rna_data_filtered <- rna_data[, common_cells]
atac_data_filtered <- atac_data[, common_cells]
write.csv(rna_data_filtered,file="RNA.SC.DATA.PROCESSED/filtered_multiomics_common_RNA_E7.5_rep1_L2.csv")
write.csv(atac_data_filtered,file="ATAC.SC.DATA.PROCESSED/filtered_multiomics_common_ATAC_E7.5_rep1_L2.csv")
# Check the dimensions to confirm the filtering
dim(rna_data_filtered)
dim(atac_data_filtered)
