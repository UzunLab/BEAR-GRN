# Load required data
RNA <- read.csv("RNA.SC.DATA.PROCESSED/filtered_common_RNA_E7.5_rep1.csv", header = TRUE, row.names = 1)
ATAC <- read.csv("ATAC.SC.DATA.PROCESSED/filtered_common_ATAC_E7.5_rep1.csv", header = TRUE, row.names = 1)

# Load barcodes
#barcodes <- readLines("/gpfs/Home/kmk7420/Multi_omics_GRN/CellOracle/mESC_DS014/intersection/common_cells_E7.5_rep1.txt")
barcodes <- colnames(RNA)
# Verify that barcodes were read correctly
print(typeof(barcodes))  # Should be 'character'
print(length(barcodes))  # Check how many barcodes are loaded
gene_names <- rownames(RNA)
peak_names <- rownames(ATAC)
print(gene_names)
print(peak_names)
# Initialize an empty list to store updated barcodes
barcode_vec <- list()
subsampled_RNA_list <- list()
subsampled_ATAC_list <- list()

# Define sizes for subsampling
sizes <- c(1000, 2000, 3000, 4000, 5000)

# Loop over different sizes to sample barcodes
for (cell_num in sizes) {
  
  # Randomly sample 'cell_num' barcodes
  set.seed(123)  # Set seed for reproducibility
  sampled_barcodes <- sample(barcodes, size = cell_num, replace = FALSE)
  
  # Check if the barcodes are sampled correctly
  print(paste("Sampled", cell_num, "barcodes"))
  print(head(sampled_barcodes))  # Print first few barcodes for confirmation
  
  # Clean barcodes: Remove quotes and replace `#` and `-` with `.`
  cleaned_sampled_barcodes <- gsub('"', '', sampled_barcodes)
  cleaned_sampled_barcodes <- gsub("#", ".", cleaned_sampled_barcodes)
  updated_barcodes <- gsub("-", ".", cleaned_sampled_barcodes)
  
  # Add "genes" as the first entry for both RNA and ATAC matrices
  #updated_barcodes <- c("genes", cleaned_sampled_barcodes)
  
  # Check if the barcodes are cleaned and updated correctly
  print(paste("Updated barcodes for cell_num:", cell_num))
  print(head(updated_barcodes))  # Print the first few updated barcodes
  
  # Subset RNA and ATAC data with sampled barcodes
  subsamples_RNA <- RNA[, updated_barcodes, drop = FALSE]
  subsamples_RNA <- data.frame(genes=rownames(RNA), subsamples_RNA, stringsAsFactors = FALSE)
  subsamples_ATAC <- ATAC[, updated_barcodes, drop = FALSE]
  subsamples_ATAC <- data.frame(genes=rownames(ATAC), subsamples_ATAC, stringsAsFactors = FALSE)  
  # Check dimensions of the subset matrices
  print(dim(subsamples_RNA))
  print(dim(subsamples_ATAC))
  
  # Optionally save the results to CSV files (commented out for now)
  # write.csv(subsamples_RNA, paste0("E7.5_rep1_", cell_num, "_cells_RNA.csv"), row.names = FALSE)
  # write.csv(subsamples_ATAC, paste0("E7.5_rep1_", cell_num, "_cells_ATAC.csv"), row.names = FALSE)
  
  # Store the updated barcodes in the list
  barcode_vec[[as.character(cell_num)]] <- updated_barcodes
  subsampled_RNA_list[[as.character(cell_num)]] <- subsamples_RNA
  subsampled_ATAC_list[[as.character(cell_num)]] <- subsamples_ATAC

}

# Final barcode list contains barcodes for each size in `sizes`
head(barcode_vec[['2000']])
head(subsampled_RNA_list[["2000"]])
list_data=subsampled_ATAC_list[["2000"]]
write.csv(list_data, "E7.5_rep1_2000_cells_ATAC.csv", row.names = FALSE)
saveRDS(barcode_vec, file = "RNA.SC.DATA.PROCESSED/E7.5_rep1_multiomics_common_subsampled_barcodes_list.rds")
# Save the subsampled_RNA_list as an RDS file
saveRDS(subsampled_RNA_list, file = "RNA.SC.DATA.PROCESSED/E7.5_rep1_multiomics_common_subsampled_RNA_list.rds")

# Save the subsampled_ATAC_list as an RDS file (if needed)
saveRDS(subsampled_ATAC_list, file = "ATAC.SC.DATA.PROCESSED/E7.5_rep1_multiomics_common_subsampled_ATAC_list.rds")
