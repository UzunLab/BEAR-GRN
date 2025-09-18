import pandas as pd
import scanpy as sc
import anndata

# --- For scRNA-seq data ---

# 1. Read the gene-by-cell CSV file.
#    (Assuming the first column contains gene names, which become the row index.)
rna_df = pd.read_csv("subsampled_RNA_all_cells_E7.5_rep1.csv", index_col=0)

# 2. Transpose the dataframe so that rows represent cells and columns represent genes.
rna_data = rna_df.transpose()

# 3. Create an AnnData object.
rna = anndata.AnnData(rna_data)

# 4. Add observation (cell) annotations.
#    Here, all cells are of type 'buffer1'. You can also add a 'domain' field to help distinguish modalities.
rna.obs["cell_type"] = "buffer1"
rna.obs["domain"] = "scRNA"  # You can name this field as you prefer.

# 5. Save the AnnData object to an h5ad file.
rna.write("Chen-2019-RNA.h5ad")


# --- For scATAC-seq data ---

# 1. Read the peak-by-cell CSV file.
#    (Assuming the first column contains peak identifiers, which become the row index.)
atac_df = pd.read_csv("subsampled_ATAC_all_cells_E7.5_rep1.csv", index_col=0)

# 2. Transpose the dataframe so that rows represent cells and columns represent peaks.
atac_data = atac_df.transpose()

# 3. Create an AnnData object.
atac = anndata.AnnData(atac_data)

# 4. Add observation (cell) annotations.
atac.obs["cell_type"] = "buffer1"
atac.obs["domain"] = "scATAC"  # Again, name this field as you prefer.

# 5. Save the AnnData object to an h5ad file.
atac.write("Chen-2019-ATAC.h5ad")
