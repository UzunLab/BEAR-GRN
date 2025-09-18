import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import celloracle
print(sc.__version__)
print(celloracle.__version__)
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [6, 4.5]
import sys
from scipy import sparse

# Assuming you have loaded 'counts' as a sparse matrix (genes x cells)
# And 'metadata' is a DataFrame with cell metadata, indexed by cell IDs
input_file=sys.argv[1]
output=sys.argv[2]
#indir = 'scRNA'
# Read the CSV files
#counts = pd.read_csv(f"{indir}/E7.5_rep1_1000_cells_RNA.csv", index_col=0)
counts = pd.read_csv(input_file, index_col=0)
#counts = pd.read_csv("data/E7.5_rep1/output_filtered_L2_E7.5_rep1/subsampled_RNA_filtered_L2_E7.5_rep1.csv", index_col=0)
#metadata = pd.read_csv("metadata.csv", index_col=0)

# Create an AnnData object
adata = ad.AnnData(X=counts.values)
#adata.obs = metadata

# Set variable names (features)
adata.var_names = counts.columns
adata.obs_names = counts.index
num_genes = len(adata.obs_names)
# Check that metadata rows match the counts columns (cells)
#assert metadata.shape[0] == counts.shape[1], "Number of cells in metadata and counts matrix do not match!"

# Create AnnData object
adata = ad.AnnData(X=counts.T)  # Transpose to get cells x genes
print(adata.var_names)
# Set cell names from metadata (assuming metadata is correctly aligned with counts)
#adata.obs_names = metadata.index  # Cells
#adata.var_names = counts.index      # Genes (make sure 'gene_names' is defined)

# Add metadata to obs
#adata.obs = metadata
print(adata)
# Save to h5ad
#adata.write(f"{indir}/E7.5_rep1_1000_csv_to_anndata.h5ad")


indata = 'scRNA'
#adata = ad.read_h5ad(f'{indata}/E7.5_rep1_1000_csv_to_anndata.h5ad')

print(adata)
ct = np.random.choice(["E7.5_rep1"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
# Only consider genes with more than 1 count
sc.pp.filter_genes(adata, min_counts=1)

adata.X = adata.X.astype('float64')

# Normalize gene expression matrix with total UMI count per cell
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')
#adata.X = adata.X.astype('float64')
#Convert 'n_counts' column to float64
#adata.var['n_counts'] = adata.var['n_counts'].astype('float64')
filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='cell_ranger',
                                              n_top_genes=num_genes,
                                              log=False)

# Subset the genes
adata = adata[:, filter_result.gene_subset]

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)


# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()
print(adata)

# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Diffusion map
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

sc.tl.diffmap(adata)
# Calculate neihbors again based on diffusionmap
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

sc.tl.louvain(adata, resolution=0.8)

# PAGA graph construction
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, show=False)
sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
outdir='RNA_DS14'
adata.write_h5ad(f"{output}_processed.h5ad")
print(adata.obsm.keys())
print(adata)

