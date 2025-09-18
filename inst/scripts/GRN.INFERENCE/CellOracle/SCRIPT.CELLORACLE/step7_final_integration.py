# 0. Import

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

input_RNA=sys.argv[1]
input_ATAC=sys.argv[2]
output = sys.argv[3]

adata = sc.read_h5ad(input_RNA)

print(adata)

print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")


# Load TF info which was made from mouse cell atlas dataset.

base_GRN = pd.read_parquet(input_ATAC)
# Check data
print(base_GRN.head())

# Instantiate Oracle object
oracle = co.Oracle()

# Check data in anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
adata.X = adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="cell_type",
                                   embedding_name="X_draw_graph_fr")

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)


outdir = 'output_GRN/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

#outdir = 'output' 
# Save oracle object.
oracle.to_hdf5(f"{outdir}{output}.celloracle.oracle")

# Load file.
oracle = co.load_hdf5(f"{outdir}{output}.celloracle.oracle")

# Check clustering data
sc.pl.draw_graph(oracle.adata, color="cell_type")

links = oracle.get_links(cluster_name_for_GRN_unit="cell_type", alpha=10,
                         verbose_level=10)

print(links.links_dict.keys())

links.links_dict["E7.5_rep1"]

# Set cluster name
cluster = "E7.5_rep1"
#if not os.path.exists(outdir):
#    os.makedirs(outdir)

# Save as csv
links.links_dict[cluster].to_csv(f"{outdir}{output}_raw_GRN_for_{cluster}.csv")


links.filter_links(p=0.001, weight="coef_abs",threshold_number=None)

# Save Links object.
links.to_hdf5(file_path=f"{outdir}{output}_links.celloracle.links")
links.filtered_links[cluster].to_csv(f"{outdir}{output}_{cluster}_final_GRN.csv")

