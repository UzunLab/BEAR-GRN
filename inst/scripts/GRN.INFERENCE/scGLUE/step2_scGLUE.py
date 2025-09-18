import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import numpy as np
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

rna = ad.read_h5ad("Chen-2019-RNA.h5ad")
print(rna.shape)
n_genes=rna.shape[1]
atac = ad.read_h5ad("Chen-2019-ATAC.h5ad")
print(atac)

print(rna.X, rna.X.data)

rna.layers["counts"] = rna.X.copy()

sc.pp.highly_variable_genes(rna, n_top_genes=n_genes, flavor="seurat_v3")

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=30, svd_solver="auto")

sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.umap(rna, color="cell_type",show=False)

print(atac.X, atac.X.data)
import scipy.sparse as sp
atac = atac[:, atac.X.sum(axis=0) > 0] 
#if sp.issparse(atac.X):
#    atac.X.data = np.nan_to_num(atac.X.data, nan=0)
sc.pp.normalize_total(atac)
sc.pp.log1p(atac)
#atac.X = np.nan_to_num(atac.X, nan=0)
scglue.data.lsi(atac, n_components=30, n_iter=10)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)
sc.pl.umap(atac, color="cell_type",show=False)
print(rna.var.head())

scglue.data.get_gene_annotation(
    rna, gtf="gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)
print(rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head())

print(atac.var_names[:5])

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
print(atac.var.head())
#rna.var.dropna(subset=["chrom", "chromStart", "chromEnd"], inplace=True)
valid_genes = rna.var.dropna(subset=["chrom", "chromStart", "chromEnd"]).index
rna = rna[:, valid_genes] 
#rna.var.fillna(0, inplace=True)
#rna.var["chromStart"].fillna(0, inplace=True)
#rna.var["chromEnd"].fillna(0, inplace=True)

guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
print(guidance)

scglue.graph.check_graph(guidance, [rna, atac])

print(atac.var.head())
print(rna.var.dtypes)

#rna.var["artif_dupl"] = rna.var["artif_dupl"].astype(str)

rna.write("rna-pp.h5ad", compression="gzip")
atac.write("atac-pp.h5ad", compression="gzip")
nx.write_graphml(guidance, "guidance.graphml.gz")
