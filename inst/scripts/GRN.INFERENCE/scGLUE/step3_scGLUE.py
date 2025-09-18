from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

rna = ad.read_h5ad("rna-pp.h5ad")
print(rna.shape)  # (cells, genes)
print(rna.var.shape)  # (features, metadata)

print(set(rna.var_names) - set(rna.var.index))  # Should be empty
print(set(rna.var.index) - set(rna.var_names))  # Should be empty

atac = ad.read_h5ad("atac-pp.h5ad")
guidance = nx.read_graphml("guidance.graphml.gz")
rna.X = rna.X.astype(int)  # Ensure integer type for the counts
atac.X = atac.X.astype(int)  # Same for ATAC data
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=False,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=False,
    use_rep="X_lsi"
)

#guidance_hvf = guidance.subgraph(chain(
#    rna.var.query("highly_variable").index,
#    atac.var.query("highly_variable").index
#)).copy()

guidance_hvf = guidance.subgraph(chain(
    rna.var_names,
    atac.var_names
)).copy()


glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)

glue.save("glue.dill")
# glue = scglue.models.load_model("glue.dill")
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
print(dx)

_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.65)

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

rna.write("rna-emb.h5ad", compression="gzip")
atac.write("atac-emb.h5ad", compression="gzip")
nx.write_graphml(guidance_hvf, "guidance-hvf.graphml.gz")
