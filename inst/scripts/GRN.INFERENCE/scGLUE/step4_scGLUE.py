import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue
import seaborn as sns
from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout

scglue.plot.set_publication_params()
rcParams['figure.figsize'] = (4, 4)

rna = ad.read_h5ad("rna-emb.h5ad")
atac = ad.read_h5ad("atac-emb.h5ad")
guidance_hvf = nx.read_graphml("guidance-hvf.graphml.gz")

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

#genes = rna.var.query("highly_variable").index
#peaks = atac.var.query("highly_variable").index

genes = rna.var_names
peaks = atac.var_names


features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])

skeleton = guidance_hvf.edge_subgraph(
    e for e, attr in dict(guidance_hvf.edges).items()
    if attr["type"] == "fwd"
).copy()

reginf = scglue.genomics.regulatory_inference(
    features, feature_embeddings,
    skeleton=skeleton, random_state=0
)

gene2peak = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
    if attr["qval"] < 0.1
)
#default 0.05
scglue.genomics.Bed(atac.var).write_bed("peaks.bed", ncols=3)
scglue.genomics.write_links(
    gene2peak,
    scglue.genomics.Bed(rna.var).strand_specific_start_site(),
    scglue.genomics.Bed(atac.var),
    "gene2peak.links", keep_attrs=["score"]
)

motif_bed = scglue.genomics.read_bed("JASPAR2022-mm10.bed.gz")
#motif_bed = scglue.genomics.read_bed("ENCODE-TF-ChIP-hg38.bed.gz")
print(motif_bed.head())

tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)
print(tfs.size)

rna[:, np.union1d(genes, tfs)].write_loom("rna.loom")
np.savetxt("tfs.txt", tfs, fmt="%s")

import os
print("running pyscenic")
os.system("pyscenic grn rna.loom tfs.txt \
    -o draft_grn.csv --seed 0 --num_workers 20 \
    --cell_id_attribute obs_names --gene_attribute name")
print("test1")
peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
    gene2peak, peak2tf, genes, peaks, tfs,
    region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
    random_state=0
)
print(gene2tf_rank_glue.iloc[:5, :5])

flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(500, 500)
flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)

gene2flank = nx.Graph([(g, g) for g in genes])
gene2tf_rank_supp = scglue.genomics.cis_regulatory_ranking(
    gene2flank, flank2tf, genes, genes, tfs,
    n_samples=0
)
print(gene2tf_rank_supp.iloc[:5, :5])

gene2tf_rank_glue.columns = gene2tf_rank_glue.columns + "_glue"
gene2tf_rank_supp.columns = gene2tf_rank_supp.columns + "_supp"
gene2tf_rank_glue.to_csv("gene2tf_rank_glue.csv", index=True)
gene2tf_rank_supp.to_csv("gene2tf_rank_supp.csv", index=True)
scglue.genomics.write_scenic_feather(gene2tf_rank_glue, "glue.genes_vs_tracks.rankings.feather")
scglue.genomics.write_scenic_feather(gene2tf_rank_supp, "supp.genes_vs_tracks.rankings.feather")

pd.concat([
    pd.DataFrame({
        "#motif_id": tfs + "_glue",
        "gene_name": tfs
    }),
    pd.DataFrame({
        "#motif_id": tfs + "_supp",
        "gene_name": tfs
    })
]).assign(
    motif_similarity_qvalue=0.0,
    orthologous_identity=1.0,
    description="placeholder"
).to_csv("ctx_annotation.tsv", sep="\t", index=False)

os.system("pyscenic ctx draft_grn.csv \
    glue.genes_vs_tracks.rankings.feather \
    supp.genes_vs_tracks.rankings.feather \
    --annotations_fname ctx_annotation.tsv \
    --expression_mtx_fname rna.loom \
    --output pruned_grn.csv \
    --rank_threshold 5000 --min_genes 1 \
    --num_workers 32 \
    --cell_id_attribute obs_names --gene_attribute name 2> /dev/null")

grn = scglue.genomics.read_ctx_grn("pruned_grn.csv")
