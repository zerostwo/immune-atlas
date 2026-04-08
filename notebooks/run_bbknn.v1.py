import anndata
import scanpy as sc
import pandas as pd

adata = sc.read("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.scanpy.h5ad")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

g1s_genes = pd.read_csv("/home/sduan/projects/sahmi/data/g1s_genes.txt", header=None)[0].tolist()
g2m_genes = pd.read_csv("/home/sduan/projects/sahmi/data/g2m_genes.txt", header=None)[0].tolist()

sc.tl.score_genes_cell_cycle(adata, s_genes=g1s_genes, g2m_genes=g2m_genes)

adata.obs["cc_diff"] = adata.obs["S_score"] - adata.obs["G2M_score"]

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, subset=True, layer="counts", batch_key="batch")
sc.pp.regress_out(adata, ["nCount_RNA", "percent.mt", "S_score", "G2M_score"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

sc.external.pp.bbknn(adata, batch_key="batch")

sc.tl.umap(adata)

import numpy as np

resolutions = np.round(np.arange(0.1, 1.01, 0.1), 2)
resolutions

for res in resolutions:
    print(res)
    key_added = f'leiden_{res:.1f}'
    sc.tl.leiden(adata, resolution=res, flavor='igraph', n_iterations=2, key_added=key_added)

# sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=1)

adata.write_h5ad("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.scanpy.processed.v1.h5ad")
