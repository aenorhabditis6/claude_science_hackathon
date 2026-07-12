"""One shared low-rank latent for both groups, plus a 2D UMAP for the visuals.

Every downstream module (compare, localize) works in this pooled PCA space so the two
groups are directly comparable. Only drivers.py goes back to raw genes.
"""

import numpy as np
import scanpy as sc


def embed(adata, n_pcs=30, n_neighbors=15, run_umap=True, seed=0):
    """Fit ONE PCA on the pooled cells; return (Z_ref, Z_query) and store the UMAP.

    Writes `obsm['X_pca']`, `obsm['X_umap']` back onto `adata`. Returns the PCA
    coordinates split by `_group` so compare.distances() can consume them directly.
    """
    sc.pp.scale(adata, max_value=10)
    # Genes with zero variance across the pooled cells become NaN after scaling; PCA
    # rejects NaN, so zero them out.
    if not np.all(np.isfinite(adata.X)):
        adata.X = np.nan_to_num(adata.X)
    n_pcs = min(n_pcs, min(adata.shape) - 1)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=seed)

    if run_umap:
        sc.pp.neighbors(adata, n_neighbors=min(n_neighbors, adata.n_obs - 1),
                        n_pcs=n_pcs, random_state=seed)
        sc.tl.umap(adata, random_state=seed)

    Z = adata.obsm["X_pca"]
    ref = adata.uns["_reference"]
    is_ref = (adata.obs["_group"] == ref).values
    Z_ref, Z_query = Z[is_ref], Z[~is_ref]
    return Z_ref, Z_query
