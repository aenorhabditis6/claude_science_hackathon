"""Where in cell-state space did the shift happen?

Leiden-cluster the shared PCA embedding, then ask, per cluster, whether the query is
over- or under-represented relative to its global fraction. Each cluster gets a 2x2
Fisher exact test (query/reference x in-cluster/out) and a log2 enrichment; p-values are
BH-corrected across clusters. Clusters where the query significantly gained or lost cells
are the localized shift, and get flagged for coloring on the UMAP.

This is the interpretable, dependency-light take on differential abundance. The heavier
standard for KNN-neighborhood differential abundance is Milo (Dann et al. 2022); we may
swap it in later, but per-cluster proportions are enough to point drivers.py at the right
cells.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact, false_discovery_control


def localize(adata, resolution=1.0, alpha=0.05, min_log2=0.5, seed=0):
    """Cluster the latent and return a per-cluster table of query enrichment.

    Adds `obs['_cluster']` (Leiden) and `obs['_shift']` in {"gained", "lost", "ns"}.
    Returns a DataFrame sorted by absolute log2 enrichment, one row per cluster.
    """
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        adata, resolution=resolution, key_added="_cluster",
        random_state=seed, flavor="igraph", n_iterations=2, directed=False,
    )

    ref = adata.uns["_reference"]
    is_query = (adata.obs["_group"] != ref).values
    cluster = adata.obs["_cluster"].values
    n_query_total, n_ref_total = int(is_query.sum()), int((~is_query).sum())

    rows = []
    for cl in pd.unique(cluster):
        in_cl = cluster == cl
        q_in = int((in_cl & is_query).sum())
        r_in = int((in_cl & ~is_query).sum())
        q_out, r_out = n_query_total - q_in, n_ref_total - r_in
        # log2 of the query:reference odds inside vs outside the cluster.
        odds = ((q_in + 0.5) / (r_in + 0.5)) / ((q_out + 0.5) / (r_out + 0.5))
        _, p = fisher_exact([[q_in, r_in], [q_out, r_out]])
        rows.append({
            "cluster": cl,
            "n_ref": r_in,
            "n_query": q_in,
            "frac_query": q_in / max(in_cl.sum(), 1),
            "log2_enrichment": float(np.log2(odds)),
            "pval": p,
        })

    tab = pd.DataFrame(rows)
    tab["qval"] = false_discovery_control(tab["pval"].values, method="bh")
    sig = (tab["qval"] < alpha) & (tab["log2_enrichment"].abs() >= min_log2)
    tab["direction"] = np.where(
        ~sig, "ns", np.where(tab["log2_enrichment"] > 0, "gained", "lost")
    )
    tab = tab.sort_values("log2_enrichment", key=np.abs, ascending=False).reset_index(drop=True)

    shift_map = dict(zip(tab["cluster"], tab["direction"]))
    adata.obs["_shift"] = pd.Categorical(
        [shift_map[c] for c in cluster], categories=["gained", "lost", "ns"]
    )
    return tab


def shifted_clusters(tab, direction=None):
    """Cluster labels with a significant shift (optionally only 'gained' or 'lost')."""
    sig = tab[tab["direction"] != "ns"]
    if direction is not None:
        sig = sig[sig["direction"] == direction]
    return sig["cluster"].tolist()
