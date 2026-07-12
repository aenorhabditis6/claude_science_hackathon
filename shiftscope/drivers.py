"""Which genes drive the shift?

Wilcoxon rank-sum differential expression (scanpy's `rank_genes_groups`) between the
query and reference cells, run on the full log-normalized matrix stashed in `adata.raw`.
Optionally restricted to the shifted clusters from localize.py so we ask "what's different
*where* the composition changed" rather than averaging over the whole population.

This is the one module that goes back to raw genes rather than the latent embedding.
"""

import scanpy as sc


def drivers(adata, clusters=None, n_top=25):
    """Top up/down genes for query vs reference.

    If `clusters` is given, DE is computed only within those `_cluster` labels (typically
    the shifted clusters). Returns {"up", "down", "table"} as tidy DataFrames with
    log2 fold change, Wilcoxon score, and BH-adjusted p-value.
    """
    ref = adata.uns["_reference"]
    query = adata.uns["_query"]

    sub = adata
    if clusters is not None:
        sub = adata[adata.obs["_cluster"].isin(list(clusters))].copy()

    sc.tl.rank_genes_groups(
        sub, groupby="_group", groups=[query], reference=ref,
        method="wilcoxon", use_raw=True,
    )
    table = sc.get.rank_genes_groups_df(sub, group=query)
    table = table.rename(columns={
        "names": "gene", "logfoldchanges": "log2fc",
        "scores": "score", "pvals_adj": "padj", "pvals": "pval",
    })
    # Rank by the Wilcoxon score (z-statistic), not raw log2fc: sorting on fold change
    # alone surfaces low-expression noise genes with huge but meaningless ratios. Split
    # up/down by the *score* sign (higher/lower in query) rather than log2fc sign, which
    # is robust even when inputs arrive pre-scaled (log2fc degenerate). Keep it to
    # significant genes when any are.
    sig = table[table["padj"] < 0.05]
    if len(sig) == 0:
        sig = table
    up = (sig[sig["score"] > 0].sort_values("score", ascending=False)
          .head(n_top).reset_index(drop=True))
    down = (sig[sig["score"] < 0].sort_values("score", ascending=True)
            .head(n_top).reset_index(drop=True))
    return {"up": up, "down": down, "table": table}
