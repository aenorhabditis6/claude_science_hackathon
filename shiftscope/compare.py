"""How far apart are two cell populations in latent space, and is it real?

Three complementary distances between the reference and query point clouds, all
computed in the shared PCA space from embed.py:

  - E-distance (energy distance)  the scPerturb / pertpy standard for perturbation
                                  strength: 2*mean(d_XY) - mean(d_XX) - mean(d_YY) on
                                  squared-Euclidean distances. This is the headline.
  - Sinkhorn optimal transport    entropy-regularized Wasserstein; geometry-aware.
  - Maximum Mean Discrepancy      RBF-kernel two-sample distance.

Significance is the **E-test** (Peidli et al., scPerturb, Nat. Methods 2024): a label
permutation test on the E-distance. We pool the two groups, compute the pairwise distance
matrix once, then reshuffle labels many times and see how often a random split reaches the
observed E-distance. We subsample per group first so it stays fast for a live demo.

Refs: scPerturb (Peidli et al. 2024); pertpy `pt.tl.Distance` / `DistanceTest`.
"""

import numpy as np
import ot
import pandas as pd
from sklearn.metrics import pairwise_distances


def _subsample(Z, n, rng):
    if len(Z) <= n:
        return Z
    return Z[rng.choice(len(Z), n, replace=False)]


def _median_sigma(A, B):
    """Median-heuristic RBF bandwidth from a pooled pairwise-distance sample."""
    pool = np.vstack([A, B])
    m = min(len(pool), 500)
    idx = np.random.default_rng(0).choice(len(pool), m, replace=False)
    d = pairwise_distances(pool[idx])
    med = np.median(d[d > 0])
    return med if med > 0 else 1.0


def _edist_from_D(D, a, b):
    """E-distance given a precomputed pooled distance matrix and two index sets."""
    s_ab = D[np.ix_(a, b)].mean()
    s_aa = D[np.ix_(a, a)].mean()
    s_bb = D[np.ix_(b, b)].mean()
    return float(2 * s_ab - s_aa - s_bb)


def edistance(A, B):
    """scPerturb E-distance on squared-Euclidean distances (0 iff identical clouds)."""
    D = pairwise_distances(np.vstack([A, B]), metric="sqeuclidean")
    nA = len(A)
    return _edist_from_D(D, np.arange(nA), np.arange(nA, len(A) + len(B)))


def etest(A, B, n_perm=1000, seed=0):
    """The scPerturb E-test: permutation p-value on the E-distance.

    Returns (observed E-distance, p-value). The pooled distance matrix is computed once,
    so thousands of permutations are cheap.
    """
    pooled = np.vstack([A, B])
    nA = len(A)
    D = pairwise_distances(pooled, metric="sqeuclidean")
    idx = np.arange(len(pooled))
    observed = _edist_from_D(D, idx[:nA], idx[nA:])

    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(len(pooled))
        if _edist_from_D(D, perm[:nA], perm[nA:]) >= observed:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return observed, p


def sinkhorn_ot(A, B, reg=0.05):
    """Entropy-regularized OT distance between two equally-weighted point clouds.

    Cost is squared-Euclidean normalized to a max of 1, so `reg` is on a comparable
    scale across datasets. This is a regularized Wasserstein, reported as a complementary
    geometry-aware distance (not the significance test).
    """
    M = ot.dist(A, B, metric="sqeuclidean")
    M /= M.max() + 1e-12
    a = np.ones(len(A)) / len(A)
    b = np.ones(len(B)) / len(B)
    return float(ot.sinkhorn2(a, b, M, reg))


def mmd(A, B, sigma=None):
    """Squared Maximum Mean Discrepancy with an RBF kernel (biased estimator)."""
    if sigma is None:
        sigma = _median_sigma(A, B)
    g = 1.0 / (2 * sigma ** 2)
    Kxx = np.exp(-g * pairwise_distances(A, A, metric="sqeuclidean"))
    Kyy = np.exp(-g * pairwise_distances(B, B, metric="sqeuclidean"))
    Kxy = np.exp(-g * pairwise_distances(A, B, metric="sqeuclidean"))
    return float(Kxx.mean() + Kyy.mean() - 2 * Kxy.mean())


def distances(Z_ref, Z_query, n_sub=500, n_perm=1000, reg=0.05, seed=0):
    """Distances between the two populations + an E-test p-value.

    n_sub caps cells per group (for speed); n_perm is the number of label shuffles.
    The p-value is (#permuted E-distances >= observed + 1) / (n_perm + 1).
    """
    rng = np.random.default_rng(seed)
    A = _subsample(Z_ref, n_sub, rng)
    B = _subsample(Z_query, n_sub, rng)

    edist, p = etest(A, B, n_perm=n_perm, seed=seed)
    return {
        "edistance": edist,
        "p_value": p,
        "sinkhorn": sinkhorn_ot(A, B, reg=reg),
        "mmd": mmd(A, B),
        "n_ref": int(len(Z_ref)),
        "n_query": int(len(Z_query)),
        "n_sub": int(min(n_sub, len(Z_ref), len(Z_query))),
        "n_perm": int(n_perm),
    }


def rank(adata, condition_key, reference, n_pcs=30, n_sub=500, n_perm=200, seed=0):
    """Rank every non-reference condition by how far it shifts the population.

    For each other value of `obs[condition_key]`, measures the E-distance (+ E-test p and
    MMD) from the reference in one shared PCA latent fit over all conditions. Returns a
    DataFrame sorted by E-distance — the multi-condition view: rank every perturbation in a
    screen, or every stimulus/disease/age group, by effect size. Expects a preprocessed
    (log-normalized) AnnData.
    """
    import scanpy as sc

    ad = adata.copy()
    sc.pp.scale(ad, max_value=10)
    if not np.all(np.isfinite(ad.X)):
        ad.X = np.nan_to_num(ad.X)
    sc.tl.pca(ad, n_comps=min(n_pcs, min(ad.shape) - 1),
              svd_solver="arpack", random_state=seed)
    Z = ad.obsm["X_pca"]
    labels = ad.obs[condition_key].astype(str).values
    Z_ref = Z[labels == str(reference)]
    if len(Z_ref) == 0:
        raise ValueError(f"reference {reference!r} not in obs[{condition_key!r}].")

    rows = []
    for q in sorted(set(labels) - {str(reference)}):
        Z_q = Z[labels == q]
        if len(Z_q) < 10:
            continue
        d = distances(Z_ref, Z_q, n_sub=n_sub, n_perm=n_perm, seed=seed)
        rows.append({"query": q, "edistance": d["edistance"], "p_value": d["p_value"],
                     "mmd": d["mmd"], "n_query": int(len(Z_q))})
    return pd.DataFrame(rows).sort_values("edistance", ascending=False).reset_index(drop=True)
