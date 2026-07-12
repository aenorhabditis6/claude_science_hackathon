"""Load an AnnData and tag two conditions as the (reference, query) pair.

The one design invariant: a condition is just a categorical column in `.obs`, and a
comparison is two values of it. Nothing here knows about immunology.

`load()` works on any local `.h5ad`. `load_demo()` streams a subsample of the Marson
CD4+ T-cell Perturb-seq screen straight from public S3 so we never download the
120-170 GB per-donor files. Because the interesting cells are scattered across the file,
we pull their rows in parallel and cache the assembled slice to disk.
"""

import os

import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

# Marson genome-scale CD4+ T-cell Perturb-seq (GSE314342 / SRP643211).
# Processed cell-level h5ads live on the CZI Virtual Cells public S3 bucket, one file
# per donor (D1-D4) x culture condition (Rest / Stim8hr / Stim48hr). Anonymous access.
DEMO_BUCKET = "genome-scale-tcell-perturb-seq/marson2025_data"
DEMO_CONTROL = "NTC"  # what we rename non-targeting-guide cells to
CACHE_DIR = os.environ.get("SHIFTSCOPE_CACHE", os.path.expanduser("~/.shiftscope_cache"))

# Small public single-cell datasets, spanning domains, to show the tool is
# condition-agnostic: control vs stimulated, disease vs healthy, cell type vs cell type.
# The scverse-hosted h5ads are raw counts (clean differential expression).
PUBLIC_DATASETS = {
    "kang": dict(
        url="https://exampledata.scverse.org/pertpy/kang_2018.h5ad",
        file="kang_2018.h5ad", condition_key="label", celltype_key="cell_type",
        note="PBMC, control vs IFN-beta stimulated (Kang 2018). label: ctrl/stim.",
    ),
    "hagai": dict(
        url="https://exampledata.scverse.org/pertpy/hagai_2018.h5ad",
        file="hagai_2018.h5ad", condition_key="condition", species_key="species",
        note="Phagocytes across 4 species, unstimulated vs LPS (Hagai 2018). "
             "condition: unst/LPS6. Subset a species with within=('species','mouse').",
    ),
    "pbmc68k": dict(
        loader="pbmc68k_reduced", condition_key="bulk_labels",
        note="68k PBMC (scanpy built-in). Compare any two cell types in bulk_labels.",
    ),
    "dong": dict(
        url="https://exampledata.scverse.org/pertpy/dong_2023.h5ad",
        file="dong_2023.h5ad", condition_key="perturbation", celltype_key="cell_type0528",
        note="PBMC: No stimulation / IFNb / IFNg / IFNb+IFNg (Dong 2023). "
             "Multi-condition — rank them with compare.rank(adata, 'perturbation', 'No stimulation').",
    ),
    # --- the "any two populations, any domain" gallery: disease, development, drug target ---
    "haber": dict(
        url="https://exampledata.scverse.org/pertpy/haber_2017_regions.h5ad",
        file="haber_2017_regions.h5ad", condition_key="condition",
        note="DISEASE: mouse small-intestine epithelium, healthy vs parasite/bacterial "
             "infection (Haber 2017). condition: Control / Salmonella / Hpoly.Day3 / Hpoly.Day10. "
             "load_public('haber','Control','Hpoly.Day10') -> antimicrobial Reg3b/g, defensins.",
    ),
    "paul15": dict(
        loader="paul15", condition_key="paul15_clusters",
        note="DEVELOPMENT: mouse hematopoiesis (Paul 2015). paul15_clusters are differentiation "
             "states; compare two fates, e.g. load_public('paul15','14Mo','2Ery') -> erythroid "
             "master TF Klf1 + Car1/2, Ermap. (Cluster names are the paper's, e.g. Mo=monocyte, "
             "Ery=erythroid.)",
    ),
    "shifrut": dict(
        url="https://exampledata.scverse.org/pertpy/shifrut_2018.h5ad",
        file="shifrut_2018.h5ad", condition_key="perturbation", stim_key="perturbation_2",
        note="DRUG TARGET: human CD8 T cells, CRISPR KO of 20 genes +/- TCR stim (Shifrut 2018, "
             "Marson lab). Rank the knockouts under stimulation: filter within=('perturbation_2',"
             "'stim') then compare.rank(adata,'perturbation','control') -> LCP2/CD3D (TCR core) + "
             "RASA2/CBLB/SOCS1 (the negative regulators that paper found enhance T cells). "
             "perturbation_2: control/stim for a simple stim-vs-rest comparison.",
    ),
    "datlinger": dict(
        url="https://exampledata.scverse.org/pertpy/datlinger_2017.h5ad",
        file="datlinger_2017.h5ad", condition_key="perturbation_2",
        note="Jurkat T cells +/- anti-CD3/CD28 TCR stim (Datlinger 2017, CROP-seq). "
             "perturbation_2: unstimulated/stimulated -> activation markers CD69, IER3, BACH2.",
    ),
    "pbmc3k": dict(
        loader="pbmc3k_processed", condition_key="louvain",
        note="The classic 3k PBMC dataset (10x / Seurat tutorial). louvain: annotated cell "
             "types (CD4 T cells, CD14+ Monocytes, B cells, ...). "
             "load_public('pbmc3k','CD4 T cells','CD14+ Monocytes') -> myeloid TYROBP/S100A8.",
    ),
}


def tag_groups(adata, condition_key, reference, query):
    """Subset to the two named conditions and tag them in `.obs['_group']`.

    Returns a fresh AnnData containing only the reference and query cells, with
    `_group` a category ordered [reference, query].
    """
    col = adata.obs[condition_key].astype(str)
    keep = col.isin([str(reference), str(query)])
    if keep.sum() == 0:
        raise ValueError(
            f"No cells match {reference!r} or {query!r} in obs[{condition_key!r}]. "
            f"Seen values include: {sorted(col.unique())[:10]}"
        )
    out = adata[keep.values].copy()
    grp = out.obs[condition_key].astype(str)
    out.obs["_group"] = np.where(grp == str(reference), reference, query)
    out.obs["_group"] = out.obs["_group"].astype("category").cat.set_categories(
        [reference, query]
    )
    out.uns["_reference"] = reference
    out.uns["_query"] = query
    out.uns["_condition_key"] = condition_key
    return out


def _preprocess_if_needed(adata):
    """Normalize + log1p + HVG, but only if the matrix still looks like raw counts."""
    X = adata.X
    sample = X[:100].data if sp.issparse(X) else np.asarray(X[:100]).ravel()
    looks_raw = sample.size == 0 or np.allclose(sample, np.round(sample))
    if looks_raw:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    # Stash the full log-normalized matrix so drivers.py can run DE on all genes later,
    # even after we subset to HVGs and scale for the embedding (standard scanpy pattern).
    adata.raw = adata
    if adata.n_vars > 2000:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
        adata = adata[:, adata.var["highly_variable"]].copy()
    return adata


def load(path, condition_key, reference, query, preprocess=True):
    """Load a local `.h5ad`, tag the (reference, query) pair, and lightly preprocess."""
    adata = sc.read_h5ad(path)
    adata = tag_groups(adata, condition_key, reference, query)
    if preprocess:
        adata = _preprocess_if_needed(adata)
    return adata


def _fetch_public(spec):
    """Return the AnnData for a PUBLIC_DATASETS entry, downloading + caching if needed."""
    if "loader" in spec:  # a scanpy built-in dataset
        return getattr(sc.datasets, spec["loader"])()
    import urllib.request
    import shutil
    d = os.path.join(CACHE_DIR, "public")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, spec["file"])
    if not os.path.exists(path):
        # A plain urlretrieve gets a 403 from the CDN (it blocks the default agent).
        req = urllib.request.Request(spec["url"], headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r, open(path, "wb") as out:
            shutil.copyfileobj(r, out)
    return sc.read_h5ad(path)


def load_public(name, reference, query, condition_key=None, within=None, preprocess=True):
    """Load one of the small public demo datasets and tag a (reference, query) pair.

    `name` is a key of PUBLIC_DATASETS (e.g. "kang", "hagai", "pbmc68k"). `within` is an
    optional (obs_column, value) to first restrict to one subpopulation — e.g.
    `within=("cell_type", "CD14+ Monocytes")` to compare ctrl vs stim inside monocytes.
    This is the same condition-agnostic contract as `load`, on ready-made datasets.
    """
    spec = PUBLIC_DATASETS[name]
    adata = _fetch_public(spec)
    if within is not None:
        col, val = within
        adata = adata[adata.obs[col].astype(str) == str(val)].copy()
    key = condition_key or spec["condition_key"]
    adata = tag_groups(adata, key, reference, query)
    if preprocess:
        adata = _preprocess_if_needed(adata)
    return adata


def _pull_rows(f, rows):
    """Read the selected CSR rows from the open remote h5ad handle.

    anndata's lazy `sparse_dataset` resolves each row to its real byte range via the HDF5
    chunk index (the chunks are scattered across the 140 GB file, so this must be
    index-aware, not linear) and fetches them. Rows are sorted so nearby cells share reads.
    Returns a scipy CSR matrix.
    """
    from anndata.io import sparse_dataset

    return sparse_dataset(f["X"])[np.sort(rows)]


def load_demo(
    query="SETDB1",
    condition="Rest",
    donor="D1",
    n_per_group=500,
    seed=0,
    preprocess=True,
    cache=True,
    verbose=True,
):
    """Stream a control-vs-one-perturbation slice of the Marson Perturb-seq screen.

    Opens the remote per-donor h5ad lazily (never downloading it whole), reads only the
    small `.obs` guide columns, picks non-targeting-control cells and cells carrying a
    guide against `query`, subsamples each group, and pulls just those rows of the sparse
    count matrix in parallel. Returns a normal in-memory AnnData tagged for ShiftScope,
    with `obs['perturbation']` in {"NTC", query}. Cached to disk so re-runs are instant.
    """
    import s3fs
    import h5py
    import pandas as pd
    from anndata.io import read_elem

    cache_path = os.path.join(
        CACHE_DIR, f"{donor}_{condition}_{query}_n{n_per_group}_s{seed}.h5ad"
    )
    if cache and os.path.exists(cache_path):
        adata = sc.read_h5ad(cache_path)
        return _preprocess_if_needed(adata) if preprocess else adata

    import time

    def say(msg, t):
        if verbose:
            print(f"[load_demo +{time.time() - t:5.0f}s] {msg}", flush=True)

    t0 = time.time()
    url = f"{DEMO_BUCKET}/{donor}_{condition}.assigned_guide.h5ad"
    fs = s3fs.S3FileSystem(anon=True)
    # readahead cache + small blocks: h5py does many scattered metadata reads, and the
    # default s3fs cache refetches whole large blocks, which is ~100x slower here.
    store = fs.open(url, "rb", cache_type="readahead", block_size=512 * 1024)
    f = h5py.File(store, "r")
    n_vars = int(f["X"].attrs["shape"][1])
    say(f"opened remote h5ad ({donor}/{condition}, {n_vars} genes)", t0)

    guide_type = pd.Series(np.asarray(read_elem(f["obs"]["guide_type"]))).astype(str)
    gene = pd.Series(np.asarray(read_elem(f["obs"]["perturbed_gene_name"]))).astype(str)
    if "low_quality" in f["obs"]:
        lq = np.asarray(read_elem(f["obs"]["low_quality"])).astype(bool)
    else:
        lq = np.zeros(len(gene), dtype=bool)
    say(f"read guide annotations ({len(gene)} cells)", t0)

    ntc_idx = np.where((guide_type.values == "non-targeting") & ~lq)[0]
    tgt_idx = np.where((gene.values == query) & ~lq)[0]
    if len(tgt_idx) == 0:
        top = gene[guide_type == "targeting"].value_counts().head(15).index.tolist()
        raise ValueError(f"No cells for perturbation {query!r}. Try one of: {top}")

    rng = np.random.default_rng(seed)
    pick = lambda idx: (
        idx if len(idx) <= n_per_group else rng.choice(idx, n_per_group, replace=False)
    )
    sel = np.sort(np.concatenate([pick(ntc_idx), pick(tgt_idx)]))

    gene_name = np.asarray(read_elem(f["var"]["gene_name"])).astype(str)
    say(f"streaming {len(sel)} cell rows...", t0)

    X = _pull_rows(f, sel)
    f.close()
    say("finished streaming cells", t0)

    # Synthesize cell ids from the row index; reading the real barcodes (`obs['_index']`)
    # costs ~48 MB over the slow link and we don't use them downstream.
    var = pd.DataFrame({"gene_name": gene_name}, index=gene_name)
    obs = pd.DataFrame(index=[f"{donor}_{condition}_{i}" for i in sel])
    obs["guide_type"] = guide_type.values[sel]
    obs["perturbation"] = np.where(
        guide_type.values[sel] == "non-targeting", DEMO_CONTROL, gene.values[sel]
    )
    obs["donor"] = donor
    obs["culture_condition"] = condition

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.uns["_dataset"] = f"Marson CD4+ T Perturb-seq (GSE314342) {donor}/{condition}"

    if cache:
        os.makedirs(CACHE_DIR, exist_ok=True)
        adata.write_h5ad(cache_path)

    adata = tag_groups(adata, "perturbation", DEMO_CONTROL, query)
    if preprocess:
        adata = _preprocess_if_needed(adata)
    return adata
