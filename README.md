# ShiftScope

**Compare any two single-cell populations, measure how the cell-state distribution
shifts, localize where, and let Claude interpret the biology.**

Built for the "Built with Claude: Life Sciences" hackathon (Anthropic × Gladstone
Institutes). ShiftScope is condition-agnostic: give it one `.h5ad` with a categorical
column in `.obs` and a `(reference, query)` pair, and it tells you *how far apart* the
two populations are, *where* in cell-state space they diverge, *which genes* drive it,
and *what it might mean*. Control vs. perturbation, disease vs. healthy, young vs. old —
same tool.

## How it works

| Module | Job |
|--------|-----|
| `io.py` | Tag two conditions as `(reference, query)` in any `.h5ad`. `load_public()` fetches ready-made demo datasets across domains; `load_demo()` streams the 140 GB Marson Perturb-seq screen from S3 without downloading it. |
| `embed.py` | One PCA fit on the pooled cells (shared latent) + a 2D UMAP. |
| `compare.py` | E-distance (scPerturb) with an E-test permutation p-value, plus MMD and Sinkhorn-OT. `rank()` orders *many* conditions by how far each shifts the population. |
| `localize.py` | Leiden-cluster the latent; per-cluster compositional shift test (Fisher + BH) → clusters that gained/lost query cells. |
| `drivers.py` | Wilcoxon DE (scanpy) inside the shifted clusters → top driver genes. |
| `interpret.py` | Claude turns the numbers into a cited cell-state → pathway → target rationale; logs the exact prompt + evidence. |
| `prioritize.py` | Rank screen hits by **strong phenotype × under-studied**: real PubMed paper counts ground the novelty axis, then Claude scores each hit and gives a validate/skip verdict → a "validate these" shortlist + 2-D figure. |
| `calibration.py` | Stress-test the detection limits on a known shift: E-test **power vs cell number** and **vs effect size** → an honest "here's the operating envelope" figure. |
| `app.py` | Gradio UI (`launch(share=True)`): a **Compare** tab (pick any dataset + pair → distances, UMAP, drivers, Claude write-up) and a **Prioritize** tab (rank the Marson screen → "validate these" shortlist + figure). |

## Methods & prior art

ShiftScope deliberately builds on the standard single-cell toolkit rather than reinventing it:

- **E-distance / E-test** — the perturbation-strength metric and permutation test from
  [scPerturb](https://www.nature.com/articles/s41592-024-02144-6) (Peidli et al., *Nat.
  Methods* 2024), also in [pertpy](https://pertpy.readthedocs.io/)'s `pt.tl.Distance`.
- **MMD & Sinkhorn OT** — RBF-kernel maximum mean discrepancy and entropy-regularized
  optimal transport ([POT](https://pythonot.github.io/)) as complementary distances.
- **Differential abundance** — Leiden clustering + a per-cluster proportion test; the
  heavier neighborhood-level standard is [Milo](https://www.nature.com/articles/s41587-021-01033-z)
  (Dann et al. 2022), whose author co-produced this dataset.
- **scanpy / anndata** for preprocessing, PCA, UMAP, Leiden, and Wilcoxon DE.

## Quickstart (runs in ~20 s, downloads a 38 MB dataset)

```bash
pip install -r requirements.txt
```

```python
from shiftscope import io, embed, compare, localize, drivers

# Control vs IFN-beta-stimulated monocytes (Kang 2018). Same tool, any condition.
adata = io.load_public("kang", "ctrl", "stim", within=("cell_type", "CD14+ Monocytes"))
Z_ref, Z_query = embed.embed(adata)
print(compare.distances(Z_ref, Z_query))     # E-distance, E-test p, MMD, Sinkhorn
print(drivers.drivers(adata)["up"]["gene"].head(10).tolist())
```

## Worked examples — real output

Because the tool is condition-agnostic, the *same* five lines run across **four different
fields** — immunology, host–pathogen disease, developmental biology, and cell-type ID:

| Domain | Comparison | dataset | E-distance | drivers ShiftScope surfaced |
|---|---|---|--:|---|
| **Immunology** | IFN-β stim vs ctrl, monocytes | Kang 2018 (human) | **216** | ISG15, IFI6, IFIT1/3, MX1, OAS1 — **interferon** response |
| **Immunology** | LPS vs unstim, phagocytes | Hagai 2018 (**mouse**) | **826** | Ccl3/4/5, Il1a, Nfkbia, Tlr2 — **NF-κB / inflammatory** |
| **Disease** | healthy vs *H. polygyrus* infection, gut | Haber 2017 (mouse) | **81** | Reg3b/g, Defa24/17, Spink4 — **antimicrobial** response |
| **Development** | monocyte vs erythroid progenitor | Paul 2015 (mouse) | **1309** | Klf1, Car1/2, Ermap, Blvrb — **erythroid master program** |
| **Cell types** | monocyte vs B cell | pbmc68k (human) | **288** | TYROBP, S100A8, FCN1 — **myeloid** identity |

Two things to notice: the E-distance is **graded** (1309 for distinct lineages, 88 for a
subtle infection response — a meaningful "how different are they" statistic), and the drivers
recover the **correct, distinct biology** every time — interferon for a cytokine, NF-κB for a
TLR ligand, antimicrobials for a parasite, an erythroid TF for a blood fate — across human and
mouse. Same five lines of code, four different fields.

### Rank many conditions by effect size

`compare.rank()` orders every condition by how far it moves the population — the core
Perturb-seq question ("which perturbation does the most?"), here on Dong 2023 (PBMCs):

```python
adata = io._preprocess_if_needed(io._fetch_public(io.PUBLIC_DATASETS["dong"]))
compare.rank(adata, "perturbation", "No stimulation")
#         query   edistance   p_value
#    IFNb+ IFNg       100.0     0.002     # combined stimulus shifts most
#          IFNb        73.4     0.002
#          IFNg        35.7     0.002     # type-II IFN weakest here
```

The same call ranks a **real CRISPR screen** — the Shifrut 2018 human CD8 T-cell screen (20
knockouts ± TCR stim, Marson lab) — and recovers the paper's biology: the TCR-signaling core
(**LCP2/SLP-76, CD3D**) shifts cells most, followed by the **negative regulators (TCEB2, CBLB,
RASA2, SOCS1)** it identified as *enhancing* T-cell function — the same class of immunotherapy
targets. No cell streaming; see the notebook.

### Prioritize hits to validate — strong phenotype × under-studied

Ranking by effect size isn't the real bottleneck; the top of any screen is crowded with
genes we already understand. The decision that costs bench time is *which strong hit is
worth chasing* — which means pairing effect size with **how under-studied** a gene is.

`prioritize.py` runs this on the **real Marson CD4+ T-cell screen** (`DE_stats.suppl_table.csv`,
11k knockdowns — no cell streaming needed). The under-studied axis is **grounded, not
guessed**: a live NCBI PubMed query counts papers on each gene in immune biology. Claude then
scores novelty and gives a verdict *with that real count in front of it*, and the exact prompt
+ inputs are logged so every call traces back to data.

```python
from shiftscope import prioritize
shortlist = prioritize.prioritize(condition="Stim48hr", top_n=20)  # needs ANTHROPIC_API_KEY
prioritize.plot(shortlist)   # phenotype (y) vs log PubMed papers (x); top-left = validate these
#       gene   phenotype   pubmed_count   priority_score   (well-known genes sink)
#     TADA2B      5260.0              1          +2.96      strong effect, ~unstudied → prioritize
#      ELOF1      4258.0              1          +1.38
#       VAV1      3575.0            564          -2.37      strong but textbook T-cell gene → skip
#      PTPRC      3226.0            789          -3.08      (CD45) famous → skip
```

The tool sinks the textbook genes (**VAV1**, **CD45/PTPRC**, hundreds of papers) and surfaces
strong-but-obscure chromatin regulators (**TADA2B**, **SGF29**, **ELOF1**, 0–4 papers) as the
shortlist. Grounded (real counts), agentic (Claude scores + reasons), and honest about it.

### Calibration — when does it work?

A distance-with-a-p-value is only useful if you know where it breaks. `calibration.py`
stress-tests the E-test on the *known* Kang IFN-β shift and reports **power** (fraction of
random draws with p<0.05) over many repeats — not one lucky call:

```python
from shiftscope import calibration
cal = calibration.calibrate(reps=20)          # ~20 s; loads Kang monocytes, runs both sweeps
print(cal["summary"])
calibration.plot(cal["cells"], cal["effect"])
# Cell number: on a subtle shift, reliable detection (power ≥ 80%) needs ~100 cells/group.
# Effect size: at n=200 the shift is still detected until ~90% of the query is reference
#              cells — the minimum detectable shift.
```

Two boundaries, honestly stated: power collapses below ~100 cells/group, and holds until the
shift is diluted to <~10% of its true size. That's the operating envelope a reviewer asks for.

## Datasets

`io.load_public(name, reference, query, ...)` ships small, ready-made datasets spanning
domains (see `io.PUBLIC_DATASETS`):

| name | domain | comparison |
|---|---|---|
| `kang` | immunology | PBMC, ctrl vs IFN-β |
| `hagai` | immunology | phagocytes ± LPS (4 species) |
| `dong` | immunology | PBMC, multi-IFN (rank) |
| `haber` | **disease** | gut epithelium, healthy vs *Salmonella* / *H. polygyrus* |
| `paul15` | **development** | hematopoiesis, any two differentiation states |
| `shifrut` | **drug target** | human T cells, 20-gene CRISPR screen ± TCR stim (rank) |
| `datlinger` | immunology | Jurkat T cells ± TCR stim (CROP-seq) |
| `pbmc68k` | cell types | compare any two annotated cell types |
| `pbmc3k` | cell types | the classic 3k PBMC (Seurat tutorial), by cell type |

Any local `.h5ad` works the same way via `io.load(path, condition_key, reference, query)`.

For scale, `io.load_demo()` targets the Marson lab **genome-scale CRISPRi Perturb-seq in
primary human CD4+ T cells** (GEO **GSE314342** / SRA **SRP643211**; ~22M cells). Each
per-donor file on the CZI public S3 bucket is 120–170 GB, so `load_demo()` never downloads
them — it opens the remote h5ad lazily and streams only the subsample it needs (NTC controls
+ one perturbed gene). Best run from Colab, where S3 latency is low.

## License

MIT.
