# ShiftScope — Claude Code project guide

ShiftScope is a Builder-track project for the "Built with Claude: Life Sciences"
hackathon (Anthropic x Gladstone Institutes). It is a general tool that compares two
single-cell populations, measures and localizes how the cell-state distribution shifts,
and uses Claude to interpret the biology.

Deadline: Mon July 13, 9:00pm ET. All code is new work from this week. License: MIT.

## Golden rules
- Write simple, readable code. Short functions, clear names, one job per file.
- Minimal error handling. Don't wrap everything in try/except; let it fail loudly during
  the hackathon so bugs are obvious.
- No premature abstraction. No classes unless a function genuinely needs to hold state.
- Colab-first: everything must run top-to-bottom in `demo.ipynb`. The app uses gradio
  with `launch(share=True)`.
- Comment the *why*, not the *what*.

## The one design invariant
The tool is condition-agnostic. Input is always:
- one AnnData (`.h5ad`) with a categorical column in `.obs` (the "condition")
- a `(reference, query)` pair naming two values of that column

Nothing may be hardcoded to immunology or Perturb-seq. Marson is only the demo dataset;
swapping the column should let it do disease vs. healthy, young vs. old, etc.
Downstream modules operate on a shared low-rank latent embedding. Only `drivers.py`
touches raw genes.

## Modules (build in this order)
io.py -> embed.py -> compare.py -> localize.py -> drivers.py -> interpret.py -> app.py

- io.py — `load(path, condition_key, reference, query)` returns an AnnData tagged with
  `.obs["_group"]` in {reference, query}. `load_demo()` fetches the Marson CD4+ T-cell
  Perturb-seq set. Normalize + log1p + HVG only if the data isn't already processed.
- embed.py — `embed(adata)` fits ONE PCA on the pooled cells (~30-50 dims). Returns
  `Z_ref`, `Z_query`, and a 2D UMAP for the visuals.
- compare.py — `distances(Z_ref, Z_query)` returns Sinkhorn-OT, MMD, and energy distance
  with a permutation p-value (subsample per group for speed). `rank(adata, condition_key,
  reference)` loops every query vs. reference and returns a sorted table.
- localize.py — Leiden-cluster the shared latent, compare per-cluster group proportions,
  return the clusters that gained/lost cells (colorable on the UMAP).
- drivers.py — Wilcoxon differential expression inside the shifted clusters; return the
  top genes per shift.
- interpret.py — send `{shifted clusters, direction, top genes, distances}` to the Claude
  API; return a short, cited rationale (cell states -> pathway -> candidate target). Log
  the exact prompt + inputs so every claim traces back to data.
- app.py — Gradio UI: pick condition -> distances, UMAP with shifted clusters highlighted,
  driver panel, Claude write-up, one-click report. `launch(share=True)`.

## Data
- Marson CD4+ T-cell Perturb-seq. Find the REAL GEO accession from the lab's preprint or
  code first. Do not invent an accession or download URL.
- Assume the data may be raw counts; keep any preprocessing in io.py and keep it minimal.

## Commands
- Install: `pip install -r requirements.txt`
- Run pipeline: open `demo.ipynb` in Colab, Run all
- Launch app: `python app.py` (or the last cell of the notebook)

## Environment
- Claude API key in env var `ANTHROPIC_API_KEY`. Use the anthropic SDK with a current
  Claude model string (check https://docs.claude.com for the latest); Sonnet is fine for
  speed in interpret.py.
- Keep requirements.txt minimal: scanpy, anndata, pot, scikit-learn, umap-learn, leidenalg,
  gradio, anthropic, numpy, pandas.

## Judging priorities (optimize here)
- Demo 30% — the app must work live on Marson and look clean.
- Claude Use 25% — interpret.py is where we win; make it genuinely agentic and cited, not
  a one-shot summary.
- Impact 25% — lead with "compare any two cell populations."
- Depth 20% — the OT / MMD math is the moat; get it right.
