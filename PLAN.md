# ShiftScope — next-phase plan (from Noga + friend feedback)

Judging weights: **Demo 30% · Claude-use 25% · Impact 25% · Depth 20%.** These two features
target the highest-leverage criteria. Build #1 first (the differentiator), then #2 (rigor).

## Where we are (recap for a fresh session)
- Pipeline complete + **validated on real public data**: `kang` (IFN-β), `hagai` (mouse LPS),
  `pbmc68k` (cell types), `dong` (multi-IFN). Correct biology each time (ISGs / NF-κB / markers).
- `io.load_public(name, ref, query, within=)`, `compare.rank(adata, key, ref)` (ranks conditions
  by E-distance), `interpret.interpret()` (Anthropic SDK — needs `ANTHROPIC_API_KEY`).
- Env: `.venv` (Python 3.12). Demo data cached in `~/.shiftscope_cache/`. Nothing committed yet.
- Marson genome-scale hit table is downloadable WITHOUT streaming cells:
  `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/suppl_tables/DE_stats.suppl_table.csv`
  (per-gene `n_total_de_genes`, `ontarget_effect_size`, `culture_condition`) = real phenotype
  strength for every perturbation. Use this as the ranking source for the prioritization demo.

## Feature 1 — Hit prioritization: "strong phenotype × under-studied"  (the differentiator)
**Why:** measuring the shift isn't the bottleneck — deciding *which hits to validate* is. Turn
ShiftScope from a metric into a shortlist. Leans hard on Claude (25%) + Impact (25%).

New module `shiftscope/prioritize.py`:
1. **Phenotype strength** — from `compare.rank(...)` (or the Marson `DE_stats` table): E-distance /
   n_DE_genes per hit gene.
2. **Literature coverage (GROUNDING — this is Noga's rigor point)** — real signal, not LLM vibes.
   Query NCBI E-utilities PubMed esearch (free, no key):
   `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term="<gene>" AND (<context terms, e.g. "T cell" OR "immune")` → `esearchresult.count`. Low count + strong
   phenotype = under-studied + promising. Cache counts to disk; be polite (≤3 req/s, add `&tool=&email=`).
3. **Claude verdict + novelty score** — give Claude `{gene, phenotype rank/E-dist, pubmed_count,
   known pathway membership if available}` → returns novelty score (1–5) + one-line verdict
   ("strong + under-studied → prioritize" vs "strong but well-characterized" vs "weak"). Log prompt+inputs.
4. **Output**: table `[gene, phenotype, pubmed_count, novelty, verdict]` sorted to surface
   strong+novel. **Money figure**: 2-D scatter — x = literature coverage (log PubMed count),
   y = phenotype strength; top-LEFT quadrant = "validate these" (strong + novel). Label points.

Demo: run on Marson `DE_stats` hits (Stim48hr) — real T-cell perturbation hits, ranked and
novelty-scored. Add a notebook section + README subsection.

## Feature 2 — Calibration / operating boundaries  (Noga's rigor → Depth points)
**Why:** honestly state *when it works*. New `calibration.py` (or a notebook section) + one figure.
- **Sensitivity vs cell number**: subsample n = {25,50,100,200,400} per group on a known shift
  (Kang monocytes ctrl vs stim); plot E-distance + E-test p vs n. Where does power collapse?
- **Sensitivity vs effect size**: dilute the query group with a fraction f of reference cells
  (f = 0→0.9); the true shift shrinks. Plot E-distance + detection vs f → minimum detectable shift.
- Deliverable: figure + one sentence: "detects shifts down to ~N cells / ~M% mixing."

## Build order
1. `prioritize.py` + PubMed helper + Claude scoring; demo on Marson DE_stats table.
2. Prioritization money-figure (2×2 scatter).
3. `calibration.py` + figure.
4. Wire prioritize into `app.py` / `demo.ipynb`; update README.

## Watch-outs
- PubMed esearch: rate-limit, cache, URL-encode the term. Handle 0-count / network fail gracefully.
- Keep the novelty score GROUNDED in the real count (show the count next to Claude's score) — that's
  what makes it defensible vs. "an LLM guessed."
- `interpret`/prioritize need `ANTHROPIC_API_KEY`; degrade gracefully (show phenotype + PubMed count
  even without a key).
