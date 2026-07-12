"""Turn a list of screen hits into a ranked "validate these" shortlist.

Measuring the shift isn't the bottleneck — deciding *which* hits to chase is. A hit is
worth a bench scientist's time when it has a **strong phenotype** but is **under-studied**:
a big effect that nobody has explained yet. This module scores hits on exactly those two
axes and asks Claude for a validate/skip verdict.

The under-studied axis is GROUNDED, not vibes: we query NCBI PubMed (E-utilities esearch,
free, no key) for how many papers mention the gene in the relevant biological context, and
show that real count next to Claude's score. Low count + strong phenotype = prioritize.

Condition-agnostic like the rest of ShiftScope: `phenotype_strength` can come from
`compare.rank(...)` (E-distance per condition) or from any per-hit table. The demo feeds it
the real Marson CD4+ T-cell Perturb-seq hit table (`DE_stats.suppl_table.csv`), where
phenotype strength is the number of downstream DE genes each knockdown causes.

Pipeline: `marson_hits()` -> `add_pubmed_counts()` -> `score_hits()` (Claude) -> `plot()`.
Or just call `prioritize()` to run all four. Degrades gracefully with no `ANTHROPIC_API_KEY`
(you still get phenotype + PubMed counts + the figure).
"""

import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

CACHE_DIR = os.environ.get("SHIFTSCOPE_CACHE", os.path.expanduser("~/.shiftscope_cache"))

# Real per-gene hit strengths from the Marson genome-scale CD4+ T-cell Perturb-seq screen.
# One row per (gene x culture_condition); n_total_de_genes = downstream DE genes a knockdown
# causes = phenotype strength. No cell streaming needed to rank hits.
DE_STATS_URL = (
    "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
    "marson2025_data/suppl_tables/DE_stats.suppl_table.csv"
)

# PubMed context for the Marson screen: how studied is this gene *in T-cell / immune biology*.
# Swap these for the biology at hand (e.g. ("neuron", "brain") for a neural screen).
IMMUNE_CONTEXT = ("T cell", "immune", "lymphocyte")


# --------------------------------------------------------------------------- hit tables

def marson_hits(condition="Stim48hr", require_ontarget=True, url=DE_STATS_URL):
    """Load the Marson DE_stats table as a ranked hit list for one culture condition.

    Returns a DataFrame `[gene, phenotype, ontarget_effect_size, n_cells]` sorted by
    phenotype (downstream DE genes), strongest first. With `require_ontarget` (default),
    keeps only knockdowns that actually worked (`ontarget_significant`) so the phenotype is
    attributable to the perturbation rather than to noise or an off-target guide.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, "marson_DE_stats.csv")
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r, open(path, "wb") as out:
            out.write(r.read())
    df = pd.read_csv(path)

    df = df[df["culture_condition"] == condition]
    if require_ontarget:
        df = df[df["ontarget_significant"] == True]  # noqa: E712 (pandas mask, not `is`)
    hits = pd.DataFrame({
        "gene": df["target_contrast_gene_name"].astype(str).values,
        "phenotype": df["n_total_de_genes"].astype(float).values,
        "ontarget_effect_size": df["ontarget_effect_size"].astype(float).values,
        "n_cells": df["n_cells_target"].astype(float).values,
    })
    return hits.sort_values("phenotype", ascending=False).reset_index(drop=True)


def from_rank(rank_table):
    """Adapt a `compare.rank(...)` result into the hit-list shape this module expects.

    So the same prioritization runs on any ShiftScope comparison, not just Marson: each
    condition's E-distance becomes its phenotype strength.
    """
    return pd.DataFrame({
        "gene": rank_table["query"].astype(str).values,
        "phenotype": rank_table["edistance"].astype(float).values,
    }).sort_values("phenotype", ascending=False).reset_index(drop=True)


# --------------------------------------------------------- literature coverage (grounding)

def _pubmed_cache():
    path = os.path.join(CACHE_DIR, "pubmed_counts.json")
    if os.path.exists(path):
        with open(path) as f:
            return path, json.load(f)
    return path, {}


def pubmed_count(gene, context=IMMUNE_CONTEXT, email=None, cache=None, pause=0.34):
    """How many PubMed papers mention `gene` in this biological context. None on failure.

    Real grounding for the "under-studied" axis. Queries NCBI E-utilities esearch (free, no
    key) and returns `esearchresult.count`. Results are cached to disk, so a gene is only
    ever fetched once. `pause` keeps us under NCBI's ~3 req/s courtesy limit.
    """
    ctx_key = "|".join(context)
    key = f"{gene}||{ctx_key}"
    own_cache = cache is None
    cache_path, store = _pubmed_cache() if own_cache else (None, cache)
    if key in store:
        return store[key]

    ctx = " OR ".join(f'"{t}"' for t in context)
    term = f'"{gene}" AND ({ctx})'
    params = {"db": "pubmed", "retmode": "json", "term": term, "tool": "shiftscope"}
    if email:
        params["email"] = email
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + \
        urllib.parse.urlencode(params)

    count = None
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "shiftscope"})
        with urllib.request.urlopen(req, timeout=20) as r:
            count = int(json.load(r)["esearchresult"]["count"])
    except Exception:
        count = None  # network / parse failure: leave unknown, don't crash the demo
    time.sleep(pause)

    if count is not None:
        store[key] = count
        if own_cache:
            with open(cache_path, "w") as f:
                json.dump(store, f)
    return count


def add_pubmed_counts(hits, context=IMMUNE_CONTEXT, email=None, top_n=None, verbose=True):
    """Add a `pubmed_count` column to a hit table (top `top_n` rows if given).

    One shared on-disk cache for the whole batch, written once at the end so a run of N
    genes does at most N esearch calls ever.
    """
    df = (hits if top_n is None else hits.head(top_n)).copy()
    cache_path, store = _pubmed_cache()
    counts = []
    for i, gene in enumerate(df["gene"], 1):
        counts.append(pubmed_count(gene, context=context, email=email, cache=store))
        if verbose:
            print(f"  [pubmed {i}/{len(df)}] {gene}: {counts[-1]}", flush=True)
    with open(cache_path, "w") as f:
        json.dump(store, f)
    df["pubmed_count"] = counts
    return df


# ------------------------------------------------------------------ Claude novelty verdict

SYSTEM = """You are a computational biologist triaging hits from a genetic screen for \
experimental follow-up. For each hit you are given its phenotype strength (how large an \
effect knocking the gene out had) and a real PubMed paper count measuring how well studied \
the gene already is in the relevant biology.

You are prioritizing for VALIDATION, so the prize is a strong phenotype on an UNDER-STUDIED \
gene: a big effect nobody has explained yet. A strong effect on a famous, well-characterized \
gene is lower priority (likely already known); a weak effect is low priority regardless.

Rules:
- Ground the novelty score in the PubMed count you were given; a low count is the evidence \
of "under-studied". Do not invent facts about specific papers.
- Do not claim a gene's function unless it is broadly known; if unsure, say so briefly.
- Reply with ONLY a JSON array, one object per gene, no prose or code fences:
  [{"gene": str, "novelty": int 1-5 (5 = strong phenotype + very under-studied), \
"recommendation": "prioritize" | "maybe" | "skip", "verdict": one short sentence}]"""


def _extract_json(text):
    """Pull the JSON array out of Claude's reply, tolerating fences, prose, and truncation.

    If the array parses cleanly, return it. If it doesn't (e.g. the response hit the token
    cap mid-array), salvage every complete `{...}` object we can — so a long shortlist still
    yields verdicts for the hits Claude did finish rather than crashing the demo.
    """
    s = text.strip()
    if s.startswith("```"):
        s = s.split("```")[1]
        if s.lstrip().startswith("json"):
            s = s.lstrip()[4:]
    a, b = s.find("["), s.rfind("]")
    if a != -1 and b > a:
        try:
            return json.loads(s[a:b + 1])
        except json.JSONDecodeError:
            pass
    objs, dec, i = [], json.JSONDecoder(), (a + 1 if a != -1 else 0)
    while i < len(s):
        j = s.find("{", i)
        if j == -1:
            break
        try:
            obj, i = dec.raw_decode(s, j)
            objs.append(obj)
        except json.JSONDecodeError:
            i = j + 1
    return objs


def score_hits(hits, context=IMMUNE_CONTEXT, model="claude-sonnet-5",
               log_dir="interpretations", max_tokens=4000):
    """Ask Claude for a novelty score + validate/skip verdict per hit. GROUNDED in PubMed.

    Expects a `pubmed_count` column (from `add_pubmed_counts`). Sends the whole shortlist in
    one call and logs the exact prompt + inputs, so every score traces back to the phenotype
    and paper count we handed over. If there's no `ANTHROPIC_API_KEY`, returns the table
    unchanged with empty score columns so the rest of the demo still runs.
    """
    df = hits.copy()
    if not os.environ.get("ANTHROPIC_API_KEY"):
        df["novelty"] = np.nan
        df["recommendation"] = ""
        df["verdict"] = "(no ANTHROPIC_API_KEY — showing phenotype + PubMed only)"
        return df

    import anthropic

    evidence = [
        {
            "gene": r.gene,
            "phenotype_strength": round(float(r.phenotype), 3),
            "phenotype_rank": i + 1,
            "pubmed_count": (None if pd.isna(r.pubmed_count) else int(r.pubmed_count)),
        }
        for i, r in enumerate(df.itertuples())
    ]
    user = (
        f"Context for the PubMed counts: papers mentioning the gene together with "
        f"{' / '.join(context)}. Phenotype strength is the effect size of the knockout "
        f"(higher = stronger), ranked 1 = strongest. Score these {len(evidence)} hits:\n\n"
        "```json\n" + json.dumps(evidence, indent=2) + "\n```"
    )

    client = anthropic.Anthropic()
    resp = client.messages.create(
        model=model, max_tokens=max_tokens, system=SYSTEM,
        messages=[{"role": "user", "content": user}],
    )
    text = "".join(b.text for b in resp.content if b.type == "text")
    scores = {s["gene"]: s for s in _extract_json(text)}

    df["novelty"] = [scores.get(g, {}).get("novelty", np.nan) for g in df["gene"]]
    df["recommendation"] = [scores.get(g, {}).get("recommendation", "") for g in df["gene"]]
    df["verdict"] = [scores.get(g, {}).get("verdict", "") for g in df["gene"]]

    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        path = os.path.join(log_dir, f"prioritize_{int(time.time())}.json")
        with open(path, "w") as f:
            json.dump({"model": model, "system": SYSTEM, "user_prompt": user,
                       "evidence": evidence, "output": text,
                       "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}, f, indent=2)
    return df


# ----------------------------------------------------------------------- priority + figure

def _priority_score(df):
    """Combine the two axes: strong phenotype + few papers = high priority. Grounded, no LLM.

    z(phenotype) - z(log10 papers): rises with effect size, falls with how well-studied.
    """
    ph = df["phenotype"].to_numpy(float)
    pc = np.log10(df["pubmed_count"].to_numpy(float) + 1)
    z = lambda v: (v - np.nanmean(v)) / (np.nanstd(v) + 1e-9)
    return z(ph) - z(pc)


def prioritize(condition="Stim48hr", top_n=25, context=IMMUNE_CONTEXT,
               require_ontarget=True, model="claude-sonnet-5", email=None, verbose=True):
    """End-to-end: Marson hits -> PubMed grounding -> Claude verdict -> ranked shortlist.

    Returns a DataFrame sorted so strong + under-studied hits ("validate these") are on top.
    """
    hits = marson_hits(condition=condition, require_ontarget=require_ontarget)
    if verbose:
        print(f"{len(hits)} hits in {condition}; scoring top {top_n} by phenotype.")
    hits = add_pubmed_counts(hits, context=context, email=email, top_n=top_n, verbose=verbose)
    hits = score_hits(hits, context=context, model=model)
    hits["priority_score"] = _priority_score(hits)
    return hits.sort_values("priority_score", ascending=False).reset_index(drop=True)


def plot(df, title="Validate these: strong phenotype x under-studied", annotate=True):
    """The money figure: phenotype strength (y) vs literature coverage (x).

    Top-LEFT quadrant = strong phenotype + few papers = the validate-these hits (shaded).
    x = log10(PubMed papers + 1); y = phenotype strength. Points colored by Claude's
    recommendation when available, else by the grounded priority score. Returns a Figure.
    """
    import matplotlib.pyplot as plt

    d = df.dropna(subset=["pubmed_count"]).copy()
    x = np.log10(d["pubmed_count"].to_numpy(float) + 1)
    y = d["phenotype"].to_numpy(float)
    xm, ym = np.median(x), np.median(y)

    fig, ax = plt.subplots(figsize=(8, 6))
    # shade the "validate these" quadrant: low coverage (left), high phenotype (top)
    ax.axvspan(x.min() - 0.3, xm, ymin=0, ymax=1, color="#2ca02c", alpha=0.06)
    ax.axhline(ym, color="0.8", lw=1, zorder=0)
    ax.axvline(xm, color="0.8", lw=1, zorder=0)

    rec = d.get("recommendation")
    if rec is not None and rec.astype(str).str.len().gt(0).any():
        cmap = {"prioritize": "#2ca02c", "maybe": "#ff7f0e", "skip": "#999999"}
        colors = [cmap.get(r, "#1f77b4") for r in rec.astype(str)]
        for label, c in cmap.items():
            ax.scatter([], [], c=c, label=label)
        ax.legend(title="Claude verdict", loc="lower right", frameon=False)
    else:
        colors = d["priority_score"] if "priority_score" in d else "#1f77b4"

    sc = ax.scatter(x, y, c=colors, s=90, edgecolor="k", linewidth=0.5, zorder=3)
    if not isinstance(colors, list) and not isinstance(colors, str):
        fig.colorbar(sc, ax=ax, label="priority score")

    if annotate:
        in_quadrant = (x <= xm) & (y >= ym)
        for xi, yi, gi, q in zip(x, y, d["gene"], in_quadrant):
            ax.annotate(gi, (xi, yi), fontsize=9 if q else 7,
                        fontweight="bold" if q else "normal",
                        xytext=(4, 3), textcoords="offset points",
                        color="black" if q else "0.4")

    ax.set_xlabel("literature coverage  —  log10(PubMed papers + 1)")
    ax.set_ylabel("phenotype strength")
    ax.set_title(title)
    ax.text(0.02, 0.88, "← validate these\n   (strong + under-studied)",
            transform=ax.transAxes, va="top", ha="left", fontsize=10, color="#2ca02c",
            fontweight="bold")
    fig.tight_layout()
    return fig
