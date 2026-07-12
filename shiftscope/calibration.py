"""When does ShiftScope actually work? Calibrate its detection limits, honestly.

A distance-with-a-p-value is only trustworthy if you know where it breaks. This module maps
the two operating boundaries that matter for a real experiment, both on a shift we *know* is
real (Kang ctrl vs IFN-β monocytes, E-distance ≈ 220):

  1. **Cell number** — how few cells per group can you have and still detect the shift? We
     measure this on a *subtle* version of the shift (the query partly diluted, below),
     because for a huge effect the answer is trivially "almost none" — the interesting
     regime is a small effect where cell count is the limiting factor.
  2. **Effect size** — how small a shift is still detectable? Dilute the query group with a
     fraction f of *reference* cells; the true difference shrinks toward zero as f→1. The
     largest f still detected is the minimum detectable shift.

Both boil down to the same move: draw n reference cells and n (optionally diluted) query
cells, run the E-test, repeat over many random draws, and report **power** — the fraction of
draws reaching p < 0.05 — not one lucky call. Output: a table per axis, a 2×2 figure, and one
plain-English sentence. Same E-distance / E-test machinery as `compare.py`, stress-tested.
"""

import numpy as np
import pandas as pd

from . import compare


def _draw(Z, k, rng):
    """Draw k rows from Z (without replacement when possible, else with)."""
    return Z[rng.choice(len(Z), k, replace=k > len(Z))]


def _power_point(Z_ref, Z_query, n, f, reps, n_perm, alpha, rng):
    """Mean E-distance and E-test power for one (cells=n, dilution=f) setting.

    Each rep draws n reference cells and an n-cell query group in which a fraction f has been
    replaced by reference cells (f=0 → the true query; f→1 → no shift). Returns
    (edist_mean, edist_sd, power) where power is the fraction of reps with p < alpha.
    """
    k_ref = int(round(f * n))      # reference cells mixed into the pseudo-query
    k_qry = n - k_ref              # remaining real query cells
    eds, ps = [], []
    for _ in range(reps):
        ref = _draw(Z_ref, n, rng)
        query = (np.vstack([_draw(Z_query, k_qry, rng), _draw(Z_ref, k_ref, rng)])
                 if k_qry > 0 else _draw(Z_ref, n, rng))
        e, p = compare.etest(ref, query, n_perm=n_perm, seed=int(rng.integers(1 << 31)))
        eds.append(e)
        ps.append(p)
    return float(np.mean(eds)), float(np.std(eds)), float(np.mean(np.array(ps) < alpha))


def sensitivity_vs_cells(Z_ref, Z_query, ns=(10, 25, 50, 100, 200, 400), dilute=0.85,
                         reps=20, n_perm=200, alpha=0.05, seed=0):
    """E-distance and E-test power vs cells per group, at a fixed (subtle) effect size.

    `dilute` fixes the shift's magnitude (fraction of the query replaced by reference cells)
    so the sweep probes the cell-count limit rather than the trivially-detectable full shift.
    Returns a tidy DataFrame [n_per_group, edist_mean, edist_sd, power].
    """
    rng = np.random.default_rng(seed)
    rows = []
    for n in ns:
        m, sd, pw = _power_point(Z_ref, Z_query, n, dilute, reps, n_perm, alpha, rng)
        rows.append({"n_per_group": n, "edist_mean": m, "edist_sd": sd, "power": pw})
    return pd.DataFrame(rows)


def sensitivity_vs_effect(Z_ref, Z_query, fracs=(0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95),
                          n=200, reps=20, n_perm=200, alpha=0.05, seed=0):
    """E-distance and E-test power vs how diluted the shift is, at fixed cell count `n`.

    The query is rebuilt as `(1-f)` real query cells + `f` reference cells, so at f=0 it is
    the true query and at f→1 the shift vanishes. Returns a tidy DataFrame; the smallest
    detectable shift is the largest `f` whose power is still high.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for f in fracs:
        m, sd, pw = _power_point(Z_ref, Z_query, n, f, reps, n_perm, alpha, rng)
        rows.append({"dilution_frac": f, "edist_mean": m, "edist_sd": sd, "power": pw})
    return pd.DataFrame(rows)


def _summary(res_cells, res_effect, dilute, n_effect, power_thresh=0.8):
    """Turn the two tables into two honest sentences about the operating range."""
    det = res_cells[res_cells["power"] >= power_thresh]
    cells_txt = (f"~{int(det['n_per_group'].min())} cells/group" if len(det)
                 else f">{int(res_cells['n_per_group'].max())} cells/group")
    kept = res_effect[res_effect["power"] >= power_thresh]
    max_dil = float(kept["dilution_frac"].max()) if len(kept) else 0.0
    return (
        f"Cell number: on a subtle shift (query {dilute:.0%}-diluted), reliable detection "
        f"(power ≥ {power_thresh:.0%}) needs {cells_txt}. "
        f"Effect size: at n={n_effect} cells/group the shift is still detected when up to "
        f"~{max_dil:.0%} of the query is reference cells — the minimum detectable shift."
    )


def calibrate(adata=None, ns=(10, 25, 50, 100, 200, 400),
              fracs=(0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95), dilute_for_cells=0.85,
              n_effect=200, reps=20, n_perm=200, n_pcs=30, seed=0, verbose=True):
    """Run both calibrations on a known shift and return tables + a plain-English verdict.

    With no `adata`, loads the default known-strong comparison (Kang ctrl vs IFN-β
    monocytes) and embeds it. Pass your own preprocessed, `_group`-tagged AnnData to
    calibrate a different shift. `dilute_for_cells` sets how subtle the shift is for the
    cell-number sweep. Returns `{"cells": df, "effect": df, "summary": str}`.
    """
    from . import io, embed

    if adata is None:
        adata = io.load_public("kang", "ctrl", "stim",
                               within=("cell_type", "CD14+ Monocytes"))
    Z_ref, Z_query = embed.embed(adata, n_pcs=n_pcs, run_umap=False, seed=seed)
    if verbose:
        print(f"calibrating on {len(Z_ref)} ref + {len(Z_query)} query cells "
              f"in {Z_ref.shape[1]} PCs")

    res_cells = sensitivity_vs_cells(Z_ref, Z_query, ns=ns, dilute=dilute_for_cells,
                                     reps=reps, n_perm=n_perm, seed=seed)
    res_effect = sensitivity_vs_effect(Z_ref, Z_query, fracs=fracs, n=n_effect,
                                       reps=reps, n_perm=n_perm, seed=seed)
    summary = _summary(res_cells, res_effect, dilute_for_cells, n_effect)
    if verbose:
        print(summary)
    return {"cells": res_cells, "effect": res_effect, "summary": summary}


def plot(res_cells, res_effect, alpha=0.05, dilute=0.85):
    """2×2 boundaries figure: E-distance (top) and detection power (bottom) vs each axis.

    Left column = vs cells/group (at a fixed subtle shift); right column = vs dilution
    fraction (at fixed cell count). Dashed line = the α significance floor; shaded band =
    power ≥ 0.8 (reliable detection). Returns a matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    c, e = res_cells, res_effect

    ax[0, 0].errorbar(c["n_per_group"], c["edist_mean"], yerr=c["edist_sd"],
                      marker="o", capsize=3, color="#1f77b4")
    ax[0, 0].set(xlabel="cells per group", ylabel="E-distance",
                 title=f"Effect vs cell number (shift {dilute:.0%}-diluted)")
    ax[0, 0].set_xscale("log")
    ax[1, 0].axhspan(0.8, 1.0, color="#2ca02c", alpha=0.08)
    ax[1, 0].plot(c["n_per_group"], c["power"], marker="o", color="#1f77b4")
    ax[1, 0].axhline(0.8, color="0.6", ls="--", lw=1)
    ax[1, 0].set(xlabel="cells per group", ylabel=f"detection power (p<{alpha:g})",
                 title="Power vs cell number", ylim=(-0.05, 1.05))
    ax[1, 0].set_xscale("log")

    ax[0, 1].errorbar(e["dilution_frac"], e["edist_mean"], yerr=e["edist_sd"],
                      marker="o", capsize=3, color="#d62728")
    ax[0, 1].set(xlabel="query diluted with reference (fraction)", ylabel="E-distance",
                 title="Effect vs dilution")
    ax[1, 1].axhspan(0.8, 1.0, color="#2ca02c", alpha=0.08)
    ax[1, 1].plot(e["dilution_frac"], e["power"], marker="o", color="#d62728")
    ax[1, 1].axhline(0.8, color="0.6", ls="--", lw=1)
    ax[1, 1].set(xlabel="query diluted with reference (fraction)",
                 ylabel=f"detection power (p<{alpha:g})",
                 title="Power vs dilution", ylim=(-0.05, 1.05))

    fig.suptitle("ShiftScope operating boundaries (known IFN-β shift)", fontweight="bold")
    fig.tight_layout()
    return fig
