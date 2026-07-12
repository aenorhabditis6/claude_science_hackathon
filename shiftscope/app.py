"""Gradio UI for ShiftScope.

Pick a dataset and a (reference, query) pair -> distances table, UMAP colored by group
and by the shifted clusters, a driver-gene panel, and Claude's cited write-up. Runs the
whole pipeline (io -> embed -> compare -> localize -> drivers -> interpret) on one click.
Launches with share=True so it works from Colab.

The public datasets (kang / pbmc68k / hagai) run in ~20 s and are the reliable live demo.
"Marson (genome-scale)" streams from the 140 GB remote Perturb-seq file — impressive but
slow unless run from a low-latency network (Colab).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gradio as gr
import scanpy as sc

from shiftscope import io, embed, compare, localize, drivers, interpret

# name -> (reference, query, within "col=value" or ""). Label starts with the dataset key.
# Spans domains so the app itself tells the "any two populations" story: immune, disease,
# development, T-cell CRISPR, cell types.
PRESETS = {
    "kang (immune: ctrl vs IFN-beta)": ("ctrl", "stim", "cell_type=CD14+ Monocytes"),
    "haber (disease: gut, healthy vs infection)": ("Control", "Hpoly.Day10", ""),
    "paul15 (development: blood cell fate)": ("14Mo", "2Ery", ""),
    "shifrut (T-cell CRISPR: ctrl vs LCP2 KO)": ("control", "LCP2", "perturbation_2=stim"),
    "datlinger (T cells: unstim vs stim)": ("unstimulated", "stimulated", ""),
    "hagai (phagocytes: unst vs LPS)": ("unst", "LPS6", "species=mouse"),
    "pbmc68k (compare two cell types)": ("CD14+ Monocyte", "CD19+ B", ""),
    "pbmc3k (classic PBMC cell types)": ("CD4 T cells", "CD14+ Monocytes", ""),
    "Marson (genome-scale, slow)": ("NTC", "SETDB1", "condition=Rest"),
}


def _umap_figure(adata):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc.pl.umap(adata, color="_group", ax=axes[0], show=False, frameon=False,
               title=f"{adata.uns['_reference']} vs {adata.uns['_query']}")
    color = "_shift" if "_shift" in adata.obs else "_group"
    sc.pl.umap(adata, color=color, ax=axes[1], show=False, frameon=False,
               title="Shifted clusters" if color == "_shift" else "group")
    fig.tight_layout()
    return fig


def _load(dataset, reference, query, within):
    within_tuple = None
    if within.strip():
        col, val = within.split("=", 1)
        within_tuple = (col.strip(), val.strip())
    if dataset.startswith("Marson"):
        condition = within_tuple[1] if within_tuple else "Rest"
        return io.load_demo(query=query, condition=condition)
    return io.load_public(dataset.split()[0], reference, query, within=within_tuple)


def run(dataset, reference, query, within, do_interpret):
    """Full pipeline for one comparison. Returns the UI outputs."""
    adata = _load(dataset, reference, query, within)
    Z_ref, Z_query = embed.embed(adata)
    dist = compare.distances(Z_ref, Z_query)

    dist_md = (
        f"### {adata.uns['_reference']} vs {adata.uns['_query']}  ·  {dataset}\n\n"
        f"| metric | value |\n|---|---|\n"
        f"| **E-distance** (scPerturb) | **{dist['edistance']:.2f}** |\n"
        f"| **E-test p-value** | **{dist['p_value']:.4f}** ({dist['n_perm']} permutations) |\n"
        f"| MMD (RBF) | {dist['mmd']:.4f} |\n"
        f"| Sinkhorn-OT | {dist['sinkhorn']:.4f} |\n"
        f"| cells (ref / query) | {dist['n_ref']} / {dist['n_query']} |\n"
    )

    shift_tab = localize.localize(adata)
    sig = localize.shifted_clusters(shift_tab)
    drv = drivers.drivers(adata, clusters=sig or None)

    fig = _umap_figure(adata)
    drivers_md = (
        "### Top driver genes (query vs reference"
        + (f", in shifted clusters {sig})" if sig else ")") + "\n\n"
        + "**Up:** " + ", ".join(drv["up"]["gene"].head(15)) + "\n\n"
        + "**Down:** " + ", ".join(drv["down"]["gene"].head(15))
    )

    interp_md = "_Interpretation off._"
    if do_interpret:
        if not os.environ.get("ANTHROPIC_API_KEY"):
            interp_md = "_Set `ANTHROPIC_API_KEY` to enable Claude's interpretation._"
        else:
            ev = interpret.assemble_evidence(adata, dist, shift_tab, drv)
            interp_md = "### Claude's read on the biology\n\n" + interpret.interpret(ev)["output"]

    return dist_md, fig, drivers_md, interp_md


def run_prioritize(condition, top_n):
    """Prioritize tab: rank the real Marson screen hits by strong-phenotype x under-studied."""
    from shiftscope import prioritize

    df = prioritize.prioritize(condition=condition, top_n=int(top_n), verbose=False)
    fig = prioritize.plot(df)
    cols = ["gene", "phenotype", "pubmed_count", "priority_score",
            "novelty", "recommendation", "verdict"]
    return df[cols], fig


def build_app():
    with gr.Blocks(title="ShiftScope") as demo:
        gr.Markdown(
            "# ShiftScope\n"
            "Compare **any two single-cell populations**, measure how the cell-state "
            "distribution shifts, localize where, and let Claude interpret it. "
            "Condition-agnostic: control vs perturbation, disease vs healthy, young vs old."
        )
        with gr.Tab("Compare two populations"):
            with gr.Row():
                dataset = gr.Dropdown(list(PRESETS), value=list(PRESETS)[0], label="Dataset")
                reference = gr.Textbox("ctrl", label="Reference")
                query = gr.Textbox("stim", label="Query")
                within = gr.Textbox("cell_type=CD14+ Monocytes",
                                    label="Restrict to (col=value, optional)")
                do_interpret = gr.Checkbox(True, label="Claude interpretation")
            run_btn = gr.Button("Compare", variant="primary")

            def _fill(name):
                r, q, w = PRESETS[name]
                return r, q, w
            dataset.change(_fill, dataset, [reference, query, within])

            dist_out = gr.Markdown()
            umap_out = gr.Plot()
            drivers_out = gr.Markdown()
            interp_out = gr.Markdown()
            run_btn.click(run, [dataset, reference, query, within, do_interpret],
                          [dist_out, umap_out, drivers_out, interp_out])

        with gr.Tab("Prioritize screen hits (Marson)"):
            gr.Markdown(
                "Rank real Marson CD4+ T-cell CRISPR hits by **strong phenotype x "
                "under-studied**. Novelty is grounded in a live PubMed count; Claude scores "
                "each hit and gives a validate/skip verdict. Top-left of the figure = validate."
            )
            with gr.Row():
                condition = gr.Dropdown(["Stim48hr", "Stim8hr", "Rest"], value="Stim48hr",
                                        label="Culture condition")
                top_n = gr.Slider(5, 40, value=20, step=5, label="Top hits to score")
            prio_btn = gr.Button("Rank + prioritize", variant="primary")
            prio_table = gr.Dataframe(label="Shortlist (validate-these on top)", wrap=True)
            prio_plot = gr.Plot()
            prio_btn.click(run_prioritize, [condition, top_n], [prio_table, prio_plot])
    return demo


if __name__ == "__main__":
    build_app().launch(share=True)
