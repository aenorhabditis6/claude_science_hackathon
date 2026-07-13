"""One call, one look: run the whole pipeline on a demo dataset and show the result.

For anyone who just wants to see ShiftScope work. `run()` loads a ready-made public dataset,
measures the shift, localizes it, finds the driver genes, and (if a Claude key is set) writes
the interpretation. Then it prints a short summary and shows the UMAP. Everything it calls is
the same code the full pipeline uses; this is only a friendly front door.

    from shiftscope.quickstart import run
    run()
"""

import os


def run(name="kang", reference="ctrl", query="stim",
        within=("cell_type", "CD14+ Monocytes"), interpret=True, show=True):
    """Run ShiftScope end to end on one comparison and show the result. Returns the AnnData.

    Defaults to control vs IFN-β-stimulated monocytes (Kang 2018), about 20 s on CPU. Swap
    `name` and the (reference, query) pair for anything in `io.PUBLIC_DATASETS`, e.g.
    `run("haber", "Control", "Hpoly.Day10", within=None)` for gut infection.
    """
    from . import io, embed, compare, localize, drivers, interpret as interp

    adata = io.load_public(name, reference, query, within=within)
    Z_ref, Z_query = embed.embed(adata)
    dist = compare.distances(Z_ref, Z_query)
    shift = localize.localize(adata)
    drv = drivers.drivers(adata)

    ref, qry = adata.uns["_reference"], adata.uns["_query"]
    print(f"\n{ref}  vs  {qry}   ({adata.n_obs} cells)")
    print(f"  E-distance {dist['edistance']:.1f}   E-test p = {dist['p_value']:.3f}"
          f"   (MMD {dist['mmd']:.3f}, Sinkhorn {dist['sinkhorn']:.3f})")
    print(f"  up in {qry}:   {', '.join(drv['up']['gene'].head(10))}")
    print(f"  down in {qry}: {', '.join(drv['down']['gene'].head(8))}")

    if show:
        import scanpy as sc
        sc.pl.umap(adata, color=["_group", "_shift"], frameon=False,
                   title=[f"{ref} vs {qry}", "shifted clusters"])

    if interpret and os.environ.get("ANTHROPIC_API_KEY"):
        ev = interp.assemble_evidence(adata, dist, shift, drv)
        print("\nClaude's read:\n" + interp.interpret(ev)["output"])
    elif interpret:
        print("\n(set ANTHROPIC_API_KEY to add Claude's interpretation)")

    return adata
