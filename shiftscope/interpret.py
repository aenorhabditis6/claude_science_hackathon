"""Ask Claude to turn the numbers into a cited biological rationale.

We hand Claude only the structured evidence ShiftScope computed — distances, the shifted
clusters and their direction, and the top driver genes — and ask for a short write-up that
walks cell states -> pathway -> candidate target. Every claim must cite a number we passed
in, so we log the exact prompt + evidence alongside the output and every claim traces back
to data.

Uses the anthropic SDK. Default model is Claude Sonnet 5 (fast, and the project's
`CLAUDE.md` says Sonnet is fine here); pass `model="claude-opus-4-8"` for max quality.
Needs `ANTHROPIC_API_KEY` in the environment (or an `ant auth login` profile).
"""

import json
import os
import time

MODEL = "claude-sonnet-5"

SYSTEM = """You are a careful computational biologist interpreting a single-cell \
comparison. You are given only structured evidence computed by an upstream tool \
(distribution distances, which cell clusters shifted, and the top differentially \
expressed genes). Write a short, decision-useful rationale for a bench scientist.

Rules:
- Ground every claim in a number from the evidence. When you cite one, name it \
(e.g. "E-distance 14.9, p=0.001" or "cluster 3, log2 enrichment +1.8").
- Do not invent genes, pathways, cell types, statistics, or citations. If the evidence \
is insufficient for a claim, say so.
- Structure: (1) is the shift real and how large; (2) where it localizes in cell-state \
space; (3) what the driver genes suggest about pathway / cell state; (4) one candidate \
mechanism or target to follow up, with the caveat that this is hypothesis-generating.
- Be concise. Flag uncertainty honestly."""


def assemble_evidence(adata, distances, shift_table=None, drivers=None):
    """Collect the pipeline outputs into one JSON-able evidence dict for Claude."""
    ev = {
        "dataset": adata.uns.get("_dataset", "unknown"),
        "reference": adata.uns["_reference"],
        "query": adata.uns["_query"],
        "n_reference_cells": int((adata.obs["_group"] == adata.uns["_reference"]).sum()),
        "n_query_cells": int((adata.obs["_group"] == adata.uns["_query"]).sum()),
        "distances": {k: distances[k] for k in ("edistance", "p_value", "mmd", "sinkhorn")
                      if k in distances},
    }
    if shift_table is not None:
        sig = shift_table[shift_table["direction"] != "ns"]
        ev["shifted_clusters"] = [
            {
                "cluster": str(r.cluster),
                "direction": r.direction,
                "log2_enrichment": round(float(r.log2_enrichment), 2),
                "qval": float(f"{r.qval:.2e}"),
                "n_query": int(r.n_query),
                "n_ref": int(r.n_ref),
            }
            for r in sig.itertuples()
        ]
    if drivers is not None:
        ev["top_up_genes"] = drivers["up"]["gene"].head(15).tolist()
        ev["top_down_genes"] = drivers["down"]["gene"].head(15).tolist()
    return ev


def interpret(evidence, model=MODEL, log_dir="interpretations", max_tokens=1500):
    """Send the evidence to Claude; return the write-up and log the exact prompt+inputs."""
    import anthropic

    user = (
        "Here is the structured evidence from a ShiftScope comparison. Interpret it per "
        "your instructions.\n\n```json\n" + json.dumps(evidence, indent=2) + "\n```"
    )

    client = anthropic.Anthropic()
    resp = client.messages.create(
        model=model, max_tokens=max_tokens, system=SYSTEM,
        messages=[{"role": "user", "content": user}],
    )
    if resp.stop_reason == "refusal":
        text = "[Claude declined to answer this request.]"
    else:
        text = "".join(b.text for b in resp.content if b.type == "text")

    record = {
        "model": model,
        "system": SYSTEM,
        "user_prompt": user,
        "evidence": evidence,
        "output": text,
        "stop_reason": resp.stop_reason,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
        path = os.path.join(log_dir, f"interpret_{int(time.time())}.json")
        with open(path, "w") as f:
            json.dump(record, f, indent=2)
        record["log_path"] = path
    return record
