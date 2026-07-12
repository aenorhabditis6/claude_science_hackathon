# ShiftScope — submission materials

Built with Claude: Life Sciences · **Builder Track** · MIT-licensed · runs in Colab.
Deadline: Mon Jul 13, 9:00 PM ET. Required uploads: 3-min video · repo/notebook · 100–200-word summary.

---

## Written summary (paste into the form) — ~175 words

ShiftScope compares **any two single-cell populations** — disease vs healthy, drug vs
control, any two — and turns the comparison into a decision. It measures how far the
cell-state distribution shifts (the scPerturb E-distance with a permutation E-test, plus MMD
and optimal transport), localizes *where* in cell-state space the shift concentrates, finds
the *driver genes*, and has Claude write a cited, auditable rationale.

Its differentiator: on Alex Marson's genome-scale CD4+ T-cell CRISPR Perturb-seq screen,
ShiftScope ranks hits by **strong-phenotype × under-studied** — grounding novelty in *live
PubMed counts*, then asking Claude for a validate/skip verdict with that evidence in hand. It
sinks textbook genes (VAV1, CD45) and surfaces strong-but-obscure candidates (the SAGA
chromatin complex) worth the bench.

The same five lines run across immunology, host–pathogen disease, development, and cell-type
identity — recovering the correct, distinct biology each time — and we map the method's
detection limits honestly. ShiftScope turns "these populations look different" into "validate
*this* gene next."

---

## 3-minute demo video — recordable script

**Format:** screen-record the **already-executed notebook** (rendered outputs + the four
figures) with voiceover — safer than live-running on camera. Loom or QuickTime is fine. Keep
it to 3:00. Show the figures full-screen when referenced. Optional: 10 s of the Gradio app.

| Time | On screen | Say (voiceover) |
|---|---|---|
| **0:00–0:20** | Title slide / notebook top | "This is ShiftScope. Point it at two groups of cells — disease vs healthy, drug vs control, any two — and it tells you how far apart they are, where they diverge, which genes drive it, and **what to validate next**." |
| **0:20–0:50** | `fig_domains.png` (the 4-field bar chart) | "The same five lines of code run across four different fields — immunology, host–pathogen disease, development, cell-type ID — and recover the correct, distinct biology every time: interferon, antimicrobials, an erythroid program. Nothing is hardcoded to any disease." |
| **0:50–1:35** | Scroll the Kang cells: distances → UMAP (`fig_kang_umap`) → drivers → Claude write-up | "Here it is end-to-end on real data: a large, significant shift — E-distance 216, p=0.001 — that localizes to specific clusters, driven by textbook interferon genes. Then Claude writes a rationale where **every claim cites a number we computed**, and we log the exact prompt." |
| **1:35–2:25** | `fig_prioritize.png` + the shortlist table | "The differentiator: on Marson's genome-scale T-cell CRISPR screen, ShiftScope ranks hits by **strong phenotype × how under-studied they are**. Novelty is grounded in a *live PubMed count* — not a guess — and Claude gives a validate/skip verdict with that count in front of it. Watch it bench the famous genes and surface strong, obscure chromatin regulators as the shortlist. That's the decision a scientist actually faces." |
| **2:25–2:50** | `fig_calibration.png`, then the Shifrut rank table | "We're honest about limits — here's exactly where detection breaks by cell number and effect size. And it scales: the same call recovers a second real CRISPR screen, and streams Marson's 140 GB dataset without downloading it." |
| **2:50–3:00** | Title / repo URL | "ShiftScope — for any two cell populations, measure the shift and let Claude tell you what to validate. Open source, runs in Colab." |

**Recording tips:** rehearse once for timing; speak to the *decision* (what to validate), not
the plumbing; pause a beat on the domains figure and the prioritize figure — those are the two
"wow" moments. If you show the app, do it live for ~10 s on the Prioritize tab.

---

## Submission checklist
- [ ] Merge PR #4 → `main` has the rendered notebook + no-restart setup
- [ ] Record the 3-min video (script above) → upload to YouTube/Loom, get the link
- [ ] Paste the written summary into the CV form
- [ ] Repo link: https://github.com/aenorhabditis6/claude_science_hackathon (public ✅, MIT ✅)
- [ ] Submit at cerebralvalley.ai before **Mon Jul 13, 9:00 PM ET**
