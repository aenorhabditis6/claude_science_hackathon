# ShiftScope — 5-minute pitch deck plan

**Event:** Built with Claude: Life Sciences (Anthropic × Gladstone). **Format:** ~5 min live.
**Judging:** Demo 30% · Claude-use 25% · Impact 25% · Depth 20%. Every slide below is tagged
with the criterion it feeds. Rule of thumb: **talk ≤ 40s/slide, one idea per slide, one
figure per slide.** Lead with impact, prove with a live demo, win on Claude + rigor.

## The arc (why this order)
Hook → what it is → **it works everywhere** (impact) → **watch it work** (live demo) →
**the differentiator: Claude picks what to validate** (Claude-use) → **we know its limits**
(depth) → **it scales** → close. The two "wow" beats are Slide 3 (four fields, one tool) and
Slide 5 (Claude turns a metric into a decision).

## Assets ready in `figures/`
- `fig_domains.png` — E-distance across 4 fields + driver genes (**Slide 3**, the Impact money shot)
- `fig_kang_umap.png` — ctrl/stim UMAP + localized shift (**Slide 4**, live-demo backup)
- `fig_prioritize.png` — the "validate these" 2×2, Claude verdicts (**Slide 5**, the differentiator)
- `fig_calibration.png` — power vs cells & effect size (**Slide 6**, rigor)

---

## Slide-by-slide

### 1 — Title & hook  ·  *Impact*  ·  ~20s
- **ShiftScope — compare any two single-cell populations, and know what to do next.**
- One line: "Point it at two groups of cells — disease vs healthy, drug vs control, young vs
  old — and it tells you *how far apart* they are, *where* they diverge, *which genes* drive
  it, and *what to validate*."
- Visual: title + a faint UMAP with two colored clouds and an arrow between them.
- Note: say the tag-line verbatim; it's the whole pitch.

### 2 — The problem  ·  *Impact / Depth*  ·  ~35s
- Single-cell comparison is the most common question in the field, done ad hoc every time:
  eyeball a UMAP, run DE, argue about it. Two real gaps: **(a) is the shift real and how big?**
  and **(b) of hundreds of hits, which do you spend a month validating?**
- Visual: a messy "before" — a UMAP with a hand-drawn "…so what?" and a 500-row DE table.
- Note: name the pain the judges feel in their own labs.

### 3 — One tool, four fields  ·  *Impact (money slide)*  ·  ~45s
- **`fig_domains.png`.** Same five lines of code across immunology, infection/disease,
  development, and cell-type ID — correct, distinct biology every time (interferon → NF-κB →
  antimicrobial Reg3g → erythroid Klf1). E-distance is *graded* and meaningful.
- Note: "Nothing here is hardcoded to any disease. Swap the column, get a new comparison."
  This is the condition-agnostic claim, proven on 4 real public datasets.

### 4 — Live demo: measure → localize → drivers  ·  *Demo (30%)*  ·  ~60s
- Run the notebook live on Kang ctrl vs IFN-β monocytes (~20s): print E-distance + E-test p,
  show the UMAP colored by group and by *shifted cluster*, and the driver genes (ISG15/MX1/…).
- Visual: `fig_kang_umap.png` as backup if live run is slow.
- Note: this is the core loop working end-to-end on real data in front of them. Keep talking
  while it runs; have the figure ready if wifi/compute stalls.

### 5 — The differentiator: Claude decides what to validate  ·  *Claude-use (25%) + Impact*  ·  ~60s
- **`fig_prioritize.png`.** On the real Marson genome-scale T-cell CRISPR screen, ShiftScope
  ranks hits by **strong phenotype × under-studied**. The novelty axis is **grounded** — a live
  PubMed count, not an LLM guess — and Claude scores each hit + gives a validate/skip verdict
  *with that count in front of it*. Watch it **sink the textbook genes** (VAV1, CD45 — hundreds
  of papers) and **surface strong-but-obscure SAGA chromatin regulators** (TADA2B, SGF29, ELOF1
  — 0–4 papers) as the shortlist.
- Note: THIS is the Claude-use win — agentic, grounded, cited, logged. "We don't just measure
  the shift; Claude turns it into a decision, and every call traces back to a real number."

### 6 — We know where it breaks  ·  *Depth (20%)*  ·  ~35s
- **`fig_calibration.png`.** Honest operating envelope: power collapses below ~100 cells/group,
  and the shift stays detectable until it's diluted to <~10% of its true size. Reported as
  *power over many random draws*, not one lucky p-value.
- Note: "A metric you can't trust is worthless — so we mapped its limits." Reviewer catnip.

### 7 — And it scales  ·  *Depth / Impact*  ·  ~30s
- Streams the 140 GB Marson Perturb-seq screen from S3 without downloading it; the same
  `rank()` recovers a *second* real CRISPR screen (Shifrut human T cells → LCP2/CD3D + the
  RASA2/CBLB/SOCS1 negative regulators = immunotherapy targets).
- Visual: the Shifrut rank table + a one-line "140 GB → streamed, never downloaded."
- Note: depth + on-theme for the Gladstone/Marson immunology audience.

### 8 — Close  ·  *Impact / Claude*  ·  ~20s
- "ShiftScope: measure the shift, localize it, find the drivers, and let Claude tell you what
  to validate — for **any** two cell populations." MIT-licensed, runs in Colab.
- Recap the Claude story in one line: **grounded, agentic, and fully auditable.**
- Team + repo link. End on the tag-line.

---

## Live-demo checklist (de-risk Slide 4)
- Pre-run the notebook once so datasets are cached; keep the browser tab open at the demo cell.
- Have `fig_kang_umap.png` + `fig_prioritize.png` on the next slides as instant fallbacks.
- `ANTHROPIC_API_KEY` set in the Colab secret before you start (Slide 5 verdicts need it).
- If offline: the whole story runs from the four figures alone.

## Design notes
- 16:9, big type (≥24pt body), one figure per slide, minimal text — you narrate.
- Domain color code is consistent (blue=immunology, red=disease, green=development,
  purple=cell types) — reuse it across slides 3–7 so the eye tracks the "many fields" story.
- Keep code off the slides except one 5-line snippet on Slide 3 to show how little it takes.
- Two "wow" beats (S3, S5) — pause a beat after each.

## Build options (pick one, I can generate it)
1. **Google Slides / Keynote** — you build from this outline + the 4 PNGs (most control).
2. **Auto-built `.pptx`** — I compile this outline + figures into an editable deck via the
   pptx skill (fastest to a first draft; you restyle).
3. **The notebook *is* the deck** — present `demo.ipynb` top-to-bottom; slides only for the
   hook, the domains figure, and the close. (Least to maintain; most "live".)

---

# Slide contents — paste-ready

Copy each block onto a slide. **On-slide** = what the audience reads (keep it terse).
**Say** = your spoken line. Numbers are the validated demo values.

### Slide 1 — Title
**On-slide:**
- # ShiftScope
- Compare any two single-cell populations — and know what to validate next.
- *Built with Claude · Life Sciences (Anthropic × Gladstone)*
**Say:** "ShiftScope takes two groups of cells — disease vs healthy, drug vs control, any two —
and tells you how far apart they are, where they diverge, which genes drive it, and what to
validate next."

### Slide 2 — The problem
**On-slide:**
- ## The question every single-cell study asks: *"how are these two populations different?"*
- Answered ad hoc every time — eyeball a UMAP, run DE, argue.
- Two gaps stay hard: **(1)** is the shift real, and how big? **(2)** of hundreds of hits,
  which do you spend a month validating?
**Say:** "Everyone eyeballs a UMAP and runs differential expression. But two questions stay
hard: is the difference real and sizeable — and which hit actually deserves bench time?"

### Slide 3 — One tool, every domain  *(Impact)*
**On-slide:** *(figure: `fig_domains.png`)*
- The **same five lines of code** across immunology, disease, development, and cell-type ID
- Correct, *distinct* biology every time; E-distance is graded and meaningful
- Nothing hardcoded — swap the `.obs` column, get a new comparison
**Say:** "Same code, four different fields — interferon, an antimicrobial infection response, an
erythroid program, myeloid identity. It's condition-agnostic by design."

### Slide 4 — Watch it work (live)  *(Demo)*
**On-slide:** *(live notebook; `fig_kang_umap.png` as backup)*
- ctrl vs IFN-β monocytes: **E-distance 216, p = 0.001** (E-test)
- Shift **localizes** to specific clusters (UMAP)
- Drivers: **ISG15, MX1, IFIT1/3, OAS1** — textbook interferon
- Claude writes the cited rationale
**Say:** "Twenty seconds, real 10x data, end to end — measure, localize, drivers, and Claude's
read." *(Run the notebook; keep talking; use the figure if it stalls.)*

### Slide 5 — From a metric to a decision  *(Claude — the differentiator)*
**On-slide:** *(figure: `fig_prioritize.png`)*
- Rank hits by **strong phenotype × under-studied**
- Novelty is **grounded in a live PubMed count** — not an LLM guess
- Claude scores + gives a validate/skip verdict *with that count in front of it*; every call logged
- **Sinks** textbook genes (VAV1, CD45); **surfaces** obscure SAGA regulators (TADA2B, SGF29, ELOF1)
**Say:** "This is where Claude earns its keep — agentic, grounded, and auditable. It turns the
shift into a shortlist a scientist can act on, and every claim traces back to a real number."

### Slide 6 — Honest about where it breaks  *(Depth)*
**On-slide:** *(figure: `fig_calibration.png`)*
- Reports **power over many random draws**, not one lucky p-value
- Needs **~100 cells/group** on a subtle shift; detectable until **~90% dilution**
**Say:** "A metric you can't trust is worthless — so we mapped the operating envelope."

### Slide 7 — From toy to genome-scale  *(Depth / Impact)*
**On-slide:** *(shifrut rank table)*
- Streams the **140 GB Marson Perturb-seq** from S3 — never downloaded
- Same `rank()` recovers a **second real CRISPR screen** (Shifrut human T cells):
  **LCP2/CD3D** (TCR core) + **RASA2/CBLB/SOCS1** (the regulators that *enhance* T cells) —
  immunotherapy targets
**Say:** "It runs on a 38 MB demo and on a 140 GB genome-scale screen, streamed — same code."

### Slide 8 — Close
**On-slide:**
- ## ShiftScope
- Measure the shift · localize it · find the drivers · **let Claude tell you what to validate**
- For **any** two cell populations. MIT-licensed, runs in Colab.
- *Claude: grounded, agentic, auditable.*
- github.com/aenorhabditis6/claude_science_hackathon
**Say:** "ShiftScope — for any two cell populations, measure the shift and let Claude tell you
what to validate."
