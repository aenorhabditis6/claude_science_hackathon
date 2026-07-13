# ShiftScope submission materials

Built with Claude: Life Sciences. Builder Track. MIT licensed. Runs in Colab.
Deadline: Mon Jul 13, 9:00 PM ET. Uploads: 3-minute video, repo/notebook, 100 to 200 word summary.

---

## Written summary (paste into the form)

ShiftScope compares two groups of single cells and helps you decide what to do about the
difference. You give it one dataset and name two conditions. That can be disease versus
healthy, drug versus control, or two cell types. It measures how far the two populations have
moved apart, using the E-distance and permutation test from scPerturb, plus MMD and optimal
transport. It finds which clusters gained or lost cells, and which genes drive the change.
Then it asks Claude to explain the biology, and every sentence Claude writes points back to a
number we gave it.

The most useful part runs on Alex Marson's CD4+ T-cell CRISPR screen. ShiftScope ranks the
hits by how strong the effect is and how little the gene has been studied. It counts how many
papers mention each gene on PubMed, then asks Claude for a keep-or-skip call using that real
number. Well-known genes like VAV1 and CD45 drop down the list. Strong but barely-studied
genes rise to the top as candidates worth testing at the bench.

The same code works across immune signaling, infection, blood development, and cell types. It
runs in Colab and is MIT licensed.

---

## Video plan (no editing needed)

You are not editing video. Make the slides, then let PowerPoint or Keynote record and export
the video for you.

**How to record with zero editing:**
1. Open the deck (`shiftscope_deck.pptx` in `figures/`, or your own).
2. PowerPoint: **Slide Show > Record**. Keynote: **Play > Record Slideshow**. Talk over each
   slide. It saves your voice and timing per slide.
3. Export to a movie file. PowerPoint: **File > Export > Create a Video** (pick 1080p).
   Keynote: **File > Export To > Movie**. That gives you the MP4 to upload. No editing tools.
4. Keep it under 3 minutes. Practice once for timing.

**Want to show it actually running (optional, looks good):**
Record your screen for about 20 seconds with QuickTime (File > New Screen Recording) or Loom,
showing the Gradio app: pick a dataset, click Compare, results appear. Drag that clip onto one
slide. PowerPoint plays it while you record the slideshow. Do not run code live on camera. If
you want to show a run, record it first, then drop in the clip.

**The script is in each slide's Speaker Notes**, so read it while you record. Seven slides,
about 25 seconds each:

1. **Title** — what ShiftScope does in one sentence.
2. **The problem** — the two hard questions: is the shift real and how big, and which hit to test.
3. **How it works** — the pipeline (Load, shared PCA, distances, localize, drivers, interpret)
   and the methods: energy distance + E-test, MMD, Sinkhorn OT, Leiden + Fisher/BH, Wilcoxon DE.
4. **One example, end to end** — Kang ctrl vs IFN-β, with the real numbers and driver genes.
5. **From metric to decision** — the Marson screen, PubMed grounding, Claude's keep/skip verdict.
6. **Does it hold up** — calibration limits, validation across four domains and two screens.
7. **Close** — the one-line pitch, the stack, the repo link.

The deck now carries the technical detail on the slides; keep the spoken lines short so the
whole thing stays under 3 minutes.

---

## Checklist
- [ ] Merge PR #4 (rendered notebook + no-restart setup) and PR #5 (these docs)
- [ ] Make the deck, record narration, export MP4, upload to YouTube or Loom
- [ ] Paste the summary above into the CV form
- [ ] Repo link: https://github.com/aenorhabditis6/claude_science_hackathon (public, MIT)
- [ ] Submit before Mon Jul 13, 9:00 PM ET
