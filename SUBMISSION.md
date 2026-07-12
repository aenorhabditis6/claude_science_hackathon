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

**Slides (about 30 seconds each, 6 slides, plain talk):**

1. **Title.** "This is ShiftScope. You give it two groups of cells and it tells you how they
   differ, where, which genes, and what to test next."
2. **Works on many problems** (show `fig_domains.png`). "The same code runs on immune
   signaling, infection, blood development, and cell types, and gets the right biology each
   time. Nothing is hard-coded to one disease."
3. **One example, start to finish** (show `fig_kang_umap.png`, then the driver list and
   Claude's note). "Control versus interferon. A clear shift, it lands in specific clusters,
   and the top genes are the known interferon response. Then Claude explains it, and every
   claim points to a number we computed."
4. **The useful part** (show `fig_prioritize.png`). "On a real T-cell CRISPR screen, we rank
   hits by how strong they are and how little they have been studied. We count papers on
   PubMed, then Claude gives a keep-or-skip call using that count. Famous genes drop. Strong,
   under-studied genes rise to the top as things worth testing."
5. **Does it hold up** (show `fig_calibration.png`). "We checked where it stops working, by
   cell number and by effect size. And the same ranking recovers a second real screen."
6. **Close.** "ShiftScope turns 'these look different' into 'test this gene next.' It runs in
   Colab, it is open source. Thanks."

---

## Checklist
- [ ] Merge PR #4 (rendered notebook + no-restart setup) and PR #5 (these docs)
- [ ] Make the deck, record narration, export MP4, upload to YouTube or Loom
- [ ] Paste the summary above into the CV form
- [ ] Repo link: https://github.com/aenorhabditis6/claude_science_hackathon (public, MIT)
- [ ] Submit before Mon Jul 13, 9:00 PM ET
