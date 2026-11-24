# Fixes: close-out

Open items for the next fix point or a running lane's fix round (the loop, Each step 7).

- **Fixed by CL2** (`a0397521c`). **Look 24, the Contacts legend entry** (the fifth slice's novice pass; CL1's brief, item 4, not landed: its code is in
  `ext-rondo/src/card.ts`'s `legendRondo`, CL2's files). The entry has no swatch, and its label doubles ("Contacts: by
  chromosome"). Fix: drop `<name>: ` when it repeats the legend's title, and give "by chromosome" a hue swatch. Kind:
  look. To CL2's fix round.
- **Fixed by the orchestrator.** **The Antibody and Complex index pages in a serif** (the author, 2026-11-24: "the index
  pages feel weird being the only ones with serif font"): their root style read `--msui-*` tokens with no fallbacks,
  undefined on an index page, so the whole `font` declaration was dropped; fallbacks added (system-ui, as the other
  five). Kind: look.
- **A legend row reads as one category and hides the whole track** (the author, 2026-11-24, Features › Structure link
  › Focus And Back, two screenshots: "there is a bug on the legend here, showing just the ATP, but also hiding
  everything"): the row "Sequence ⊙ ATP" lists one swatch while the track draws its letters and eight ligands; its eye
  empties the track; the header truncates to "S…" with a stray square. Kind: bug (maybe a design gap). To FX29.
- **A bridge's ribbons missing on the canvas after its side moves** (the author, 2026-11-24, Features › Across
  sections › One Side Zoomed In, a screenshot: "there is a draw error on the canvas version of this"): the canvas leaves
  ribbons out over part of the zoomed side (Pp01's at the whole genome; 1.55–1.58 Mb of Pp02 after a Ctrl+wheel, the
  orchestrator's reproduction without FX25 and CL8), the SVG twin draws them; it stays. Kind: bug. To FX28 (bisect,
  the class fixed, a canvas-against-SVG sweep over the across-sections stories).
- **Back from deferred** (the author, 2026-11-24, the plugin off: "maybe we can finish FX25 too then?"); was **Deferred** ("ok let's just defer the fixes"): FX25 stopped three times on the
  `cc-plugin-privacy` plugin's refusals (`gaps.md`). Its six items below (TP53-201's lane, the variant card's title, stub
  labels, F9's canvas pressed slot, the square's failure text, the zoom hatch) carry past this close-out, to the next
  plan's first fix lane. Its worktree (`.claude/worktrees/agent-a98a4bf365b8a5323`, `seq/fx25` at `37832382d`) keeps
  items 1 and 2 uncommitted and unreviewed: item 1's tests shown failing on the parent by the author's own run (4 failed,
  149 passed); item 2's not. Kept, never deleted without the author's word.
- **Fixed by FX25** (`9040adc8d`). **TP53-201 is first only in the SVG** (the orchestrator's look at CL4, Features › Genome formats › GFF3 › Tp53): the
  canonical transcript is bold in both copies, but on the canvas it is the third lane (under TP53-203 and TP53-227) and
  first in the SVG. Answer F10 says first in its gene's stack. Likely the TODO's "the canvas and the SVG of one spec
  lay out differently". Kind: bug. To a fix lane, with a twin test comparing the two copies' lane order.
- **Fixed by FX25** (`9040adc8d`). **A variant's card title is its raw id** ("20001:601:f95ba8e6>C"; the orchestrator's look, Examples › Genomes ›
  Views › Samples): before CL4, not its change. The title names the variant (its kind and position, or its ID). Kind:
  look. To a fix lane.
- **Fixed by FX25** (`9040adc8d`). **Stub labels drawn over each other** (the author, 2026-11-24, Features › Genome Ring › Kinds Of Links, the ring
  zoomed to chr2 from the square: "also a visual bug", a screenshot): two stubs at one angle draw their labels on the
  same spot ("chr11 · 335" and "chr1 · 3…" overprinted into one smear), and their wedges overlap. A stub's label needs
  the ring's label collision rule (stagger or drop the weaker with a "+n"), as the chromosome labels have. Kind: bug.
  Next fix point (a fix lane, with a test on the zoomed ring's label boxes not overlapping).
- **Fixed by FX25** (`9040adc8d`). **F9's pressed slot is not drawn on canvas headers** (the orchestrator's look at CL10, Features › Structure link ›
  Color By This Track; CL10's review): `pressed` reaches the DOM header parts (`aria-pressed`) but not the canvas
  header's slots (`headers-draw`), so the story the decision came from shows no pressed icon. Kind: bug. To a fix lane:
  the canvas slot drawn pressed at rest, one at a time, a test on the canvas draw list.
- **The pressed state refreshes only on a header run** (CL10's report, in the viewer spec): a host whose colouring
  changes outside a command is not seen until the next layout, hover change or header run. Kind: followup (TODO at
  the close).
- **Fixed by FX25** (`9040adc8d`). **The square's failure text found by its source alone** (CL3's size trim, `ext-rondo/src/viewer.ts` `HoleIssue.of`
  with no section): the square passes no section, so two viewers with one source address and different failures show
  the first viewer's text on both squares (the class of CL2's fixed hole bug). Pass the square's own ring's section or
  viewer. Kind: bug. To FX25, with the two-rings test extended to their squares.
- **Fixed by FX26.** **The legend's eye chunk lands after a test file ends** (the orchestrator, CL7's merge: the quick lane exit 1 on 12
  `EnvironmentTeardownError`s, every test passing, green on rerun; CL7's lane saw it once in themes-stories): the eye's
  chunk is loaded on a legend's first mount, cached per worker, so whichever story test mounts a legend first may end
  before it lands. Kind: test flake (the class: a lazy chunk not awaited by the file that starts it). Await the chunk
  where a legend is mounted in tests (or the mount returns a promise the story helpers await); a shuffled run
  (`--sequence.shuffle.files`) as the proof.
- **Fixed by FX26.** **A ring nobody follows narrows a hidden window** (CL12's profile, Features › Genome Ring › The Genome At The Fit):
  Cmd/Ctrl+wheel takes the zoom and moves the section's viewport (0–3.03 G to 568 M–2.99 G) while the ring draws the
  whole genome, against `layouts.md` (Rings: "its zoom and pan are idle"). Kind: bug, in every option of
  `drift/ring-cmd-wheel-at-fit.md`. To FX26; whether the idle wheel is taken is the drift's.
- **Fixed by FX27** (`938396a33`: the bound is the lookup's growth across 10⁴ to 10⁷). **The 1M-range lookup's ratio fails on a quiet machine too** (the orchestrator, 2026-11-24 05:38, load 2.2, `pnpm
  test:heavy packages/seq/` alone at `77c2d3e98`: every other heavy file green): `core/src/data/_test/intervals.test.ts`
  "the lookup over 1,000,000 ranges within 1.5× of 10,000", best 0.555, 0.708, 0.677 ms against 0.368, 0.400, 0.385 (1.51×,
  1.77×, 1.76×). Kind: a bound, not a regression (it was marginal under load since GATE 1). The TODO's rule: measure it at
  a larger scale (a lookup's growth from 10⁴ to 10⁶ to 10⁷, log-like or not) rather than raise it. To the next fix lane.
