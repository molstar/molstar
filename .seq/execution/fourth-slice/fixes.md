# Fixes: the author's findings and the looks' notes, collected

Findings from use that need no new design collect here, and a fix lane (FX1, FX2, …, fast unless the cause is in a
model) takes the open ones at each fix point in [plan.md](plan.md): before stage 1, after each stage's gate, and before
M4. A finding that needs a design goes to a stage 0 style draft instead; one that blocks the author's use jumps the
queue. Each line: what is seen, where (the story), where it came from, its kind (bug, look, missing design).

## Open

- **Ring letters off their cells** (look): on a ring the letters drift from the cells they name, unevenly spaced and
  crowding in places (UUUAG near 10, GGGAG near 20), not centred on their cell's angle and the track's mid radius. Where:
  tRNA › Views › Chords (and every ring with letters). The author, 2026-11-20. Check the placement against the cell's
  centre angle and the glyph's centre (not its baseline start), in the canvas and the SVG.
- **The Reverse track reserves an empty label lane** (look): Plasmid › Ring › Detail with a strip. Third slice's DS look.
- **Dense leaders when a ring group opens** (look): Plasmid, a group of more than 6 labels hovered. Third slice's DS look.
- **The Antibody's IMGT label reads as a second line** (look): Antibody, the ruler's numbering label. Third slice's look.
- **A moved label of an item in a lower lane is dropped** (bug): leaders exist for lane 0 only (LX). Third slice's LX
  review.

- **Alignment Rows And Columns collapses when the frame's size changes** (bug; look if the cause is only the story's
  sizing): Features › Query › Alignment Rows And Columns. When the Storybook zoom slider changes the frame's size, the
  plot collapses: at a small zoom the rows shrink into the top of a much larger, empty canvas; at a middle zoom the rows
  end two thirds down with an empty well below; at 100% the last column's letters are cut at the right edge under the
  end numbers (N, (, A half drawn). The status line under it says "The alignment's rows have not arrived" while the
  rows are drawn. The author, 2026-11-20. Check: the plot box against the host's size on resize (grids-and-rows.md: the
  plot box is fixed); the canvas's backing size against its CSS size at a changed devicePixelRatio or zoom; the right
  edge's clip against the end-numbers column; the status text's settle.
- **Two frames of about 100 ms in every zoom of the 10,000 × 300 alignment** (bug): Globins › Rows › Scale at dpr 2. A
  frame callback outside the viewer's ticks, most likely the rows' image tiles rebuilt at a new level. Build them
  across frames and keep the old level drawn meanwhile (performance.md, Levels). From the third slice's M3
  (`third-slice/measurements.md`).
- **A label that does not fit its box disappears whole; an ellipsis inside first** (look, a rule added to layouts.md
  "Labels inside if they fit, else beside"): p53's domain track at a mid zoom ("Nuclear localization", "C-terminal",
  "Tetramerization" shown, then the boxes empty a step later, C-terminal cut at the edge). The author, 2026-11-20:
  "should we allow ellipsis (up to a degree) in the boxes instead of the text disappearing as a whole?" Proposed order:
  whole inside; else whole beside when the lane has room; else inside cut with an ellipsis while at least 4 characters
  (or the whole label, if shorter) and the ellipsis fit (4 by default, configurable: a layer param beside the
  label's, `labels.minChars`, decided by the author, 2026-11-20: "4 characters as default is fine, should be
  configurable"); else dropped. The cut keeps the hover card and the address
  whole; the canvas and the SVG cut alike (text-metrics port). A label cut at the section's edge is placed within it
  (C-terminal), never clipped mid-glyph.
- **The vertical overview went missing from the stories and examples** (bug; first in FX1): the rows rail was to be added
  beside the vertical overview (the author: "also have the thin slider"), but it replaced it. Globins › Rows › Folded
  and the other alignments now show only the rail. The author, 2026-11-20: "at the very least we should have stories for
  this in features/navigation & use it in some examples too". Fix:
  - Features › Navigators: one story for each vertical form side by side (the rows rail, the vertical overview with its
    four handles as in `packages/sequences`, both together);
  - the vertical overview declared again in the examples that had it (Globins Family, Folded and Scale; any other
    alignment with more rows than fit);
  - the rail stays the default only where no navigator is declared (navigators.md, the rows rail);
  - check that the overview's picture and handles still draw and drag (`render/navigator.ts`, `mount-navigator.ts`).
- **A folded clade's summary row is not named to assistive technology** (bug, only if the author decided item 4 of the
  third slice's `drift/ec-contradictions.md` as recommended: name it by its summary label, "HBA +12"). Third slice's EO.
- **A disabled toolbar button's dimming not looked at** (look): `@molstar/ui`'s Storybook and Prism, light and dark
  (computed at 1.7-1.9:1, not seen). Third slice's SF leftover.
- **The frame-time test's max bound flaked once** (bug in a test): `render/_test/frame-time.test.ts`, 18.8 ms against
  16, passing on the rerun; best of three or a percentile. Third slice's gate.
- **A cursor card is replaced when a zoom re-hits the hover** (bug): the card should stay while the hit is the same.
  Third slice's HC leftover.
- **A ruler naming its chain stays labelled; an empty track label is still drawn as empty text** (look): the Complex and
  the Antibody's against view (rulers); core `tracks.ts` (the label). Third slice's HC leftovers.
- **An inline `strand` of null breaks core's model** (bug): `strand.get is not a function`; a null reads as no strand.
  A section over a source sequence needs the sequence named after the source's system, else its length is unknown and
  nothing lays out, and nothing says so: an issue naming it. Third slice's SC leftovers.
- **`seq.resolve.length.unknown` before open** (look): quieter in the viewer, which says it while a source opens.
  Third slice's CY leftover.
- **Antibody: IMGT 112.1 never on the ruler at the closest zoom** (three insertion labels a residue apart); short CDR
  labels at the fit view (H1 by Kabat, L1 and L2 by IMGT); the against view's gaps faint at the fit (look). Third
  slice's DA leftovers.
- **Not looked at past the lanes' own screenshots** (look, for the use gate's look): the Antibody's insertions at the
  closest zoom and the rail's row (CF), the rows layer (RX).

## Fixed

(none yet; each line moves here with its lane and commit)
