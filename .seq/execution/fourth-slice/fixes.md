# Fixes: the author's findings and the looks' notes, collected

Findings from use that need no new design collect here, and a fix lane (FX1, FX2, …, fast unless the cause is in a
model) takes the open ones at each fix point in [plan.md](plan.md): before stage 1, after each stage's gate, and before
M4. A finding that needs a design goes to a stage 0 style draft instead; one that blocks the author's use jumps the
queue. Each line: what is seen, where (the story), where it came from, its kind (bug, look, missing design).

## Open

- **Hover cards too loose** (look): every view's card (Globins family: "Row, columns, index"; a ring's features). The
  author, 2026-11-20: "the hover cards should be more compact in general". Tighter row spacing and padding, the key and
  value columns closer; the card's title once (the Globins card repeats the row's name as title and first row).
- **A hover card with many items lists every one** (missing design): a ring of 200 kb, the card at a position lists
  20+ features down the page. The author, 2026-11-20: "if there is too many items, should truncate? (app level
  overridable ... we should have ports anyway)". Proposed: the default card shows the first n items (a `hover.card`
  param, by count or by the card's height) and "+k more"; the full list stays in the hover record; a host replaces the
  card through the existing `hoverCard` port. State it in viewer-models.md (or ux-foundations.md's hover) and its page
  with the fix.
- **The ring's card lists features that do not cover the hovered position** (bug, to confirm): the same card at
  131,892 lists features whose ranges wrap the origin (199,690-200,000, 1-24), none containing 131,892. Check the
  hit or the card's items for wrapping ranges (an `end < start` range read as its complement, or the whole ring).
  Seen in the author's screenshot, 2026-11-20.

- **The vertical map's resize handle scrolls the view instead of resizing** (bug): the Globins family view (the
  alignment's rows by species, its vertical map with the window and the thin slider beside the rows' labels). Dragging
  the window's vertical resize anchor (the bars on its top and bottom edges) scrolls the rows instead of changing the
  window's height (the rows' zoom). The author, 2026-11-20: "the vertical resize anchor doesn't work on this and
  instead scrolls the view". Check: the anchor's hit target against the map's pan target, and the gesture's
  precedence (navigators.md).
- **A double-click on the vertical map's window, and on the thin slider, resets the vertical zoom** (missing design,
  small): the same view. The author, 2026-11-20: "double click on the window (and perhaps the 'thin slider') should
  prolly reset the vertical zoom too?". As the horizontal navigator's double-click resets (navigators.md); state it in
  the design and its page with the fix.

- **The core folder's quick-lane time over 2 s** (waste or nature): `packages/seq/core` at 2.7-3.4 s before and after
  MR; bigWig 1.0 s, hover-time 0.7 s, mmCIF 0.7 s (MR's review). Fix the waste, move what is slow by nature with a
  quick representative.
- **The viewer's quick-lane folder over 2 s** (waste): `packages/seq/viewer` at 2.62-2.68 s on a quiet machine
  (load 2.3) before and after HP1, about a third imports and a fifth jsdom's environment. Measured at HP1's merge. Check
  which files set up jsdom or import the full viewer per file; share a fixture per file, or move a slow-by-nature test
  with a quick representative kept.

- **SZ's viewer split, four small findings** (low): `viewer.test.ts:116` checks the `_viewport`/`_renderer` shims, not the frame's own calls; `viewer-motion.ts`'s viewport assembly (`window`, `all`, `_bandsFor`) belongs in `viewer-viewports.ts`; `ViewerFrame`'s mutable public fields written by `useBackend`; `view-follow.ts`'s shared out-records. SZ part two's review.

- **The Antibody's IMGT label reads as a second line** (look): Antibody, the ruler's numbering label. Third slice's look. Open: not part of FX2; whether the ruler's name or the `labels` lane makes the second line is not yet read.
- **Alignment Rows And Columns collapses when the frame's size changes** (bug; look if the cause is only the story's
  sizing): Features › Query › Alignment Rows And Columns. When the Storybook zoom slider changes the frame's size, the
  plot collapses: at a small zoom the rows shrink into the top of a much larger, empty canvas; at a middle zoom the rows
  end two thirds down with an empty well below; at 100% the last column's letters are cut at the right edge under the
  end numbers (N, (, A half drawn). The status line under it says "The alignment's rows have not arrived" while the
  rows are drawn. The author, 2026-11-20. Check: the plot box against the host's size on resize (grids-and-rows.md: the
  plot box is fixed); the canvas's backing size against its CSS size at a changed devicePixelRatio or zoom; the right
  edge's clip against the end-numbers column; the status text's settle. **FX1** (`51682eebc`) fixed the collapse and the status line. **Still open (FX2 did not touch it):** the last column's letters cut under the end numbers at 100%, not reproduced headlessly; needs a look in the browser.
- **A disabled toolbar button's dimming not looked at** (look): `@molstar/ui`'s Storybook and Prism, light and dark
  (computed at 1.7-1.9:1, not seen). Third slice's SF leftover. Open: a look, not a fix lane's work.
- **Antibody: the against view's gaps faint at the fit; short CDR labels at the fit view** (look; H1 by Kabat, L1 and L2 by IMGT). Third slice's DA leftovers. FX2 fixed 112.1 (`ed8050a46`); these two were not part of FX2 and wait for the use gate's look.
- **Not looked at past the lanes' own screenshots** (look, for the use gate's look): the Antibody's insertions at the
  closest zoom and the rail's row (CF), the rows layer (RX). Open: for the use gate's look.

## Fixed

- **Decided designs for FX2** (the author, 2026-11-20; `drafts/fx-designs.md`): built, designs and Design pages in step (FD2). Fixed by FX2: `1ae3224c0` (`navigator.rowsRail`), `4fa5d5d1e`, `14a8a20ec`, `5d64807f3` (a narrow frame's order, `navigator.overview.minPlotPx`, `minPx`), `6f1952ea2`, `5f7b13c17` (`rings.cluster.fanMax`), `fbf3839b8`, `c0e55a47f` (the issues).
- **The "2 hidden" chip overlaps the vertical overview's top** (look): Globins › Rows › Folded at 1280 px. Fixed by FX2: `c3129dd89`.
- **The Reverse track reserves an empty label lane** (look): Plasmid › Ring › Detail with a strip. Fixed by FX2: `9b1d27080`.
- **Dense leaders when a ring group opens** (look): Plasmid, a group of more than 6 labels hovered. Fixed by FX2: `6f1952ea2`, `5f7b13c17`.
- **The vertical overview went missing** (bug): the rows rail replaced it. FX1 (`68ab75210`) declared it again; both together fixed by FX2: `1ae3224c0`.
- **A ruler naming its chain stays labelled** (intended): fixed by FX2 as one line in the spec, `dca490334`; the empty label by FX1 `e6487e46e`.
- **An inline `strand` of null; an unnamed source system** (bug): null by FX1 `666e3464c`; `seq.resolve.system.unnamed` by FX2 `fbf3839b8`.
- **`seq.resolve.length.unknown` before open** (look): held back while the address is pending. Fixed by FX2: `fbf3839b8`, `c0e55a47f` (`viewer.modelIssues`).
- **Antibody: IMGT 112.1 never on the ruler at the closest zoom** (look): Fixed by FX2: `ed8050a46`, `d228a4d2d`.
- **Ring letters off their cells** (look): on a ring the letters drift from the cells they name, unevenly spaced and
  crowding in places (UUUAG near 10, GGGAG near 20), not centred on their cell's angle and the track's mid radius. Where:
  tRNA › Views › Chords (and every ring with letters). The author, 2026-11-20. Check the placement against the cell's
  centre angle and the glyph's centre (not its baseline start), in the canvas and the SVG.
  Fixed by FX1: `ef7db02cf`.
- **A moved label of an item in a lower lane is dropped** (bug): leaders exist for lane 0 only (LX). Third slice's LX
  review.
  Fixed by FX1: `1e638edd4`.
- **Two frames of about 100 ms in every zoom of the 10,000 × 300 alignment** (bug): Globins › Rows › Scale at dpr 2. A
  frame callback outside the viewer's ticks, most likely the rows' image tiles rebuilt at a new level. Build them
  across frames and keep the old level drawn meanwhile (performance.md, Levels). From the third slice's M3
  (`third-slice/measurements.md`).
  Fixed by FX1: `ac41ad4f5`, `4562016a3`, `ae15f49f6`, `b615e3137`.
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
  Fixed by FX1: `7651d491a` (`labels.minChars`, a `RenderConfig` key).
- **A folded clade's summary row is not named to assistive technology** (bug, only if the author decided item 4 of the
  third slice's `drift/ec-contradictions.md` as recommended: name it by its summary label, "HBA +12"). Third slice's EO.
  Fixed by FX1: `4c4470716`.
- **The frame-time test's max bound flaked once** (bug in a test): `render/_test/frame-time.test.ts`, 18.8 ms against
  16, passing on the rerun; best of three or a percentile. Third slice's gate.
  Fixed by FX1: `b022e2330` (the third slice, best of three).
- **A cursor card is replaced when a zoom re-hits the hover** (bug): the card should stay while the hit is the same.
  Third slice's HC leftover.
  Fixed by FX1: `64f5d5ca0`, `1bdc2b9f0`.
