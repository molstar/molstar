# Navigator window look: one control row, one stroke, one readout place (an amendment to navigators.md)

Lane DL, the fifth slice, 2026-11-22. A look design before a fix lane. It amends [navigators.md](../../../design/navigators.md)
("The window", "The zoom strip", "Light and dark", the matrix and Configuration) and
[genome-coordinates.md](../../../design/genome-coordinates.md) ("The zoom strip over several systems"); F25 (the
fourth slice's novice finding) is decided here with it.

## Intent

The author, 2026-11-21, on a ring that switches chromosomes: "i don't understand why chr6 gets a 'box' and the other
'sub-section'"; "the 'subsection visual' is really ugly, with different line thicknesses". And 2026-11-22, the BED
cytobands story: "it's better if the label stays outside in top right instead of jumping around all the time. maybe we
can do the same thing for the 'vertical maps' as well, where it appears 'above' the map instead."

What a reader must be able to say after one look: where the control is, which mark is "where I am" and which is "what
I drag", and where to read the numbers. Today none of the three is stable.

## What is drawn today, and why it reads as two designs

- **The strip comes and goes with the window's width.** Over blocks, `zoomFlips` takes the window's own width
  (`core/src/render/navigator.ts:657`, `px = this.overBlocks ? this.windowPx : …`, `windowPx` set at `:967`) and turns
  the zoom strip on under `navigator.zoomStrip.minPx` (48, `window.ts:28`), off past it plus 8 px of hysteresis
  (`navigator.ts:94`, `:659`). chr6 whole is about 85 px on the genome row (a box with grips); chr13 whole is under 48
  (a strip). The trigger is a pixel count no reader can see. On one system it reads the reachable span instead (a fixed
  property of the section), which is the rule navigators.md states.
- **Four stroke weights.** The window's outline 1 CSS px, 2 when lit (`window.ts:247`); the context mark an outline of
  1 CSS px with no fill (`navigator.ts:1142–1168`, `drawNavigatorWindow` with `FlagOutline`); the connectors one
  device pixel (`g[6] = 1 / dpr`, `navigator.ts:1126`), so half the window's line at DPR 2; the idle outline 1 px at
  0.4 (`window.ts:36`). The lit 2 px also makes the box jump a pixel under the pointer.
- **The mark and the window share an ink but not a role**: the mark is drawn as a window without its fill
  (`seq-window`, full strength), so the context's "where" mark is as loud as the control.
- **The readout moves** (`viewer/src/readout-place.ts:82–122`, `placeReadout`): top right, else top left, else beside
  the windows, else above or below the part, chosen per frame from the window boxes and the pointer
  (`viewer/src/navigator.ts:376–378`, `mount-navigator.ts:89–95`). It jumps as the window or pointer crosses it.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| IGV: the genome or chromosome ideogram with its red box, the locus box above the tracks | the context row as a quiet "where" mark; the readout in one place, read at rest | a red, loud box (here the window is chrome, never louder than the data) |
| UCSC: the ideogram above the browser image, the position box at the top, unmoving | the numbers in a fixed line above the picture | — |
| JBrowse 2: the overview scale bar with the region highlighted and polygon lines to the view | the funnel as a thin, muted tie between two scales | its filled polygon (it covers the ideogram's bands) |
| This codebase: the zoom strip (navigators.md), the floored window, idle drawn quiet (EC-110) | the strip's parts, gestures and record; the floor | the window-width trigger over blocks |

## Decision

### 1. One rule for the second row: by what the section reaches, never by the window's width

- **A navigator has a control row and, sometimes, a context row above it.** The control row holds the window with its
  grips; the context row holds only a quiet mark of where the control row lies in the whole. A reader learns one thing:
  *the grips are always on the bottom row*.
- **When there are two rows** (one rule, readable from the spec, never from a pixel count):
  1. **Over blocks** (`navigator.system` an assembly): **always two rows**. The context row is the genome (each
     chromosome's bands, click to jump, sweep across blocks: unchanged); the control row is the hull of the blocks the
     window touches (one chromosome when the view is within one; the genome's ruler at the fit). So chr6 and chr13
     look the same: both the chromosome below, a quiet mark above.
  2. **One system whose reachable span is narrower than `navigator.zoomStrip.minPx`** (a region on its chromosome, p53's
     6 kb on chr17): two rows, as built. The span is the section's, so this never flips with a zoom.
  3. **Otherwise one row**: the window on the whole axis, its grips there; a window narrower than the floor is floored
     (16 px, centred, its hairline at the truth), never handed to a strip.
- **The height never changes with a gesture**: in cases 1 and 2 the part's height is the context plus the strip, fixed
  for the navigator's life (no hysteresis needed: `ZoomHysteresisPx` goes). Why always over blocks: the genome's and the
  chromosome's scales are two (genome-coordinates.md's reason), and a toggle by width is a rule the reader cannot see.
- **The control row at the fit is the genome's ruler** (case 1 with the whole genome in view): Mb ticks under the
  ideogram, not a blank well (the empty-well rule's "a strip with no content draws its system's ticks"), its window at
  the fit (2. below).

### 2. A window at the fit is one quiet frame (F25)

- **At the fit** (the window covers the whole control row, whether or not zoom-in is possible) the window is drawn as
  the idle frame: an outline in `seq-window` at the quiet alpha, no fill, no grips, no wash. It is not a drag target;
  **the whole row is a sweep**: a drag draws the region to zoom to (`navigator.zoom-to`), a double-click zooms there
  (`navigator.zoomPx`), the cursor `pointer` (the `outside` part's). Why: at the fit a box with grips at the row's own
  edges says "drag me" and does nothing visible (F25: "nothing"); a sweep is the one gesture that moves the view from
  there.
- **Idle** (cannot move at all) stays as built: the same quiet frame, no sweep, the default cursor.

### 3. One stroke, two inks, one fill

| Mark | Stroke | Ink | Fill | Grips |
|---|---|---|---|---|
| the window (control row), at rest | `navigator.strokePx` (1 CSS px), snapped to the device grid | `seq-window` | `seq-window` at 0.12 | `seq-grip` tabs |
| the window, lit (hover, drag) | the same 1 px (no 2 px step) | `seq-window` | `seq-window` at 0.2 | the lit grip(s) in `seq-window`, one px wider (as built) |
| the window at the fit, idle | 1 px | `seq-window` at 0.4 (quiet) | none | none |
| the context mark (where the control row lies) | 1 px | `seq-window` at 0.4 (quiet) | none | none |
| the funnel (two lines, mark's ends to the control row's ends) | 1 px | `seq-window` at 0.4 (quiet) | — | — |
| the wash outside the reach | — | — | `seq-wash` | — |

- **One stroke weight** for every line of the part: the window, the mark, the funnel, an idle frame, a floored window's
  hairline at the truth stays a device pixel (it marks a sub-pixel truth, not an outline). Lit is said by the fill and
  the grips, never by a thicker line (the box no longer jumps a pixel).
- **Two inks from one token**: the control at full `seq-window`, every "where" mark at the quiet alpha (`IdleAlpha`,
  0.4, made a `RenderConfig` value). So the loudest mark of the part at rest is the data (bands, a miniature's marks),
  then the window, then the context: the defaults audit's "the loudest mark at rest", kept.
- **The funnel kept, thin and muted** (it ties a mark 16 px wide to a row 600 px wide; without it the two rows read as
  two navigators), from the mark's lower corners to the control row's upper corners, never filled.
- **The control row has no frame of its own**: its edges are the picture's; the gap (`ZoomGapPx`, 6) holds the funnel.

### 4. The readout in one place, outside the part

- **Above a horizontal navigator, right-aligned at its top-right corner**, in a line of its own
  (`navigator.readout.linePx`, 16, the `msui-text-xs` leading). **Above a vertical navigator** (the vertical overview),
  right-aligned to its right edge, in the line above its picture. It never moves with the window or the pointer, so it
  never covers the control in use.
- **The line is the part's own**: the mount lays it out inside the part's box, the canvas under it, so a story's or a
  page's top edge cannot clip it and nothing needs room above the part. A toolbar line (the viewer's, when the host
  places the navigator under it) is not borrowed: two owners of one line would fight on width.
- **What it says**: at rest, the window's region (`chr17:20,145,469–27,700,018 of chr17`; several leading windows:
  the one last acted on); while hovering, what is under the pointer (position, band); while a gesture runs, the
  window's region as it moves, every frame. The value never truncates; the system's name gives way first (as built).
  Always shown at rest, so the line is never a blank band and the part never changes height (Gaps, 1).
- **A vertical overview is narrow (96 px, 48 in a narrow frame)**: its readout takes the short form, the rows only
  (`401–440 / 10,000`), right-aligned; past the part's width it may extend left over the label-free gap above the rows
  band, never over data (Gaps, 3). The rows rail and the ruler rail drawn inside the view have no readout (as built).
- **A ring**: a ring's window is not a navigator part and has no readout (the hover card names a position; the hole
  holds the title). A navigator part bound to a ring section (a strip over the plasmid's system; Rondo's square beside a
  genome ring) follows the horizontal rule: above, right-aligned.
- **Removed**: `placeReadout`'s candidate search, `windowBoxes`, `pointerBox`, `GripReach`, `PointerReach`
  (`readout-place.ts`, its test `viewer/src/_test/readout-place.test.ts` replaced by one on the fixed place), the
  mount's per-frame `readoutAt` (`mount-navigator.ts:91–95`) and its "no `overflow: hidden`" exception (`:36`): the
  readout is inside the part now. FX5's "beside the window, else above or below" goes with them.
- **A bug for the fix lane, not a design question**: the readout lags a gesture (the window at 21 Mb while it reads
  `1–83,257,441`). The text is set only from a hover's hit: a move with a press held returns before `say`
  (`viewer/src/navigator.ts:468–472`), and wheel, keys and moves from elsewhere never call it. Fix: while a window moves
  (any source), the readout reads the moving window's text from the frame each draw; on release the hover's text.

## Defaults: what a reader sees

Over blocks, a chromosome in view (chr6 or chr13 alike):

```
                                              chr13:1–114,364,328 of chr13   <- readout line, fixed
 ▐▌▐▐▌▌▐▌ ▐▌▐▐▌▐ ▐▌[▒]▐▌▐ ▐▌▐▐ ▐▌ … ▐▌▐▐▌         <- genome: quiet mark (1 px, 0.4), no grips
                   ╱  ╲                           <- funnel: 1 px, 0.4
 ╱                       ╲
▕▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▐▌▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▌▐▐▌▐▌▐▐▌▐▌▌▐▏   <- chr13 whole: at its fit, one quiet frame
```

Zoomed into chr13 (the control row's window has its grips; the genome mark floored to 16 px):

```
                                        chr13:48,300,000–49,100,000 of chr13
 ▐▌▐▐▌▌▐▌ ▐▌▐▐▌▐ ▐▌[▒]▐▌▐ ▐▌▐▐ ▐▌ … ▐▌▐▐▌
 ╱                 ╱  ╲                       ╲
▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▐▌▐▌▐▐▌▐▐ ▌▐[███]▌▐▐▌▐▌▐▐▌▐▌▌▐▌▐▐   <- window: 1 px seq-window, 0.12 fill, ▐ ▌ grips
░░░░░░░░░░░░░░░░░░░░░░░░░░░    ░░░░░░░░░░░░░░░░   (wash outside the window)
```

One system, the whole axis reachable (a protein's strip): one row, a box with grips; zoomed in deep, a floored box:

```
                                                    p53:120–180 of 393
▕───────────────[▒▒▒]──────────────────────────────▏  floored 16 px, a hairline at the true 60 residues
```

A cyclic axis, the window over the origin (the plasmid's strip, pUC19c, 2,600–60 of 2,686): two pieces, a grip at
each true end only (`noStart`, `noEnd`: as built), the wash between them:

```
                                              pUC19c:2,601–60 of 2,686 (across the origin)
▕██▌░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▐██▏
  ^ end grip                          start grip ^
```

The vertical overview: the readout above it, short form, right-aligned:

```
            401–440 / 10,000
           ┌────────────┐
           │▓▓▓▓▓▓▓▓▓▓▓▓│
           │▓▓[      ]▓▓│  <- window: rows and columns, four grips
           │▓▓▓▓▓▓▓▓▓▓▓▓│
```

## The interaction matrix's cells

Changed cells only; the rest stand as navigators.md gives them.

| Row | Cell | Was | Now |
|---|---|---|---|
| a navigator | zoom in / out | drag a grip; sweep from outside | the same; **at the fit, a drag anywhere on the control row sweeps** (no grips there) |
| a navigator | pan | drag the body | the same; **at the fit: none** (the frame is not a target) |
| a navigator | focus | click outside centres; double-click | the same; at the fit, double-click zooms to `zoomPx` there |
| a navigator over blocks | focus | click on a chromosome jumps (context row when the strip shows, else the row) | click on the context row's chromosome jumps, **always the context row** |
| a window | zoom, pan | drag `start`, `end`, `body` | the same off the fit; at the fit quiet, not a target |
| the readout | (no row) | — | not a target (`pointer-events: none`, as built); Enter at the readout (navigators.md's focus cell) still centres on the readout's position |

Negations checked: a window narrower than the floor (floored, no strip: case 3); a window across two blocks (the
control row their hull); the fit on a single system with nowhere to zoom (idle, no sweep); a reversed section (the
control row in display order, as built); a follower's window (moves its leader, as built).

## Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `navigator.zoomStrip.minPx` | `ViewerConfig` | 48 | **now read on the reachable span only**, never the window: a span too small to drag gets the strip | any px |
| `navigator.zoomStrip.overBlocks` (new) | `ViewerConfig` | `always` | the visible rule: over a genome the chromosome is always the control row | `narrow` (today's: by the window's width, with hysteresis), `never` (one row, floored window) |
| `navigator.zoomStrip.heightPx` | `ViewerConfig` | 16 (over blocks at least the ruler's height, as built) | — | px |
| `navigator.zoomStrip.funnel` (new) | `RenderConfig` preset | `muted` | ties two scales; quiet | `none` |
| `navigator.strokePx` (new) | `RenderConfig` preset | 1 | one weight for every line of the part | px |
| `navigator.window.quietAlpha` (new; `IdleAlpha` today, fixed) | `RenderConfig` preset | 0.4 | the context mark, the funnel, an idle or fitted frame | 0–1 |
| `navigator.readout.place` (new) | `ViewerConfig` | `above` | one place, never over the control (the author) | `none` (a host draws `readoutParts` itself) |
| `navigator.readout.linePx` (new) | `RenderConfig` preset | 16 | the `msui-text-xs` leading | px |
| `navigator.readout.rest` (new) | `ViewerConfig` | `window` | where am I, read without a hover | `none` (shown only while hovering or moving; the line stays reserved) |

## Tensions

- **Two rows over a genome at all times cost 16–22 px** that the old rule saved on a large chromosome. Taken: a fixed
  height and one place for the grips beat the saving; `overBlocks: narrow` keeps the old for a host short of room.
- **The fit's sweep** removes "drag the box" at the fit; a reader who drags at the fit gets a zoom, not a pan (there is
  nothing to pan to). Taken.
- **A readout line always reserved** costs 16 px a navigator even with `rest: none`. Taken: a line that comes and goes
  moves everything below it.

## Self-check

(a) The trigger is now visible: over a genome, always; on one system, by the section's reach. (b) Inverse: leaving a
chromosome for the whole genome keeps both rows (the control row becomes the genome's ruler at its fit). (c) Symmetric:
the vertical overview's readout above it as the horizontal's is above. (d) Plural: several windows on one control row
(multi-locus) each drawn by §3, the readout the last acted on. (e) Cyclic: two pieces, grips at true ends. (f) Export:
navigators omitted by default (EC-138, unchanged); with `include` the window an outline at `strokePx`, no readout line
(chrome state). (g) Forced colors: `CanvasText` for every line at one weight, no wash, no alpha (as built). (h) DPR 1
and 2: one CSS px snapped, so the window and the funnel match at both.

## Consequences

- `core/src/render/navigator.ts`: `zoomFlips` keyed on `overBlocks` and the reach only (`:646–670`), the hysteresis
  removed; `drawContext` and `drawMark` at the quiet alpha and `strokePx` (`:1096–1168`); the fit drawn as a quiet frame
  with the sweep as its part (`:1071`'s flags gain a `fit` case).
- `core/src/render/window.ts`: lit by fill and grips only (`:247`); `IdleAlpha` to config.
- `viewer/src/readout-place.ts` removed; `viewer/src/navigator.ts` readout from the frame while a window moves;
  `mount-navigator.ts` lays out the line above the canvas.
- navigators.md: "The zoom strip" and its Defaults picture, "Light and dark" (the stroke table), the matrix and
  Configuration rows above; genome-coordinates.md: "The zoom strip over several systems" (always), its Defaults 8.
- Tests: a twin test (canvas, SVG with `include`) that every line primitive of a navigator has one width; a frame test
  over chr6 and chr13 (two rows, grips on the bottom row only, equal part height); F25: at the fit no grip primitive and
  a drag sweeps; the readout's box fixed across a drag and its text equal to the window's every frame.

## Gaps for the author

1. Show the window's region in the readout at rest (not only on hover)? **Recommend yes**: the line is reserved anyway,
   and "where am I" without a hover is IGV's and UCSC's locus box.
2. Always two rows over a genome (chr6 and chr13 alike), rather than a strip only for small windows? **Recommend yes**
   (`overBlocks: always`); `narrow` stays for a host short of height.
3. On the narrow vertical overview, a short readout (`401–440 / 10,000`) allowed to run left past the part over the
   gap above the rows band? **Recommend yes**; the fix lane measures that the gap is data-free in the rows stories.
4. Keep the funnel (thin, 0.4) rather than drop it? **Recommend keep**: without it a 16 px mark and a full-width row read
   as two navigators; `funnel: none` for a host.
5. At the fit, make the whole control row a sweep (no grips, drag zooms)? **Recommend yes** (F25).
6. Drop the 2 px lit stroke (lit by fill and grips only)? **Recommend yes**: one weight, no pixel jump under the pointer.
