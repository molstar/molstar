# Track sizing: content height, a cap, then its own scroll (an amendment to layouts.md and ux-foundations' heights)

Lane DL, the fifth slice, 2026-11-22. A short design before a fix lane. It amends
[layouts.md](../../../design/layouts.md) ("A track's height by default fits its content", "A band may scroll inside its
section", Export, Rings, Configuration), [ux-foundations.md](../../../design/ux-foundations.md) section 3 (row heights)
and [grids-and-rows.md](../../../design/grids-and-rows.md) ("The band"), and closes `packages/seq/TODO.md`'s "Track
heights and labels beside for a source's items".

## Intent

The author, 2026-11-21, Features › Genome formats › GFF3 › TP53: "can we do a better job at the track sizing for the
transcripts in general? what if that section needs a scroll?"; earlier, blank lanes above and below and the first row
unlabelled. A track of many lanes should be as tall as what it holds, no taller, never a blank band, and when it holds
more than a screen can give it, it scrolls inside itself with its name and the tracks around it in place.

## What is drawn today

- **Inline data** fits its content up to five lanes (`maxLanes: 5`, `core/src/render/config.ts:374`;
  `lanesHeight`, `sizing.ts:60–63`: 22 … 78 px), then compacts (lanes squeezed to the band, `ranges.ts:155–159`) and
  past that a density row takes what does not fit (`ranges.overflow: extra-lanes`, `ranges.ts:260`).
- **A source's data** has no content (`rangesContent` returns undefined for data not known whole,
  `ranges-content.ts:56`, `:133`), so its track takes the fixed `heights.ranges` (36 px, `config.ts:360`) and its lanes
  are centred by the most it has shown (`seen`, `ranges.ts:264–269`). The TP53 story works around it with a hand height
  (`height: 760` for 47 transcripts, `examples/src/genomes/gff.ts:65`): taller than the page, blank where the region's
  models need fewer lanes, and wrong for the GTF's count.
- **The gene line drops a label whose room runs past the band's left edge** (`ranges-gene.ts:22`): the model starting at
  the region's start, often the first lane, is unlabelled.
- **The rows band** already scrolls (`bandTop`, `rowScale`; layouts.md "A band may scroll inside its section") with its
  rail (`rail.ts:36–46`, `SlimRailPx` 12): the mechanism exists; ranges do not use it.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| IGV: collapsed, expanded, squished display modes; expanded grows a lane a row | expanded as the default (a lane a row at full height); squished as an opt-in overflow | the data panel's one shared scrollbar (here each capped track scrolls on its own, the ruler and other tracks fixed) |
| JBrowse 2: a track's height its own, scrolling inside past it; a max height past which it says so instead of drawing | a cap, then an inner scroll; a plain note when a hard bound is reached | a fixed default height for every track (here content first) |
| UCSC: dense, squish, pack, full; pack falls back to dense past a count | the ladder of modes by name; the fallback to a density row past a hard bound | the mode switched silently by count (here the bound is said) |
| This codebase: content heights in lane steps (decided 2026-11-19), the rows band's scroll and rail, `text.maxHeight` and scroll chaining | all three, joined: content height, the text box's cap and chaining, the rows' rail | — |

## Decision

### 1. Content height, as known

- **A lanes track (ranges, the gene line, variants laid in lanes) is as tall as its lanes** (`lanesHeight`: padding
  4 + 14 a lane), with no five-lane ceiling: `maxLanes` stops being the content cap (§2 replaces it).
- **A source's lanes are known as they arrive**: the content is the most lanes the layer has packed over the data read
  for the section's reachable region (its tiles at the fit's level, as the navigator reads), grown as tiles settle,
  never shrunk by a pan or a zoom (the decided "heights from the whole data stay put"); reset only when the data's
  identity changes (another file, another region `range`). Why the reachable region: it is the source's "whole" for
  this view, bounded, and what a pan can reach.
- **Growth is batched to the settle**, not per tile: the track grows at most once per settled read (a store `settle`,
  data-store.md), so tracks below move once, never tile by tile. Data does not animate (ux-foundations section 7).
- **Pending is never empty** (EC-169): before the first tile, one lane with the loading placeholder; never a blank band
  of a guessed height. **Settled with no items**: the empty-well rule (`track.empty`, fixes.md) collapses it to its
  note row.

### 2. A cap

- **`track.maxPx`, a `RenderConfig` preset value, 240 px** (16 lanes and padding). Not a share of the view: a viewer in
  a page sized by its content has no view height to take a share of (a circle), a share changes the track as the window
  resizes (a pan never resizes; a resize should not either), and an export has no view at all. A px number is what
  `text.maxHeight` (600) already is. A host fitting its own panel sets it, or answers heights through the
  `trackSizing` port (unchanged).
- **A track's own `height` is its cap and its height**: under it, lanes past the height scroll (the TP53 story drops
  `height: 760` and gets content to the cap).
- **Tracks of one kind share their height** (as decided): the tallest content among them, capped; each scrolls on its
  own when its lanes exceed it.
- **`track.overflow`** (a spec param of the track, its default in `RenderConfig`): `scroll` (the default, §3),
  `compact` (the band at its cap, lanes squeezed: IGV's squished, today's step), `density` (the band at its cap, the
  rest in a density row: today's `extra-lanes`).
- **A hard bound** (`SpecBounds`, `ranges.maxLanes`, 1,000): past it the last lane is the density row and the track's
  header says "+n items in density" (JBrowse's "said, not hidden"). Why: lanes are layout work and draw lists per frame.

### 3. Past the cap, its own scroll

- **The lanes scroll inside the track's band**, by the band scroll layouts.md already defines (`bandTop` in the
  section's viewport, by band; the layer reads `ctx.bandTop`). The track header, its label in the label column and the
  tracks above and below stay fixed. Item labels are on their lanes and scroll with them.
- **Wheel**: its vertical delta over the band scrolls the lanes while they overflow, its horizontal delta pans the
  section (as over the rows band); at the top or bottom the vertical delta chains to the section, then the page (the
  text box's rule). A track that shows all its lanes takes no vertical wheel.
- **Keys**: the keyboard cursor moving to an item reveals its lane (the band scrolls the least to show it);
  Page Up and Page Down scroll the band under the cursor a band's height (the rows' rule, EC-155: only while a cursor
  is shown).
- **A slim rail**, the rows rail's geometry upright (`SlimRailPx`, 12 px at the band's right edge beside the right
  gutter), drawn only while the lanes overflow; its window the lanes in view; pan only (lanes have no zoom: no grips).
- **A cut edge says so**: a 1 px `msui-border` line at the band's top or bottom where lanes are cut (the scroll's
  affordance at rest, quieter than the rail).
- **Drag past the edge** (EC-150): a drag-select leaving the band autoscrolls it by time, as the rows band.
- **History and sessions** (the one rule): a band's scroll is kept as the top lane's first item key and restored by
  it, clamped when the lanes change.

### 4. The first lane labelled, no blank lanes

- **The block of lanes starts at the band's top padding** (`laneGeometry` already does past two lanes of slack,
  `ranges.ts:150–154`); centred only when it is one or two lanes in a taller fixed band.
- **A label that has no room before its model is pinned inside at the data's left edge** (a floating label on its lane,
  over the model's hairline, with the label backdrop), never dropped while its lane is shown. Why: the first model of a
  region starts at its edge; JBrowse floats labels the same way.

## Defaults: what a reader sees

TP53's region from Ensembl's GFF3, 47 transcripts, cap 240 px:

```
 ruler   |7,661,000         7,668,000         7,675,000         7,687,000|
 Genes   | WRAP53 ━━━━━━━━━━━━━━━━━━━            TP53 ━━━━━━━━━━━━━━━━━━━━|
         ─────────────────────────────────────────────────────────────────  (cut edge: none at top)
 Transcr.|TP53-201 ▮▮━━▮━━━▮━━▮━━━━━━━━━━━━━━━━━━━━━━━━━━━━▮▮▮          ┃│
 (47)    |TP53-202 ▮▮━━▮━━━▮━━▮━━━━━━━━━━━━━━━━━━━━━━━━▮▮▮              ┃│ <- rail: window 16 of 47 lanes
         |TP53-203  ▮━━▮━━━▮━━▮━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━▮▮         ││
         |   …  (16 lanes in view)                                       ││
         ═════════════════════════════════════════════════════════════════  (cut edge: 1 px, lanes below)
 next track (fixed)
```

A pending source, then settled:

```
 Transcr.| ░░░░░░░ loading ░░░░░░░ |        (one lane, the placeholder)
 Transcr.| 9 lanes, 134 px, no rail |       (settled: content height under the cap; tracks below moved once)
```

## The interaction matrix's cells

| Row | Cell | Now |
|---|---|---|
| a track or layer | pan | **new**: the wheel's vertical delta over an overflowing lanes band scrolls its lanes; drag the slim rail; Page keys with the cursor in it |
| a track or layer | edit | `height` by patch is its cap (as before, a height); `overflow` by patch (`scroll`, `compact`, `density`) |
| an element (an item) | focus | double-click / `address.focus` also reveals its lane in the band |
| the selection | extend | Shift+arrows reveal the lane the cursor lands in |
| a range | select | a drag past the band's edge autoscrolls the lanes (EC-150) |
| a leading section | pan | the wheel's vertical delta reaches the section only once a lanes band under the pointer is at its end (chaining) |

## How it meets the other rules

- **The empty-well rule** (fixes.md): settled and empty, the note row; pending, the one-lane placeholder; the sweep's
  "draws at least one mark besides its chrome" holds for a scrolled band (its lanes in view).
- **Export** (layouts.md Export, EC-125): a figure shows data, so `renderSvg` and PNG draw **every lane** of a capped
  track at its content height (as a text section draws every line), the live scroll untouched; `export.lanes: 'view'`
  draws the band as on screen with its cut edges. Bounded: past `export.maxTrackPx` (SpecBounds, 4,000 px) the track is
  drawn to it and the rest said in the track's header ("+n lanes not drawn"). No rail in an export (chrome).
- **The rows layer**: the rows band takes the same rule when it has no `height`: its content (rows × row height) to
  a cap (`rows.maxPx`, 480: an alignment wants more than features), then its existing scroll and rail. `height: fill`
  in a fixed-height view stays as built. The canvas and the SVG read the same number (the twin sweep, fixes.md).
- **Stacked blocks** (`arrange: 'stack'`): each block's band of a track is sized by that block's lanes, capped and
  scrolled on its own (one `bandTop` a block band); in `row` one band across the blocks.
- **Text** (EC-96): lanes per line compacted, no inner scroll; the text section scrolls as a whole (unchanged).
- **A ring** (EC-108): **no inner scroll**. A ring's track is radial; a scroll inside an annulus has no direction a
  reader can read, and growing it eats the hole. On a ring a lanes track takes `min(content, rings.trackMaxPx)`
  (the radial width, `RenderConfig`, 48 px) and then compacts (`overflow: compact`, forced), then the density row; the
  ring's linear detail (its follower) is where the lanes are read whole. The track's header (on the linear view) and
  the ring's hover card name the count.

## Edge cases

- **EC-169** pending never drawn as empty: the one-lane placeholder (§1).
- **EC-193** lane heights are the viewer's; a track height in px is asked for: `track.maxPx` and a track's `height`.
- **EC-226** lane order: lanes stay sticky across scroll and growth (an item never changes lane on a pan).
- **EC-150** autoscroll past an edge, by time: §3.
- **EC-155** keys: Page keys only with a cursor shown.
- **EC-125** export never changes the live viewport: §Export.
- **EC-108**, **EC-96**: rings and text, above.
- **EC-137** a size channel shrinks a mark, never its lane: unchanged under scroll.
- **EC-new-1** a tile adds lanes while the band is scrolled: lanes appended below, `bandTop` kept, the rail rescales;
  the item under the pointer does not move.
- **EC-new-2** a filter (`layer.filter`) leaves fewer lanes: the height shrinks (data identity changed), `bandTop`
  clamped to the new content.
- **EC-new-3** a linked hover or a selection on an item in a lane scrolled out: a tick on the slim rail at its lane
  (as the navigator's marks); no scroll on a hover, a scroll on focus.
- **EC-new-4** two tracks of one kind, one overflowing: one height, one rail; their lanes never scroll together.
- **EC-new-5** a track's `height` smaller than one lane: one lane drawn, the rest scrolled (never a band under a lane).

## Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `track.maxPx` (new) | `RenderConfig` preset | 240 | a cap independent of the view's height; matches `text.maxHeight`'s kind | px; `none` (content to the hard bound) |
| `track.overflow` (new) | spec param (track), default in `RenderConfig` | `scroll` | read every lane at full height | `compact`, `density` |
| `maxLanes` | `RenderConfig` | 5 → retired as the content cap | replaced by `track.maxPx` | — |
| `ranges.maxLanes` (new) | `SpecBounds` | 1,000 | layout work per frame is bounded | any |
| `rows.maxPx` (new) | `RenderConfig` preset | 480 | an alignment's rows want room | px |
| `rings.trackMaxPx` (new) | `RenderConfig` preset | 48 | a ring's radial width cannot scroll | px |
| `export.lanes` (new) | spec param of the view; `renderSvg` option | `all` | a figure shows data | `view` |
| `export.maxTrackPx` (new) | `SpecBounds` | 4,000 | an export's size is bounded | px |
| the rail's width, the cut edge | fixed | 12, 1 px | the rows rail's geometry | — |

## Tensions

- **A track grows once its tiles settle**, moving what is below it once per settle. Taken over a guessed height (blank
  or cut). A host that cannot accept any move sets a `height`.
- **Two scrolls in one view** (the page and a band) can trap a wheel. Chaining at the ends and no wheel when nothing
  overflows (the text box's rule) are the answer; a trackpad's momentum may still stop at a band's end once.
- **The author's five-lane cap of 2026-11-19 is revisited**: compaction becomes opt-in, scroll the default.

## Self-check

(a) Inverse: content shrinking (a filter) shrinks the track, a pan never does. (b) Symmetric: the rows band and the
lanes band, one rule, one rail. (c) Plural: tracks of one kind share a height, scroll apart. (d) Layouts: tracks
(scroll), text (no), ring (no, compact), stacked blocks (per block). (e) Export: every lane, bounded. (f) Pending,
empty, settled: three different pictures. (g) Keyboard and accessibility: the band announces "lanes 1–16 of 47" on
settle as the rows rail announces rows.

## Gaps for the author

1. The cap as px (`track.maxPx`, 240) rather than a share of the view? **Recommend px**: a share is circular in a
   content-sized page and resizes on a window resize.
2. Retire the five-lane content cap, scroll by default and compaction opt-in? **Recommend yes**.
3. An export draws every lane of a capped track (bounded at 4,000 px)? **Recommend yes**; `export.lanes: 'view'` for
   the on-screen crop.
4. No inner scroll on a ring: compact past `rings.trackMaxPx`, the linear detail for the lanes? **Recommend yes**.
5. Float a label with no room before its model at the data's left edge, rather than drop it? **Recommend yes**.
6. The rows band without a `height` takes the same rule with its own cap (480)? **Recommend yes**, with the twin sweep.
