# Layouts: text, blocks and rings

## Intent

The third slice's "done when": "The section layout interface carries four layouts (tracks, text, blocks, rings) and the
dot plot as its own section, with no layout-specific code in the viewer." Its known risk: "The layout interface has
carried one layout (tracks); the text layout and rings are its real test, and a change to it touches the draw lists'
key and the hover index." [ux-foundations.md](ux-foundations.md), section 0, already promises it: "Layouts share one
internal interface: tracks, text, rings and dot plots take a section of the resolved view, the settings and a size,
and give geometry, hit-testing and draw lists".

Today the interface is implicit: `layoutView` stacks sections, each one `AxisTransform` (positions to x) over one data
box, and the viewer does the same arithmetic again (an offset from a pointer's x, a pan by pixels, a zoom about an x,
the section under a y). This design names what a layout must give so that layer types, emphasis, hits, the key, the
viewer's navigation and export read the layout through one interface, then fills it three more ways: **text** (the
p53 text figure, the Globins ESPript figure, EC-94 to EC-100), **blocks** (the Complex's chains, EC-101, EC-102,
EC-115) and **rings** (the Plasmid ring with its linear detail, EC-108, EC-114, a feature across the origin, EC-9 to
EC-14). The dot plot (two axes) is its own section ([grids-and-rows.md](grids-and-rows.md)), and is out of this design except
where the interface must leave room for it.

## Precedents

- **The previous implementation** (read as a source): three geometries, `linear`, `wrapped`, `circular`, declared per
  renderer. Rings were "the same frame unrolled and mapped onto the circle": every renderer drew linearly and a backend
  wrapper mapped rects to annular sectors, points to polar places, cut long lines so they bend; hits came back through
  the inverse (`toUnrolled`). Taken: the unrolled frame and the warping writer, so no layer type knows rings. Not taken:
  wrapped text as a geometry each renderer opts into, with a table of bands per layer kind (underline lanes, strip
  lanes, letters' band) in the layout: that is layout code that knows layer kinds. Its text view's rules (whole groups
  of ten, short last group on its own line, inline labels in line breaking, lanes compacted per line, empty bands
  collapsed, frames with open sides, chips, rows in text) are kept as behavior (EC-94 to EC-100).
- **ESPript**: alignment figures: rows named left of every line, the row's last residue number right of it, group
  numbers of a reference row above, conservation shading, secondary structure above the block, frames around regions.
  Taken: the figure's parts as the text layout's gutters and the ruler and frame layers. Not taken: PostScript-only
  output, a figure with no interaction.
- **Jalview's wrapped view**: an alignment wrapped to the window's width, names and start/end numbers per line,
  annotation rows under each line, scrolled vertically, scrolled to a selection made in a linked structure; font size is
  the zoom. Taken: zoom as the letter size, scroll as lines, numbers per row, the scroll to a linked selection. Not
  taken: residue editing in the view.
- **UniProt's sequence text**: groups of ten, a number per group or per line. Taken: the defaults of the p53 text figure.
- **SnapGene and Benchling**: a circular map, a linear map and a sequence view over one model, one selection across
  them; features as arrows on rings, labels outside with leaders, the name and length in the hole; Benchling's "set
  origin" (a rotated linear view). Taken: the ring with a linear detail, labels outside, the title in the hole, the
  rotated linear view as a range over the origin. Not taken: editing; a map that zooms itself.
- **Circos and pyGenomeViz**: concentric tracks, ideogram arcs with gaps between blocks (chromosomes, contigs), chords
  between positions, ticks outside; Circos scales segments (a zoomed region takes more of the circle). Taken: blocks as
  arcs with gaps on one ring; chords as section-scope layers through the same warp. Not taken: static figures only;
  scaled segments (a tension, below).
- **Angular Plasmid**: a declarative plasmid map (tracks, markers, labels, scales as elements). Taken: a ring described
  as data (rings are tracks, marks are layers). Not taken: markup per mark.
- **This codebase**: `AxisTransform` (pieces as display units, cuts at piece edges: the frame rule), `paintRanges` (the
  emphasis over a box), `followRange` with `arc: 'shortest'` (EC-12), a cyclic follower's wrap (`follow.cyclicArc`),
  `TextMeasure` in `@molstar/core/util/text-measure` (pretext's prepare and layout split, proposed for lines in
  [rendering.md](rendering.md)).

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| Select residues by click and drag across wrapped lines (Jalview, SnapGene's sequence view) | built: a drag reads the position under the pointer in reading order |
| Shift+click to extend across lines; arrows and Shift+arrows, up and down by a line | built: the cursor's steps through the layout |
| Change the font size (Jalview), the line width follows the window | built: zoom is the letter size; a resize reflows |
| The wrapped view scrolls to a selection made in a linked view (Jalview, from a structure) | built: a text section reveals the first line of a selection projected onto it from elsewhere (`reveal.linked`) |
| Click a sequence's name to select the row (Jalview) | planned ([grids-and-rows.md](grids-and-rows.md): a row selection) |
| Copy the selected sequence as text (SnapGene, Benchling, Jalview) | planned: `selection.copyLetters`, what is shown (decided, below) |
| Hide and show annotation rows (Jalview) | built: `track.hide`, the band leaves every line |
| Click a feature on the map, it is selected in every view (SnapGene, Benchling) | built: hit through the ring's inverse; selection is projected as everywhere |
| Click or hover a feature's label outside the ring (SnapGene) | built: a label is hit as its item |
| Drag on the map to select a region, across the origin (SnapGene, Benchling) | built: an arc, two runs over the origin (EC-9) |
| Move the visible region of the linear map on the circular one (Benchling, the previous implementation's window) | built: the ring's window, dragged and resized |
| Set a new origin (Benchling) | built: a section's `range` over the origin, and the ring's `rotation`; by patch, with a command |
| Zoom into the circular map itself (Circos zoomed arcs and scaled segments) | refused: zoom on a ring zooms its window; an arc of the ring is a section `range` (a standing tension) |
| Rotate the map by dragging (none of the precedents makes it a gesture) | refused: a drag on data selects; rotation is a param, by patch |
| Edit the sequence in the text or on the map (SnapGene, Jalview) | refused: data is not edited in the viewer |
| Toggle restriction enzymes, primers (SnapGene) | not here: they are layers (hide a layer) |
| Drag a selection past the edge to scroll (text editors, Jalview) | planned (text here, rows in [grids-and-rows.md](grids-and-rows.md)): clock-timed, faster further past |
| Scroll a text by touch, a scroll grip on its edge (mobile readers) | refused this slice: no touch gesture of its own (a standing tension); the wheel, keys and Space+drag scroll |

## Decision

### The author's decisions

The author, 2026-11-18: "layouts.md: ok as recommended"; in general: "a lot of this is "configurable behavior" and
where appropriate, the code should allow the configuration via spec ports etc."; on sizes: "don't overdo it on the
limits (they were sort of self-imposed and estimate), but obviously try to minimize the amount of code". On the review of the decisions: "ok" to its recommendations.

1. Blocks are a cut every layout takes (`block` children, `arrange` in tracks), not a fourth `layout` value.
2. A layout per section, sections stacked, so a ring and its linear detail share one view; side by side later.
3. A text section exports every line (the figure), not the lines in view.
4. A text section in a host page caps at `text.maxHeight` (600 px) and scrolls inside, taking the wheel only while it
   overflows and giving it back at its ends.
5. A collapsed intron is zero width on every layout, rings included (no `rings.pieceGap`), a hairline marking it.
6. This design owns the viewport record (per-band `bandTop` and `rowScale`, the grid's `y`) and the y path in the
   emphasis passes; a history step over a changed layout is restored by position, never silently dropped.
7. No size cut order: the code is kept minimal and measured.
8. A drag on a ruler, the author: "i think drag on a rule should pan. or rather .. the original solution solved this by
   showing the rail as well ... maybe that's the way to go?" Adopted ([navigators.md](navigators.md)): the ruler rail; a positions
   ruler pans, an alignment's ruler pans in its lower 40% and selects columns above, a grid's axes pan.
9. The author, on the remaining questions: "others as recommended": `selection.copyLetters` copies what is shown
   (reading order, complement as drawn, shown positions, rows as FASTA lines), an option for the stored letters;
   `export.window` is an optional spec param of the view, `renderSvg`'s `window` overriding it; figure-meaning settings
   (`group`, `numbers`, `copies` and the like) are spec params.

### A layout cuts a section into strips and warps them

- **A strip** is a run of the section's display units drawn straight: one `AxisTransform` (its pieces, its window, x
  across a width), a stack of **bands** (one a track, with its height), and a **warp** from the strip's own frame to the
  page. A **layout** is a rule that cuts a section's display units into strips and places each with its warp:

  | Layout | Strips | Warp |
  |---|---|---|
  | tracks | one a block in view (one for a section without blocks) | a translation |
  | text | one a line in view (and a chip a one-position block) | a translation |
  | rings | one a block, unrolled along the outer ring | polar: x to an angle, y to a radius |

  Why: a layer type then draws one way, over a straight axis and a box, in every layout, as it does today; the
  previous implementation's unrolled ring proved the warp for rings, and a wrapped line is the same strip translated.
  A layout knows strips, bands and warps, never layer kinds. Layers declare no geometry kind (`linear`, `lines`,
  `polar`): the warp is the layout's.
- **Layer types keep their contract** ([layout-and-draw-lists.md](layout-and-draw-lists.md)): `layout` is called once
  a strip, with `ctx.box` the band in the strip's frame and `ctx.axis` the strip's transform; `draw` writes into the
  strip's warping writer; `hit` is called with the point brought back into the strip's frame.
- **The lay cache is bounded by the strips in view**: a layer's last layout is kept per strip, by the strip's first
  position (a line, a block), for the strips in view only; a strip that leaves the view lets its lays go. A scroll keeps
  the lays and sticky lanes of the lines still in view (lanes are per line in text, EC-96); a reflow changes every
  line's first position, so its lays are dropped and each line's lanes packed fresh (lanes stay sticky across a scroll
  or an arriving region, not across a reflow).
- **Additions to the layer contract**, each generic:
  - `used?(lay): number`: the height the layer used in its band on this strip, 0 when it drew nothing there. Text
    collapses a band no layer of it used (EC-96); a `values` layer answers its height (a gap in its data is data). Why
    a number from the layer: only the layer knows its lanes, and the layout must not.
  - `scope: 'strip' | 'section'` (`strip` by default): a `section` layer is laid out once a section with the placement
    as a reader (`ctx.place`, below) and draws in page coordinates: what crosses strips (a pair's arc from one line to
    another, a chord through a ring's hole, the pairs lane). Why: a mark with two ends on two strips belongs to neither.
  - `ctx.labels`: a label sink (`ctx.label(text, start, end, element)`) a layout collects into: on rings and linear
    tracks, every labeled item's label with its element key, inside its mark when it fits, else beside (below).
  - `ctx.span`: the union of named tracks' bands on the strip, for a layer with `across` (frames).
- **Frames across tracks** (EC-100) are a layer type of strip scope with `across` (two tracks of its section): its
  outline covers `ctx.span`; an end cut by the strip's edge (a line break) is open, from the axis's cut flags, as every
  cut end is (the frame rule); its title is drawn where the item really starts.
- **A track's height by default fits its content, in steps** (decided, the author, 2026-11-19, Antibody's CDR tracks: a
  72 px band holding one lane of marks at its top "feels empty"; "combination of both … smart layout so the heights
  don't vary too much"): a track with no `height` is as tall as the lanes its layers need (ranges: the packed lanes of
  the whole data, not of the view, so a pan never resizes it), rounded up to a step of one lane, between one lane and a
  cap (then lanes compact, as now); its marks centred in the band. Tracks of one kind in a section share the tallest
  height among them (two CDR tracks are equal), so neighbours do not vary lane by lane. A track's `height` overrides.
  Why: a default fills its space and reads as designed; heights decided by the whole data stay put as the view moves.
  The numbers (decided, the author, 2026-11-19, the defaults audit `drafts/defaults-audit.md`, "Accept as drafted"):
  a lane 14 px (a 12 px mark, a 2 px gap), 4 px padding top and bottom, so one lane 22 px (the ruler's band), two 36,
  at most five (78 px), then lanes compact; an empty track the one-lane minimum; a label lane counts as a lane when a
  label goes beside its mark. Values, letters, rulers and pairs keep their fixed heights (they have no lanes). **The
  algorithm is the spec's** (the author, 2026-11-19: "configurable in the spec … something like 'layout algorithm'"):
  the section's `heights` param names the sizing, `content` (the default, above), `uniform` (every track the tallest
  content height in the section), `fixed` (each track its layer type's fixed default, the second slice's look); an
  extension may register more (`ext:` names, as layouts are registry entries); a track's `height` overrides any. **A
  host overrides it by a port** (the author, 2026-11-19: "configurable/overrideable via a port too?"): the composition's
  `trackSizing` port (a factory, as every port, [viewer-composition.md](viewer-composition.md)), given each section's
  tracks with their content heights and the named sizing's answer, returns the heights; absent, the named sizing
  stands. Why: an application fitting the viewer into its own panel (a fixed total, a shared row grid) decides
  heights from what only it knows, without a registry entry per app.
- **Defaults decided from the audit** (the author, 2026-11-19, `drafts/defaults-audit.md`): the label column is
  measured to the widest track label plus its gap, between 80 and 160 px, an ellipsis past that (never a hard clip);
  sections 20 px apart (`sectionGap`, was 8), tracks 4; the text layout's `lineGap` 10 (was 8) and its `bandGap` 3 (was
  2). On rings, more than six labels within 12° are one grouped label ("13 sites"), opened on hover or zoom.
- **A band may scroll inside its section** (asked by [grids-and-rows.md](grids-and-rows.md)): a band has its own vertical scroll
  (`bandTop`, clamped to its content's height) and row zoom (`rowScale`, px a row), both in the section's viewport by
  band (below), beside the section's `top`. The band's layer reads them (`ctx.bandTop`, `ctx.rowScale`) and owns its y
  inside the band (the rows layer places its rows), and they are part of the section's key. So an alignment's rows
  scroll while the ruler and logo above and a summary below stay fixed in the same section. The wheel over such a band
  scrolls it first (then, at its ends, the section, then the page).
- **Glyph runs gain a per-letter y, angle and vertical scale** (optional arrays beside the per-letter x): one change of
  the draw lists' public form (versioned, in the changelog), shared by rings (letters along the ring), the rows layer's
  letters and the logo (letters stacked by height). Why one: three uses need letters off one baseline at one size, and
  three forms would mean three backends' code.

### The interface

Internal, in `seqx-core/render` (`layout/`), until an application needs one ([extensibility.md](extensibility.md)).

```ts
interface SectionLayoutType<Kept> {
  readonly name: string;                                   // 'tracks', 'text', 'rings'
  /** Cuts the section into strips for its box and viewport, into prev's buffers. */
  place(section: ResolvedSection, ctx: PlaceContext, viewport: SectionViewport, prev?: Kept): Placement & Kept;
  /** Whether two viewports place alike: the part of the section's key this layout reads. */
  sameKey(a: SectionViewport, b: SectionViewport): boolean;
  /** Gestures to viewports, pure: what the viewer calls instead of doing arithmetic. */
  readonly nav: SectionNav;
  /** The layout's own marks a pass: label column, gutters, numbers, the ring's guides, labels and leaders, a window. */
  chrome(placement: Placement, pass: PassName, out: PassWriter, emphasis?: Emphasis): void;
}

interface Placement {
  readonly box: Box;                    // the section's region on the page
  readonly contentHeight: number;       // all of it (text: every line), for scrolling and export
  readonly strips: Strips;              // struct of arrays: transform, system, block, bands' tops and heights, warp, cuts
  /** What is under a page point, into a reused record: strip, band, the point in the strip's frame, kind. */
  at(x: number, y: number, out: PlacePoint): boolean;
  /** The offset under a point, nearest when between strips or in a gutter (drags, zoom anchors). */
  offsetAt(x: number, y: number): number;
  /** Where an offset is drawn: its strip and page point (section-scope layers, reveal, the cursor). */
  place(offset: number, out: PlacePoint): boolean;
}
```

- `PlaceContext` holds what `LayoutContext` holds of the view (the size, the device pixel ratio, the config, the
  tokens, the text port, the model, the store, the overlay) and the section's box. A layout never reads layer
  parameters.
- **One frame, a viewport per section** (`FrameInput.viewports`, by section address; owned here): replaces the one
  `viewport`, and settles [layout-and-draw-lists.md](layout-and-draw-lists.md)'s open item "one window per leading
  section". The viewport carries `top`: vertical scroll is a new input of core's layout and of its key (core has none
  today; the viewer scrolled by itself).
- **Follow reads the leader's viewport, not its transform**: `followAxis` today carries `leader.visible()`, which on a
  ring is the whole axis; it reads the leader's viewport from `FrameInput.viewports` (its window on a ring, its window
  in tracks, the lines' range in text) and carries that. Tested both ways: a linear detail following a ring's window,
  and a ring following a linear section.
- **The view stacks its sections** (each placement's box under the one before), whatever their layouts; the layout is
  a section's (`layout` on the section, the view's by default), so the Plasmid view stacks a ring and its linear detail.

### Placement and hit-testing as its inverse

- `at(x, y)` is the inverse of the warp then of the cut: the page point to a strip (tracks: the block's box; text:
  the line by binary search over the line tops, a typed array; rings: `atan2` and the radius), to a band by its top,
  to the strip's frame. `hitView` becomes: the section under the point (`at`), the band's layers topmost first, each
  `hit` with the strip-frame point; a hit is the same record as today, plus `block` (its index in the section) when the
  section has blocks; `system` is the strip's. A layer's hit box is in its strip's frame, so every drawn mark keeps its
  hit box through the warp (EC-123).
- **The kinds of point** `at` gives: `data` (a band), `label` (a ring's outside label), `gutter` (text's row labels and
  numbers), `window` (a ring's window, with its part), `gap` (between groups or lines: no position for the hover; the
  nearest for a drag), `empty`. They are the placement's, not new hit targets: a label hits as `item` (by the element
  key, extent and hull the label sink was given, never a second hit at an anchor, EC-108); a gutter as `row` (the rows'
  selection, [grids-and-rows.md](grids-and-rows.md)); a window as the one `window` target with its `part` (`body`, `start`, `end`), shared with the navigators.
  Why: a target names what a binding acts on, and these act on an item, a row and a window.
- **Hit-testing is checked as the inverse of placement** by property tests over random sections in each layout: every
  position in view is placed once (`place`), and `at` of its cell's center gives it back; `offsetAt` is monotone along
  the reading order (text) and the angle (rings).

### Emphasis per strip

- The view's emphasis code (`drawSectionEmphasis`, `paintRanges`) takes a strip, not a section: the selection's tint
  and outline, the highlights, the projected hover and the cursor are painted over the strip's **band column** (all its
  bands, `emphasis.selectionShape: 'column'`) or each band (`tracks`), in the strip's frame, through its warp, on the
  strip's system. So the selection is one look in every layout ([ux-foundations.md](ux-foundations.md), section 6): a
  band on a line in text, an annular sector over every ring on a ring.
- **Cut ends follow the strip**: a selection or an item cut by a line break is square at the break and squircle where
  it really ends (EC-100); on a ring the seam at the origin is no cut: an item over the origin is two pieces whose
  square ends meet, one arc (EC-9).
- **The y path** (asked by [grids-and-rows.md](grids-and-rows.md)): a layer with a second axis (`ctx.axisY` on a grid strip, the rows
  layer's own rows) declares `axes: 2`; the emphasis passes then also call its `draw` with the emphasis records' parts on
  its y system (a row's band, a block's rect, the crosshair's rows), which it paints through its own y transform. The x
  painting stays the strip's. Why the layer's: only it knows its y transform; the strip's painting knows x.
- The emphasis passes walk the strips in view with indexed loops through reused writers: nothing allocated warm, as
  today.

### The key, the caches and the hover index

- **The key per section**: the renderer keeps, per section, its placement, its lays and its key (its viewport through
  its layout's `sameKey`: tracks start, end, top and each band's `bandTop` and `rowScale`; text the anchor, the scale, top; rings rotation and size only, never
  the window; and the section's reads); the view's parts (size, device pixel ratio, tokens, text, overlay, config)
  stay global. A frame lays out only the sections whose key changed; an unchanged section's kept placement and lays
  are drawn again into the frame's lists. The lists stay one a pass for the frame (no public change). Why not lists a
  section: a backend and every reader of the lists would change for a cost (drawing kept lays) that is a draw, not a
  layout.
- **The ring under a moving window**: a pan of the linear detail changes the ring's viewport only (its window), so the
  ring is not laid out again and its labels are not placed again; when the frame rebuilds for the detail, the ring's
  base passes are drawn from its kept lays. The ring's frame target (8 ms, below) covers that draw; when nothing but
  windows and emphasis changed (a window dragged with no follower relaid yet, a hover), only the emphasis passes are
  drawn.
- **The `window` pass** (`render/window.ts`, here; [navigators.md](navigators.md) extends it): first of the emphasis passes,
  under the highlights, on the overlay canvas, drawn by a layout's chrome from the section's viewport: a ring's window
  now, the navigators' windows next. A change of the draw lists' public form (versioned, in the changelog). Why a pass,
  not the highlight pass: a window is not data emphasis, and a highlight channel must stay above it. **The records**:
  the recorder writes an empty pass as its header, so the commit adding the pass regenerates every committed record,
  its only diff the `pass window` line (checked by a diff in that commit's review); the tracks restatement then keeps
  them byte-identical.
- **The hover index is the placement's `at` then the layer's `hit`**: O(log lines) in text, O(1) plus O(log rings) on
  rings, O(blocks) in tracks; the hover's 2 ms target and its 0 bytes warm hold in every layout (a quick representative
  per layout, below).

### Navigation: one viewport record, the layout's arithmetic

- **The viewport stays one record**, owned here: `{ start, end, top, width, pending, scale?, bands?, y? }` in offsets,
  as today plus `bands` (by band: `bandTop` and `rowScale`, a band that scrolls its own rows) and `y` (`{ start, end }`,
  a grid's second axis when it has its own window). Every design names these fields as here. Each layout reads it its
  way:

  | Layout | `start`, `end` | `top` | `scale` | `bands`, `y` |
  |---|---|---|---|---|
  | tracks | the window of display units | vertical scroll | derived (`width / (end − start)`) | `bands` |
  | text | `start` the anchor (the position at the top line, fractional within its line); `end` derived: the last position the lines in view hold at full band heights | derived from the anchor | px a position: the letter size | none (a text shows every row) |
  | rings | the window (what a linear detail follows) | none | unused (the ring shows its whole axis) | none |
  | grid ([grids-and-rows.md](grids-and-rows.md)) | x's window | none | derived | `y` unless y follows x |

- **`SectionNav`** gives every move the viewer makes, pure: `pan(v, dx, dy)`, `zoom(v, factor, at)`, `fit(v)`,
  `reveal(v, start, end)` (a focus, a link, the cursor), `step(position, dir)` (the cursor: left, right, up, down,
  home, end), `clamp(v)`. The `ViewModel` keeps groups, follows, history, intents and clamps as today and calls the
  section's nav for the arithmetic: `clampRange` becomes the tracks nav's `clamp`; text clamps its anchor and size,
  rings their window (wrapping, never clamped at the total). The input asks the placement for the section and the
  offset under a point (`at`, `offsetAt`), never `dataX` or `v.width`. Why: today's offset-from-x and pan-by-px code in
  `ViewModel`, `input`, `select` and `navigate` is the tracks layout's, and would be wrong in text and on rings.
- **Saved state**: a section's viewport is saved as positions plus its `scale` (text) and each band's top row key and
  `rowScale`, never offsets (the `AddressModel` reads and sets the session as addresses). The view history keeps whole
  viewport records with their layout and axis signature; **a step over a changed axis or layout** (one rule, in every
  design that touches history) is restored by position (the centre and the zoom; text its anchor and scale; a band its
  top row), clamped to what is drawn, never silently dropped. The language moves to `0.4.0` with no migration
  (pre-release, core `_spec` Versions).
- **A drag on a ruler** ([navigators.md](navigators.md), the ruler rail): on a positions axis the ruler row is the rail's grab zone,
  so a drag there pans (ux-foundations section 8: a drag on an axis pans); on an alignment's columns only the ruler's
  lower 40% pans (`navigator.rail.grabShare`), and a drag above it selects columns; a grid's axes pan. In text and on
  rings the ruler is drawn in the strips and takes no rail (a navigator there is declared).
- **Copying the selection** (`selection.copyLetters`, decided): copies what is shown: reading order, letters
  complemented as drawn, only the positions shown (not a collapsed intron's), rows as FASTA lines (`>label` then the
  row's letters); its option `stored` copies the stored letters of the selection's span instead.
- **Text with rows** (the ESPript figure): ↑ ↓ step the rows of a line block and, past its last row, cross to the
  next block's first row; a drag selects the rows between its ends × the columns between them in reading order
  ([grids-and-rows.md](grids-and-rows.md), the rows' selection).
- **A drag past a scrolling box's edge** (a text section, a rows band) scrolls it toward the pointer, faster the
  further past, timed by the clock (the animation port's frame time), never by the frame count; the drag's selection
  extends as content moves under the held pointer. Planned (text here, rows in [grids-and-rows.md](grids-and-rows.md)).
- **Text**: the wheel scrolls lines (pan is vertical); Cmd/Ctrl+wheel, pinch, `-`, `=` and Shift+double-click step the
  letter size along the type scale, the position under the pointer kept at its y; `0` returns to the preset's size; a
  resize reflows with the anchor's line kept at the top (`viewport.onResize` is the anchor in text); up and down move
  the cursor a line. The size never goes below the letters' threshold (a text shows letters, EC-116) nor above
  `viewport.maxPxPerPosition`.
- **The text box in a host page**: a text section's height is its `height` param, else its content's up to
  `text.maxHeight` (600 px), past which it scrolls inside its box. The wheel scrolls its lines only while its content
  overflows the box and the pointer is over it, and gives the wheel back to the page at its top and bottom ends (scroll
  chaining), so a host page keeps its scroll; a section that shows all its lines takes no wheel. Keys scroll it while the
  viewer has focus.
- **Rings**: the wheel and Shift+Alt+Left/Right move the window round the ring (the cursor's Left/Right go round too,
  wrapping at the origin, EC-154); Cmd/Ctrl+wheel, pinch, `-`, `=` and Shift+double-click narrow and widen it about the
  pointer's angle (the window's centre with keys); `0` fits it to the whole; double-click focuses an element: the window
  takes its shortest arc (EC-12). A ring nobody follows has no window: its zoom and pan are idle, as a navigator's
  (EC-110).
- **A cyclic section's viewport wraps**, leader or follower, in every layout: its end may pass the total, its pieces
  split at the origin (`follow.cyclicArc: wrap` made general). It gives the rotated linear view (a plasmid panned
  across its origin) and the ring's window over the origin.
- **Blocks** (stacked): a zoom is about the block under the pointer (EC-115); the cursor crosses block gaps (EC-154).
- **The story `view` channel** ([stories.md](stories.md)), field by field of the viewport record: tracks and rings
  `start`, `end` as positions (a ring's window round its shortest arc); text its anchor as a position and `scale`
  (the anchor's line kept at the top while the size tweens); `top`; each band's `bandTop` as its top row's key and
  `rowScale` (interpolated in display rows); a grid's `y` as positions. A section's `rotation` (a spec param) tweens
  between steps with the ring's outside labels hidden, placed once at its end (labels are placed only when the ring is
  laid out; a tween must not place them every frame).

### Follow between a ring and a linear section

- **The ring's window is its viewport**, so a linear section follows a ring with today's follow (`follow: ring`, a
  path of identity or any map), reading the ring's viewport (above): its range is the ring's window carried along the
  path. Dragging the window pans the ring's viewport; a grip resizes it with the other end held; a zoom or pan asked of
  the linear follower moves the window through the follow (`follow.followerMoves: lead`). That is EC-114's "a ring
  leads only through its window", with no follow code of its own.
- **The symmetric case holds too**: a ring that follows a linear section (`follow` on the ring) has the leader's range
  as its viewport, the same window. Either way the picture and the gestures are the same; which leads is the spec's.
- **The window's inverse**: a change of the ring's `rotation` moves where the window is drawn (the ring's key changes,
  it is laid out again) but not the window's offsets, so the follower does not move.
- A window over the origin wraps (two pieces on the linear detail, `follow.cyclicArc: wrap`); a focus on a feature
  across the origin takes its shortest arc on the ring and its two pieces on the detail (EC-9 to EC-12).
- The window is drawn as the navigators' window (a sector over the rings in `seq-window`, a grip at each end on the rim
  in `seq-grip`, idle when nothing follows): one window mark and its geometry in `render/window.ts`, its hits the
  `window` target with a `part` (`body`, `start`, `end`), its drags one `window` behavior (a `BehaviorProvider` in the
  viewer's spec), built here and extended by [navigators.md](navigators.md). The behavior gives every window (the ring's
  included) a slider with `aria-valuetext` in its section's labels, set when the viewport settles, never per frame;
  Esc during a window drag cancels it, the view back where the press found it. Shared by every window (the rules and
  their reasons in [navigators.md](navigators.md), "The ruler rail"): grips as tabs just outside the ends (resize about the other end),
  lit per part; a floor of 16 px with an inward end drag on a floored window keeping the view; a hit of at least 12 px
  (24 px by touch); idle drawn quiet and not a target; drawn grown with `grab`/`grabbing` while dragged or hovered;
  Shift+arrows resize the focused window.

### The tracks layout, restated

- One strip a block in view, translated; the label column, the section's band, the clip a track as today. Its draw
  lists are byte-identical before and after the restatement (after the `window` pass's header): the committed records
  are its guard. Why: the interface is proven on the layout that exists before the new ones lean on it.

### Blocks

- **A block** is a piece of a section on a system of its own: `block` children of a section (`of`, `range`, `label`),
  or a section `of` a complex: a block a unit, polymers first, then ligands (EC-101), when the complex kind exists. A
  section without blocks is one block. Each strip carries its block's system (`AxisTransform.system` is the strip's,
  not the section's), so a hit (with its `block`), a selection's part and a layer's path are per block: the resolved
  view gives a layer a path a block.
- **The scroll group's key** of a section with blocks is its blocks' systems, lengths and pieces in order: only equal
  block lists group.
- **Follow and blocks, this slice**: a section with blocks on more than one system neither follows nor is followed
  (`seq.render.follow.blocks`, the reason in the issue: which block's map a follow would carry is undecided); blocks on
  one system (pieces of one sequence, copies merged) follow as one axis. A standing tension, below.
- **Tracks**: `arrange: 'row'` (side by side along one axis, one viewport over the blocks laid end to end with a gap
  between; up to 24 blocks share the gaps, which shrink when zoomed out so positions keep three quarters of the width,
  EC-102), or `'stack'` (a block a row of tracks, aligned at their first position, one zoom, EC-102, EC-115).
- **Text**: every block starts a line, its label inline on its first line only (EC-95); one-position blocks (a ligand)
  are chips that flow like words after the chains (EC-99), each a strip one position wide, its label's width.
- **Rings**: a block an arc, a gap between blocks (`rings.blockGap`); a whole cyclic block takes the full circle with
  no gap.
- **A one-position block is a chip in every layout** (decided, the author, 2026-11-19: "chips is a prolly good way to
  do these"): a ligand, an ion, a glycan or a drug standing as one position takes a fixed slot in screen px, its
  label's width (the text chip's rule: its label and `×N`, 5 px either side), never a position's share: zoom does not
  scale it, so it reads at the fit view and is never clipped. Tracks: the chips follow the polymers' blocks (`row`), or
  share one short row of their own (`stack`); each track holding something at the position draws one pill there (the
  letters' band the label, a feature band its mark), the others leave it blank. Rings: a fixed-angle arc a chip after
  the polymers' arcs. The row's x is not linear over a chip, as over a gap (EC-102): the window, the ruler and
  `windowExtent` skip chips as they skip gaps. Hit, hover, selection and links read a chip as its one position. Why:
  a one-position share is a sliver at the fit view (the Complex's ATP), and the base's floor (a 2 px mark, unlabeled)
  says nothing; one rule for text, tracks and rings. **Where a ligand binds is the host's**: its contacts are an
  optional `pairs` layer a builder or an application adds (the Complex and Globins examples add it), never drawn by
  default (the author: "contacts should be optional and configurable by the app ... we are building a library, not a
  one-off viewer"). Two genomes on one ring wait for the Genomes milestone (EC-210).
- **Copies** (EC-65): `copies: 'merge'` (the default) draws one block for a copy group, labeled `×N`, hovered for
  every member; `'each'` a block a unit. Why `merge`: copies hold one molecule's sequence, so one block says it once and
  its room and drawing do not grow with the copy count (EC-67); `each` is the writer's when copies differ in what is
  shown. Toggled by patch (one undo step).
- A block's range may wrap a cyclic origin; a unit's may not (EC-13).

### Text

- **Lines of whole groups** (EC-94): positions in groups of `group` (10) with `text.groupGap` between groups; a line
  holds as many whole groups as fit the width (or `max_width`); a short last group goes on a line of its own where a
  full one would not fit. Breaking is arithmetic on cell widths plus measured inline labels: inline group labels take
  part in breaking, continuation lines have none (EC-95). Lines are laid out lazily, only those in view, their tops a
  prefix sum kept for scrolling ([performance.md](performance.md), text layout).
- **Groups are a property of the transform** (`AxisTransform` `groups: { size, gap }`): a cell skips the gaps, a span
  bridges them, so an underline or a selection runs across the space between groups and is cut only by a line break.
  Why: an item split at every group reads as many items; the gap is spacing, not a boundary (UniProt's and the previous
  implementation's text).
- **Groups are bridged**, fixed: an item or a selection is cut by a line break only. Why: no second use for cutting at
  every group.
- **Pieces in text** (the axis design's reversed pieces and collapsed introns): groups count display units along the
  line, across piece edges, so lines stay aligned; a piece edge does not start a line (`text.pieceBreak: 'none'`; `line`
  starts one) and is drawn inside a line as a boundary (a hairline with its caret, EC-127); the numbers per row print
  each line's labels, so a jump in numbering shows there. A reversed section reads in its display order, its letters as
  the axis design says.
- **Bands under the letters**: a track a band on every line, in track order, its height its own (a text preset sets
  thin ones: a `ranges` track 4 px reads as underlines, a `values` track 3 px as a strip); empty bands collapse per line
  (`text.collapseEmpty`, EC-96). The letters' band is the letters layer's track. No band kinds by layer type: a text
  figure's look is its tracks' heights and the preset.
- **Numbers per row** (`numbers`: `start` by default, `end`, `both`, `none`): the first or last position on the line
  in the band's numbering, printed in a gutter (start: left, after the label column; end: right); `start` because a
  line is read from its start, so its first number places what follows without counting; group labels
  (`group_labels`: `above`, `inline`, `none`) number each group (EC-94). A line with no position of a row has no
  number for it. A ruler band's line numbers sit on the ruler's label baseline (decided, the author, 2026-11-19): one
  row of numbers across the line, start, ticks and end. A position carries one label: a group's end number falls on a
  ruler tick of the same numbering, the ruler draws the tick without its label.
- **Alignment figures** (EC-97, the ESPript figure): a rows layer ([grids-and-rows.md](grids-and-rows.md)) is a band of text rows on each line; each
  row's label in the label column of every line, its last residue on the line in its own numbering in the end gutter;
  group labels number the reference row. The label column fades a long label (`@molstar/ui`'s rule), where the
  previous implementation cut it with an ellipsis: changed. Shading is a theme ([themes.md](themes.md)), secondary structure a track
  above (its small layer, [themes.md](themes.md)), frames the frame layer.
- **Centered and capped lines** (`align`, `max_width`) move margins, hits and the selection with them (EC-98).
- **Reflow** on resize is the layout's, keyed by the width; the anchor holds the reader's place.
- **A linked selection is revealed** (`reveal.linked`, `text`; one `reveal` group in `ViewerConfig` with
  `reveal.cursor`, today's `cursor.reveal`, and `reveal.opensClades` from [grids-and-rows.md](grids-and-rows.md)): when a selection is made
  in another section or a structure and projects onto a text section with none of it in view, the text scrolls its
  first line into view, as a move by hand (no history step). Why text only by default: a text hides most of its content,
  a tracks section shows the selection's projection where it is.
- **A bound**: a text section lays out at most `text.maxLines` lines (100,000); more is `seq.render.text.lines`, asking
  for a `range`. Why: a chromosome as text is no figure, and the line tops are kept for the whole section.

### Rings

- **Rings outside in** (EC-108): the first track the outermost ring, each ring its track's height, a hole of at least
  `rings.hole` of the radius; the diameter fits the section's box (its `height`, else the width, capped at
  `rings.maxDiameter`); the origin at `rotation` (degrees clockwise from the top, 0), positions clockwise; a reversed
  section (the axis design's direction) runs counterclockwise.
- **Pieces on a ring**: a collapsed piece edge (an intron collapsed) is zero width, as on every layout, its hairline
  marking it (decided: no ring gap); the seam of a whole cyclic piece is no break.
- **The unrolled frame** is the outer circumference wide, so the density steps and a layer's reduction per device
  column become a reduction per angular column at the outer ring (EC-108); inner rings are finer than their pixels,
  never coarser.
- **The polar warp** (`WarpWriter`, allocation-free, buffers reused): a rect becomes an annular sector path (a
  squircle's corners warped with it), a path's lines are cut every few device px so they bend, a glyph run gets its
  per-letter places and angles (`rings.letters`: `upright`, or `tangent`), an image is refused on a ring
  (`seq.render.ring.image`: dense matrices on rings wait for the Genomes milestone), a clip rect is applied in the
  strip's frame before warping. Why clip before: a clip by path would change the draw lists' form for one case.
- **Labels inside if they fit, else beside** (EC-108; decided by the author, 2026-11-19, the plasmid's ring and linear
  views reading differently): one rule for rings and linear tracks. A labeled item's text is drawn inside its mark
  when it fits there (the mark's length less a padding, at the track's text size); else it goes beside: on a ring
  outside, round the ring with a leader out along its angle across the ruler; on a linear track in a label lane just
  above the track at the item's start, with a short leader when it had to move. The layout gives layers a label sink
  (`ctx.labels`): every labeled item's text, extent and element key; placed when the ring or the strip is laid out
  (never on a window move or a hover), sorted by angle or position, spread a line apart, upright; those a side or lane
  cannot hold are dropped by the same priority in both (the layer's order, then the longer item). The same text style
  (size, weight) in both; the ink is the contrast colour against the mark's fill inside, the text colour beside
  (decided, the author, 2026-11-19: "contrast inside, text beside"). A placed label's hit is its item's. Why: one feature read two ways
  in two views of one map is two rules to learn; inside first keeps the label on its mark where there is room.
- **The hole holds the title**: the section's label and the system's length (`rings.center: 'title'`).
- **A feature across the origin** is drawn as one arc (two pieces meeting at the seam), selected as two runs (EC-9),
  focused by its shortest arc (EC-12).

### Export

- `renderSvg` and PNG draw what the layouts place: a text section all its lines (the figure: `contentHeight`), never
  only those in view, unless the export's viewport says so; a ring whole, with its labels, the window only with
  the view's `export.window` (an optional spec param, false: an exported map is a figure), `renderSvg`'s `window`
  option overriding it. The size's height may be `auto` (each section its
  content's height). Data is wanted for every line drawn, waited for as today.

### Language and the resolved view

- `section.layout`: `tracks`, `text`, `rings` (the view's `layout` by default; the view's enum grows the same).
- Text params on the section: `group` (10), `numbers` (`start`), `group_labels` (`above`), `align` (`start`),
  `max_width`, `height`. Ring params: `rotation` (0), `labels` (`outside`; `none`). Blocks: `block` children (`of`,
  `range`, `label`), `arrange` (`row`), `copies` (`merge`). `LanguageVersion` moves `0.3.0` to `0.4.0`.
- The resolved view gives each section its layout and its params, its blocks (system, pieces, label) and each layer a
  path a block. A change of the public resolved view (versioned, in the changelog).
- What the figure means is the spec's (group size, numbers, rotation); spacing is configuration, in the presets
  ([ux-foundations.md](ux-foundations.md), section 0).

### Size

The author: "don't overdo it on the limits (they were sort of self-imposed and estimate), but obviously try to minimize
the amount of code". No size cuts a feature; the code stays minimal (one strip contract, one window module shared with
the navigators, no private copies). Estimate (minified), to be measured: the interface and the tracks
restatement 2 kB, text 7 kB, the polar warp 4 kB, ring labels and leaders 3 kB, the window 2 kB, blocks 2 kB: about
20 kB.

## The interaction matrix's cells

The rows and columns of [ux-foundations.md](ux-foundations.md), section 8. "As tracks" means the default binding
reaches the same command, the layout giving the arithmetic.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a text section** | — | — | — | — | — | Cmd/Ctrl+wheel, pinch, -, =, Shift+double-click: the letter size about the pointer's position; 0 the preset's | wheel (while it overflows), Space+drag (no Space toggle on its key-up after a drag): scroll lines; Shift+Alt+up/down a line | no: hide its tracks | `group`, `numbers`, `max_width` by patch (`section.set`, previewed) |
| **a ring section** | — | — | — | — | — | Cmd/Ctrl+wheel, pinch, -, =, Shift+double-click: the window about the pointer's angle; 0 the whole; idle with no follower | wheel, Space+drag, Shift+Alt+Left/Right: the window round, across the origin | no: hide its rings | `rotation` by patch (`section.rotate`, previewed); `range` sets a new origin |
| **a follower** (a ring's linear detail) | as a leading one | as a leading one | as a leading one | as a leading one | as a leading one | moves the ring's window (`lead`) | moves the ring's window (`lead`) | as a leading one | as a leading one |
| **an element** in text | click | Shift+click across lines | Cmd/Ctrl+click | double-click: its first line revealed, the size kept | double-click on it when focused | — | — | no: hide a track | no: data |
| **an element** on a ring, or its label | click (the label as its item) | Shift+click | Cmd/Ctrl+click | double-click: the window takes its shortest arc (EC-12) | double-click on it when focused | — | — | no: hide a track | no: data |
| **a range** in text | drag in reading order across lines | Shift+drag | Cmd/Ctrl+drag adds, Alt+drag removes | Enter at the cursor | Alt+Left | — | — | — | — |
| **a range** on a ring | drag along the ring (an arc; over the origin two runs, EC-9) | Shift+drag | Cmd/Ctrl+drag, Alt+drag | Enter | Alt+Left | — | — | — | — |
| **the selection** | — | Shift+arrows (up and down a line in text; round on a ring) | Cmd/Ctrl+Shift+arrows | Enter | Alt+Left | — | a linked selection scrolls a text to it (`reveal.linked`) | — | copy its letters: `selection.copyLetters`, Cmd/Ctrl+C in `view.keys` (what is shown); the selection's bookmark command (planned, with `refs.remove` and labels) |
| **the view** | — | — | — | 0 fits every section its way | Alt+Left / Alt+Right | — | — | — | undo and redo |
| **a track or layer** | a row by its label in the gutter: planned ([grids-and-rows.md](grids-and-rows.md)) | — | — | planned ([grids-and-rows.md](grids-and-rows.md)) | — | — | — | `track.hide`: its band leaves every line, its ring goes | heights by patch: a band on every line, a ring's width |
| **a block** | — | — | — | — | — | stacked: about the block under the pointer (EC-115) | — | no: a `range` or the blocks list by patch | `copies` merge or split by patch (EC-65) |
| **a bookmark or link** | — | — | — | runs its action: text reveals, a ring moves its window | Alt+Left | — | — | — | — |
| **a navigator**: the ring's window (the `window` target) | — | — | — | a click on the window is the ring's own (select), never a move | Alt+Left | drag its `start` or `end` part: the other end held | drag its `body` | no window when nothing follows | — |

## Configuration

Whose, by the rule the axis design states (spec param, `RenderConfig` preset, `SpecBounds`, `ViewerConfig`, a port, a
registry entry, an option of a call, or fixed with its reason).

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `section.layout` | spec param | the view's `layout` | a view mixes a ring and its detail | one layout a view |
| `group` (text) | spec param | 10 | the reading unit of sequence text (EC-94) | any; 0 for none (Jalview's lines) |
| `numbers` (text) | spec param | `start` | a line is read from its start | `end`, `both`, `none` |
| `group_labels` (text) | spec param | `above` | a number a group reads without counting | `inline` (EC-95), `none` |
| `align`, `max_width` (text) | spec param | `start`, none | the whole width | `center`; a px cap (EC-98) |
| `height` (text) | spec param | none (the content, to `text.maxHeight`) | the writer may fix a box | px |
| `rotation` (rings) | spec param | 0 (origin at the top) | the plasmid maps' convention | any angle |
| `heights` (section) | spec param | `content` | tracks fit their data in steps, one kind one height (decided); `uniform`, `fixed`, `ext:` sizings; `track.height` overrides; the host's `trackSizing` port over all |
| `labels` (rings) | spec param | `outside` | EC-108 | `none` |
| `arrange` (blocks) | spec param | `row` | one axis, one scroll | `stack` |
| `copies` | spec param | `merge` | a copy group's sequence said once, its cost one block's (EC-65, EC-67) | `each` |
| `text.groupGap`, `text.lineGap`, `text.bandGap` | `RenderConfig` preset | half a cell; 8, 2 (`figure` 18, 3) | a space between words; a line block reads as one when the line gap is the larger | per preset |
| `text.collapseEmpty` | `RenderConfig` preset | true | an empty band separates what belongs together (EC-96) | false: uniform lines, `end` exact |
| `text.maxLines` | `RenderConfig` preset | 100,000 | a figure, not a genome; line tops kept per section | larger, at its memory |
| `text.maxHeight` | `RenderConfig` preset | 600 px | a text in a host page scrolls inside its box past it | px; none (the whole figure) |
| `text.pieceBreak` | `RenderConfig` preset | `none` | lines stay aligned; the boundary marks the edge | `line` |
| `rings.hole`, `rings.maxDiameter` | `RenderConfig` preset | 0.35 of the radius; 720 px | room for the title, chords readable; a map, not a wall | 0 to 0.8; px |
| `rings.blockGap` | `RenderConfig` preset | 1° | a gap between blocks reads as a boundary (Circos) | degrees |
| `rings.letters` | `RenderConfig` preset | `upright` | readable at any angle | `tangent` |
| `rings.center` | `RenderConfig` preset | `title` | the name and length in the hole | `none` |
| `rings.labelGap`, `rings.leaderBend` | `RenderConfig` preset | 8 px, 12 px | leaders fan out over the ruler, never through one point | px |
| `rings.window` | `RenderConfig` preset | `followed` | a window means only what a follower shows | `always` |
| `reveal.linked` | `ViewerConfig` (the `reveal` group, with `reveal.cursor` and `reveal.opensClades`) | `text` | a text hides most of itself (Jalview scrolls to a linked selection) | `all`, `none` |
| the window's drags | `BehaviorProvider` (`window`, shared with the navigators) | in the default spec | a host may leave it out | — |
| key of `selection.copyLetters` | `ViewerConfig` (`view.keys`) | Cmd/Ctrl+C | the platform's copy | any |
| the ring's window in an export | spec param of the view (`export.window`, optional), `renderSvg`'s `window` overriding it (decided) | false | an exported map is a figure | true |
| groups bridged; a collapsed piece edge zero width on rings | fixed | — | no second use; one rule for breaks on every layout (decided) | — |

## Tensions

- **Data decides a text line's height** (empty bands collapse, a frame's title row only where a region starts), but
  lines are laid out lazily and data arrives by region. Options: uniform lines (no collapse; scroll exact, figures
  looser); collapse with estimated heights corrected as data lands, the anchor keeping the reader's line still
  (recommended); compute every line's height up front (reads all data). Open: the text viewport's `end` is an estimate
  (lines at full height), so a follower of a text section may show a little more than the text does.
- **A viewport field read three ways** (`start` a window, an anchor, a ring's window): one record keeps the
  `ViewModel`'s groups, follows and history layout-neutral, at the cost of a field whose meaning is the layout's.
  Recommended: one record, the meaning in the layout's `_spec` section; the alternative, a viewport type a layout,
  multiplies the model's code by the layouts.
- **Kept lays drawn again versus lists a section**: a frame that relays one section draws every other section's kept
  lays again (one list a pass). Accepted while the ring's frame target holds; if it does not, lists per section (a
  public-form change, every backend) is the next step.
- **Follow across blocks of several systems** is refused this slice (`seq.render.follow.blocks`): a Complex's chains
  in one section cannot lead a detail. Options: a follow names a block; a follow carries each block along its own map.
  Open until a view needs it: the Complex's PAE is a grid section of its own whose x follows CDK2's section and y
  cyclin A's (one section a chain), so it needs no follow across blocks.
- **Scaled segments** (Circos): one scale a ring, so a region cannot take more of the circle than its length. A zoomed
  ring would put the viewport back into the ring's key. Open: the Genomes milestone decides.
- **The window is shared with navigators**: built in `render/window.ts` here, extended by the navigators; its look and gestures
  are this design's, [navigators.md](navigators.md) follows them.
- **Touch scroll of a text has no gesture**: a touch drag on data selects (the grammar), so a text on a phone scrolls
  only by a two-finger pan read as a wheel. Standing until a touch preset is designed.
- **Scroll groups across layouts**: two sections on equal axes group today; a text and a tracks section on one system
  would then share a range whose meaning differs. Decided: groups only among sections of one layout; across layouts,
  follow.

## Self-check

The claim: "this design covers what a user does with layouts: text, blocks and rings".

1. **Assume it holds.** Every layout is a cut of one linear order of display units into straight strips with a warp;
   layer types, emphasis, hits, the key and export read strips; the viewer navigates through the layout's nav.
2. **What it rests on, negated.** (a) *Every mark lies on one strip.* Negated: marks that join two places (a pair's
   arc between two lines, a chord through a ring's hole, a frame across tracks). (b) *Hit-testing is placement's
   inverse, one point one position.* Negated: points that name no position (gutters, labels outside a ring, the hole)
   and positions in two places or none (copies merged, a ligand chip, a cyclic seam). (c) *The viewport is a window of
   display units the layout draws.* Negated: a layout that shows everything and windows something else (the ring), or
   whose window is set by scrolling and letter size (text); and the window's inverse, the ring moving under a fixed
   window (its rotation). (d) *One leader a section moves the rest.* Its symmetric case: a ring following a linear
   section; its plural: a section of several systems (blocks) under a follow. (e) *What a section shows is chosen in
   it.* Negated: a selection made elsewhere that the text hides (Jalview scrolls to it). (f) *One scale a section.*
   Negated: Circos's scaled segments.
3. **The negated world.** (a) gives section-scope layers with a placement reader, and the frame's `span`; (b) gives the
   kinds of point (`label` hit by its item's key, `gutter` as a row, `gap` for nearest) and merged copies hovered for
   every member; (c) gives the per-layout reading of the viewport, follow reading the leader's viewport and not its
   transform, the window out of the ring's key, and rotation changing the ring's key but not the window's offsets;
   (d) shows the same window either way, so no ring-specific follow, and shows that a multi-system section has no one
   map to follow along; (e) gives `reveal.linked`; (f) has no answer short of a zoomed ring. It also surfaces the
   plural of leaders (per-section viewports in the frame), the inverse of reflow (the reader's place kept by the
   anchor), and a missing view: the rotated linear view, which the cyclic viewport's wrap gives.
4. **Held both.** Integrated in the Decision: `scope: 'section'`, `ctx.place`, `ctx.span`, the label sink; the kinds of
   point; `FrameInput.viewports` and follow from the viewport; the cyclic wrap for leaders; the window in its pass and
   out of the ring's key, with per-section caches; the symmetric follow; `reveal.linked`; from (c)'s plural, a second
   axis (a band's `bandTop` and `rowScale`, a grid's `y`) in the one viewport record and the y path in the emphasis
   passes; from (c)'s inverse, a history step over a changed layout restored by position, never dropped. Planned in [grids-and-rows.md](grids-and-rows.md): section-
   scope pairs and chords, the row by its gutter label. Standing tensions: a text line's height decided by
   data while lines are lazy (the text `end` an estimate); follow across blocks of several systems refused this slice;
   scaled segments and a zoomed ring refused, the Genomes milestone to decide.

## Consequences

- **`seqx-core`**:
  - `render/layout/`: `SectionLayoutType`, `Placement`, strips as a struct of arrays, the translate and polar warps
    with `WarpWriter`; `AxisTransform` with groups and a system a strip; the tracks layout restated on it first
    (records unchanged but for the `window` header), then `text` (lines, gutters, numbers, chips, collapse, frames) and
    `rings` (sectors, labels and leaders); `render/window.ts` (the window mark and its hit targets); blocks in all three;
    `FrameInput.viewports` with `top`, follow from the leader's viewport, the key and caches per section; the `window`
    pass; glyph runs with per-letter y, angle and vertical scale; a band's own scroll and row zoom (`ctx.bandTop`,
    `ctx.rowScale`) and the viewport's `bands` and `y`; the y path for `axes: 2` layers; the layer contract's `used`,
    `scope`, `ctx.span`, `ctx.place`, `ctx.labels`; `seq.render.text.lines`, `seq.render.ring.image`,
    `seq.render.follow.blocks`; `renderSvg`'s `window` option.
  - The language (`0.4.0`) and resolved view: `section.layout` and its params, `block`, `arrange`, `copies`, a path a
    block.
  - Tests: property tests over random sections a layout (every position placed once, `at` the inverse of `place`,
    `offsetAt` monotone in reading order and angle, the emphasis covering exactly the selection's cells, Canvas2D and
    SVG lists equal); a detail following a ring's window and a ring following a detail; an unchanged section not laid
    out again (its lays the same objects) while another moves; drawing against the data at the fit view, across a tile
    edge, across a line break and across the origin; reused buffers compared with fresh ones per layout; records for
    the p53 text figure and the Plasmid ring; the tracks records unchanged as the restatement's guard.
  - Targets ([performance.md](performance.md) gains the rows; each with its path and test): **text reflow** (the path:
    line breaking and line tops, then the lines in view laid out), the p53 figure and a 35,000-residue sequence under
    8 ms a frame, flat in the sequence's length (35,000 against 3,500: at most 1.2×); **a ring's frame** (the Plasmid,
    four rings, about 30 labels: a pan of its linear detail, the ring drawn from kept lays) under 8 ms, labels placed
    only when the ring is laid out; **hover** under 2 ms in each layout; draw and hover 0 bytes warm through the warp. A
    quick representative of each under `itWithin` in the stage that builds it.
  - `_spec`: Versions (`0.4.0`); Render gains "Section layouts" (the interface, strips, warps, kinds of point), "Text
    layout", "Rings layout", "Blocks", "Window"; "Tracks layout" restated as one strip a block, follow from the
    viewport; "Draw lists" (the `window` pass, glyph runs), "Layer types" (`used`, `scope`, `span`, `labels`), "Lanes
    and Hits" (hits through the placement, `block`), "The key" (per section, caches); Resolved view (layouts, blocks,
    paths a block); Configuration.
- **`seqx-viewer`**: `ViewModel` through each section's `SectionNav` (the viewport record's meanings, `scale`, the
  cyclic wrap for leaders, `clamp` per layout, scroll groups within one layout); `input` reads the placement (`at`,
  `offsetAt`, the section under a point), the wheel taken by a text section only while it overflows; the hit target
  `window` with its `part` in the bindings (labels hit as `item`, gutters as `row`); the `window` behavior provider (the
  window's drags); the `reveal` group (`reveal.cursor` moved from `cursor.reveal`, `reveal.linked`);
  `selection.copyLetters` with its key in `view.keys`; history steps restored by position across a layout change.
  Tests: the nav's arithmetic per layout in Node; a twin per story (a text and a ring composed with their linear
  detail, on a fake clock); the matrix cells above as binding tests; back across a layout patch landing on the same
  positions. `_spec`: ViewModel (viewports by layout, nav, history across a layout change), Input, Bindings,
  Configuration.
- **`seqx-examples`**: the p53 text figure; the Plasmid ring with its linear detail and the rotated linear view
  (on rings); the Globins ESPript figure (after text and the rows layer); the Complex's chains as blocks and the
  ligands as chips. Each with drawing checked against the data.
- **The storybook**: a Design page "Layouts" (strips and warps, the kinds of point, the matrix cells), each layout's
  stories with canvas and SVG side by side and their headless twins.
- **Order of work**: the interface and the tracks restatement first (the hot file `render/tracks.ts`; rings and blocks
  start after that commit merges), then text; beside it `FrameInput.viewports`, follow from the viewport, the
  per-section caches, `render/window.ts` and its pass, the warp, rings and blocks. The ESPript figure needs the rows
  layer ([grids-and-rows.md](grids-and-rows.md)).
