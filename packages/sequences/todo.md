# Seq\*: to do

What is not built yet. An item leaves this file in the commit that builds it (or moves to "In progress" when work
starts); a new deferral comes in with what triggers it and where its context is. Phase ids are explained in
`_spec/dev-phases.md`; contracts in `_spec/`.

## In progress

Nothing.

## Next

The move into Mol\* v6's repository is done. Next is Prism over Seq\* (`.v6/seq/plan.md` §3, designed in
`.v6/seq/prism-sequence.md`): the Prism viewer's sequence control, tooltip and legend replacing the default primitives
(the chrome's next look, below), then the plan's §4 follow-ups (whether `viewer/react` stays).

## Deferred features

Each builds on the coordinate model (`seqspec:model#coordinate-maps`, `seqspec:model#pieces-and-layout-maps`) with the
biology at application level, never in the core.

**Input as streams, reconciled with Mol\*'s** (David, 2026-10-14). Seq\* has its own input layer (`canvas/dom-input.ts`,
the dispatcher, its `Binding`), and a host sees interaction only through the context's states and events (`hover.live`,
`selection.current`, `events.menu`); a right-click reaches it only when the `contextMenu` behavior is in the spec.
Reconcile with Mol\*'s: build the DOM input on `InputObserver` (`@molstar/core/util/input/input-observer`, rxjs
streams), one binding vocabulary with Mol\*'s (`@molstar/core`'s bindings), and public interaction streams on the viewer
as Canvas3D's `interaction` has them (hover, click, drag, context menu: the hit with section, track, item, position and
coordinates, the pointer's place, the modifiers), configurable by params; behaviors stay the built-in reactions on the
same streams (input delivers, behaviors react). Under the performance contract (the input path's heap guards), with the
spec (`interaction.md`, `viewer.md`) and tests. By Opus, design first. Context: `.v6/seq/plan.md` §4. _Trigger:_ a PR of
its own, after the one bringing Seq\* and Prism's sequence prototype.

**Tracks in the text view, as lanes** (David, 2026-10-14). The wrapped geometry (`compact`, `wrapped`) draws letters;
other layers either draw as underlines (ranges and features) or draw nothing while still taking their lane's height (a
line track: Prism's pLDDT). Draw most layer types as lanes under each line: ranges and features (secondary structure,
domains, sites) as underline lanes with consistent heights; continuous data (pLDDT, hydrophobicity) as a thin colored
strip, one color per position from the layer's theme; letters as now. Not in text: dot plots, pairs, trees, multi-row
alignment tracks (they need tracks), which then reserve no lane. Under the performance contract: allocation-free lane
drawing, a heap guard per new lane kind with its positive control, tests, `_spec/viewer.md` (the wrapped geometry) and
the generated reference. By Opus, design and implementation. Context: `.v6/seq/plan.md` §4; Prism's Text view shows
letters only until then. _Trigger:_ a PR of its own, after the PR that brings Seq\* and Prism's sequence prototype (that
one only reviews the current renderer changes).

**An item coordinate system's own letters.** A letters track on a transcript's block, derived lazily from the source's
streamed letters through its map (as the codon label reads them); until then a `letters` layer on an item block draws
nothing. Context: `seqspec:sources#item-coordinate-systems`. _Trigger:_ a story that needs a letters track on a
transcript's block.

**Items named before they are seen.** A layer whose data waits for items is noted (`svs.ref.pending`, info, once per
source with the first refs), but a block or a target naming one still waits quietly, the note gives no nearest key (a
ref out of view is no typo until the source has sent all it will, and no source says so yet: a finished stream or a
paged source read whole would; cache the search per source and held-key count), and naming an item does not ask for its
region (requests are by region; an item is found by id only in a region held). _Trigger:_ a view that names a transcript
without a layer reading its source, or a report of a silent typo.

**Per-row targets naming streamed items.** An annotation whose rows' `target` (or `b_target`) name dynamic nodes cannot
wait for them yet (`svs.ref.derived_unknown`, saying so). _Trigger:_ a story whose rows target several transcripts.

**`map` from the genome-file providers.** The gmod providers send no `map`: a BigBed's thick start/end and blocks (or a
GTF's CDS rows) would give a transcript's codon map. _Trigger:_ the genome browser demo showing a protein on a GENCODE
transcript.

**Data streamed against item coordinate systems (B2).** A source whose own `coordinate_system` names item coordinate
systems (`<source ref>/<item id>`): a protein-level stream of per-residue scores for every transcript on a chromosome,
asked in each transcript's residues and drawn along its exons through its map. Context:
`seqspec:sources#item-coordinate-systems`. _Trigger:_ a story that needs per-residue scores for every transcript on a
chromosome.

**SeqQL as command and link targets (QV1).** A `{query}` target on `seqv:select`, `seqv:highlight` and `seqv:focus`, a
`query` argument in the links, and the viewer's `queries` registry kind (which the Mol\* adapter's
`molstarQueryExtension` would take). Context: `seqspec:viewer#commands`, `seqspec:seqql#in-documents`. _Trigger:_ a host
or a story that needs to select, highlight or focus by query (an agent's reply, a description link).

**A dot plot's y with its own window: the next step** (David, 2026-10-13). A plot's y stays its own section's axis
(decided), but with a window of its own: a vertical pan and zoom on the plot move y alone, so an off-diagonal block
(PAE's barnase × barstar interface, which is what readers zoom into) is shown by itself while x keeps the tracks'
window. The plot's y transform reads that window instead of x's; the drag, wheel and keyboard on the plot move it.
_Trigger:_ a PAE or contact-map story that needs an off-diagonal block on its own.

**A dot plot's y from another section, and synteny.** A `dot_plot` draws y on its own section's axis, showing the same
window as x (`seqspec:format#layer-types`): zoomed in it shows the residues in view against each other, never an
off-diagonal stretch alone. A `y` param naming another section of the view would take that section's axis and viewport
(placement of side b on its blocks, the y transform from its axis view, the surface following its viewport, the plot's
vertical drag and zoom moving it), so two sequences or assemblies are each zoomed on their own. Linear synteny is pairs
between blocks of two stacked sections, drawn as ribbons, with the axis links keeping both in step. A pairwise mapping's
segments as data (a provider giving its segments as `pairs`) would draw an alignment's path without precomputing pairs.
_Trigger:_ a story that compares two assemblies or two long sequences, or a PAE read off its diagonal while zoomed in.

**Read pileups.** Each read is a row coordinate system with a gapped map onto the reference (a CIGAR is a gapped map),
packed into lanes as `ranges` with per-read letters; coverage is the count pyramid. Reads stream from BAM or CRAM
through a stream provider; mismatches are a `rules` measure (letter differs from the reference). _Trigger:_ a story that
shows the evidence behind a variant, and a BAM/CRAM provider in `stream/gmod`.

**PCA.** Not a sequence axis: a scatter of a keyed set (samples) by two numeric fields of a key table, with selection
linked to the genotype matrix's states (brushing filters them). It needs a 2-D point layout outside sections, or a host
chart wired to Seq\*'s selection. _Trigger:_ a population-structure story, or linked brushing between samples and the
genotype matrix.

**Hi-C.** A contact matrix over two copies of one axis: a `dot_plot` over a `grid`, its 2-D levels (2^L cells a sample,
as the plot's tiles reduce zoomed out) streamed by region and level from `.hic` or `.mcool` through a provider (a grid
is held whole today, a block pair allocating at most 16M floats, `DotPlotLimits.cells`; past it the pair is left out),
often drawn as a triangle rotated 45° (the painter already rotates). _Trigger:_ 3D-genome stories.

**Render props through `style`: curated keys, no bag.** A section's `style` takes a preset name or an inline preset
whose keys are a curated, semantic, snake_case set mapped onto the viewer's `RenderProps`; a document gets no general
render-props bag. Pixel sizes depend on where a document is shown (a story in a column, a hosted page, an export at
another DPR), so the keys name intent (a density, a size step, a color role) and the viewer decides the pixels for the
place; and `RenderProps` names are internal, so a bag would make every internal name format. A key enters the curated
set when a story needs it and its meaning holds in every geometry. _Trigger:_ a host that needs a render prop no curated
key covers; it adds a key, not a bag.

**`.svsx` carrying stream tiles.** An archive that holds a streamed document's tiles, so it opens offline. Context:
`seqspec:sources#archives`.

**Indels: a range across an insertion, kept in pieces.** Through a partial map a range spans its first to last paired
position (David, 2026-10-13: the hull stays the default), so a wild-type domain on a design spans a grafted loop. A
`ranges` param keeping its pieces (a hairline across the insert) would draw it split. _Trigger:_ a story that needs the
split.

**Boundaries at a block's ends.** An end boundary's glyph (a terminal insertion's) is drawn on the block's edge, half
past it: at the data area's edges it is cut in half, and there its outer half does not hover (the point is off the data
area; in the gap after a block it does). An inset for end glyphs would close it. In wrapped text a boundary is an
underline under the position after it, so one at the end draws nothing, and a caret is drawn on the line's edges only.
_Trigger:_ a story with terminal insertions read closely, or a wrapped text with boundaries.

**Carets: navigator ticks and composed maps.** A channel whose marks are carets only (an insert seen from the other
side, a deletion's junction) puts no tick on the navigator, which reads the primary from runs; and a boundary through a
map of order 2 on the forward leg is not yet checked to land as a caret (the review's case; one map is tested).
_Trigger:_ a reader looking for a selected insert on a long axis, or a design mapped through two maps.

**Match-state numbering on an alignment's ruler.** Insert columns are drawn as such (lower case on a veiled cell, a mark
on the ruler), but the column ruler still counts them; HMMER's and Jalview's match-state numbering (inserts unnumbered)
is a format change (a ruler `numbering` over columns). Meanwhile a `ruler` with `of` the reference row skips them.
_Trigger:_ an alignment story read in match states.

**`emphasis` on a track whose letters are another polymer's.** An `emphasis` layer emphasizes the block's own letters
even when its track's `letters` layer shows another's (`of`): on the wild type's axis a design's substitution emphasized
shows the wild type's letter. Case 15 colors the cell under uncolored letters instead. Viewer only
(`renderers/emphasis.ts`). _Trigger:_ a story that emphasizes on an `of` letters track.

## Follow-ups

**Ligands and other short units in the text view** (David, 2026-10-14, from Prism's sequence strip). In the wrapped
(`compact`) text a structure's non-polymer units come after the chains with spacing that reads as unrelated: an empty
block, then a label such as "B HEM CHA" far below the last line. Lay short units out as part of the text (the unit's
name and its residue, close to the chains, the blocks' rhythm kept) and decide how many of them a strip shows. Trigger:
Prism's sequence view beyond its first slice.

**Dot plots, the rest** (2026-10-13). The keyboard cursor on a plot moves across only (no partner); a square plot taller
than its box's room (a wide screen) is capped there (so it shows whole once scrolled to), its cells flatter than square;
the axes are blocks only (an alignment's columns are not drawn on a plot); a pair's reversed end draws its segment
mirrored, but a pairs row has no strand of its own for an inverted repeat. Two dot plots in one track: the first layer
draws the channels' rows and the drag's rectangle (`DotPlotState.first`), the first with a grid on the hovered block
pair its cell's outline. Blocks not shown whole in one direction (positions, a split layout) have no pyramid: a
zoomed-out sample reads its cells one by one, which the samples budget does not weigh. _Trigger:_ a reader of a plot in
a story, or a host embedding one; a large grid over such a block.

**Large sets, the rest of `headerHover`'s design** (2026-10-13). Built: the behavior (feature layers' shown items,
continuous layers' row band, `maxItems`, `forward`), the size before resolving (`hover.size()`, `lociSize`, the DOM
details), item loci resolved on read, `linkViews`' size-gated reactions. Not yet: a streamed layer's resident items (its
header lights nothing positional; its items through the streams' regions, re-described as they arrive); a header's label
saying "N items" past `maxItems`; the default selection and highlight channels gated by size as behavior params (a
host's large selection lights whole today). _Trigger:_ a streamed feature track whose header a reader hovers, or a host
lighting more than 100,000 items. Decided (David, 2026-10-13): past `maxItems` a header lights the row band alone; the
overlay does not resolve a large set per window in view (the channels' marks follow the view, not the scroll).

**The chrome's next look, in Mol\*** (after the canvas-chrome pass, 2026-10-13). Left for the iterations in Mol\* v6:
the canvas draws a track's label at x 10 (`draw/composer-headers.ts`) and the chrome's text lines up with it
(`skin.test.ts` holds `--seqv-label-x` to the composer's), where Prism's inset would be 8 (a drawing change with the
pixel probes; then one shared constant); the remaining text glyphs as icons (⋯, ‹ ›, ▶ ❚❚, ✕, "Colors ▸", "‹ Colors");
floating the bar and the player over the canvas as Prism's cards (they stay in flow, in the glass's tone; the legend
unboxed under the canvas); motion (menus dropping in, the legend's reveal); a host's composition helpers beyond
`SeqViewerRoot`/`createDomUI` (a DOM part per shown section, as React's `shownSections` map does by hand). _Trigger:_
the move into Mol\* v6 and its Prism components replacing the default primitives.

**A codon text track.** A residue's three bases drawn in its cell on a protein's axis, only where a cell holds three
letters (about 20 px or more), nothing below that. _Trigger:_ a variant story where the reader must see the reference
and alternate codon side by side.

**`domain_group` in templates.** A template track's generated rows declaring one shared domain group (EN1).

**A data hover's row band on wrapped lines and rings.** The `hoverRow` band is drawn on linear sections only: a wrapped
line's row would need the hovered line (the hover carries no entry), and a ring's would be an annulus. _Trigger:_ a
reader of a text view or a plasmid map asking where the pointer's track is.

**The hover card from the keyboard.** The keyboard cursor opening the hover card of what it is on, with its actions
reachable from the keyboard.

**The ideogram navigator** (ST5):

- the bands as document data (a `bands` source, in snapshots, another assembly's);
- a lighter wash outside the window;
- band names on hover.

**Prefetch** (ST6):

- the lead sized by pan speed;
- the lead's cost for items (up to 16 regions a side, so about twice the traffic a pan);
- hints from snapshots beyond the next and the previous.

**A genotypes preset.** `b.preset.genotypes(view, region, { matrix, samples, order, filter })`: a matrix track ordered
and filtered by a sample table keyed by state (KS2's story writes the generic params by hand). _Trigger:_ a second host
or story that builds a genotype matrix.

**Sizes in the header tooltip.** A track header's tooltip lists its items and ramps, not its size legends (EN2b).

**The whole unit from an own row's label.** With `UI.OwnRowHover` `residue` (the default), a chain row's bar hovers a
residue; hovering the row's label could hover the whole chain (HR2).

**Item bins by category.** A provider's `bins` for items may carry `field` and `counts` (a count per vocabulary value),
for density colored by category; the viewer reads only `count`.

**The ring window's double click** (CG3). The window has the navigator's drags, and the ring's `+` and `-` zoom its
section, but not the navigator's double click (the full view, or a zoom there). _Trigger:_ a ring shown without its
linear detail's controls, or a reader asking for it.

**Selecting text across blocks.** A drag in a text view selects within the block it started in (onto another block's
lines it stops at its own block's first or last position), so a text of several chains selects one at a time. The
design: an end block in the pending state (`PendingSelection.endBlock`), drawn per line on each block between, and
`rangeLoci` emitting one `seq-loci` element per unit in display order (the first block from the start, the blocks
between whole, the last up to the end). _Trigger:_ a text story of several chains read across.

**Selecting text by touch.** On touch a drag on the text selects too (the input records carry no pointer type), so a
reader scrolls a box text only by its inset and its right strip (24 px), and a text sized to its content only outside
its canvas: a limitation. The design: a `pointerType` int on `PointerInput`; wrapped surfaces with
`touch-action: pan-y`, so a touch drag scrolls the page or the box; a long press (about 450 ms within the slop) starts a
selection drag, and a tap selects as a click. _Trigger:_ a touch host, or a reader on a tablet.

**The Mol\* adapter** (moving to Mol\*'s repository):

- `molql` returning one group per chain, not one group (`seqspec:libraries#the-mol-adapter`);
- the query result line names units by `label_asym_id`; chain naming is to be solved on the host's side;
- the Mol\* hover sometimes needs a second move after the pointer was over the 3D canvas.

**Nightingale's linegraph in one color.** A `values` line colored by its values (case 8's pLDDT bands) is drawn in one
color there and reported (`svs.render.theme`), as Nightingale gets no new features. If that changes: one linegraph
series per palette color, each holding only its color's positions.

**A strand through a minus-strand codon map** (2026-10-14, VL-08). A `secondary` strand's arrowhead points to the
C-terminal end in the block's own order: placed on its gene through a codon map on the minus strand (residues descending
as block indices ascend), it points the wrong way, because a prepared feature keeps no orientation of its placement.
_Trigger:_ a protein's secondary structure shown on a minus-strand gene (case 21's style).

**Core's stage plan for the viewer's assembly** (2026-10-14, VS-25). The graph stages what `compile` does by its own
lists (`node-kinds.ts`: `computedAlignment`, `readsSet`, `answersNamed`, the three assemblies' groups), checked against
`compile` over core's cases (`viewer/_test/graph-parity.test.ts`). Core exporting its stage plan (ordered kinds and
predicates) for both to drive would take about 50 lines out of the viewer and 15 in. _Trigger:_ after the move to
Mol\*'s repository, or a new stage in `compile`.

**A hover card over an item a reassembly dropped** (2026-10-14, LB-29). A reassembly of the same document leaves the
card open, its actions on the item it opened for; a data edit that removes that item leaves a card whose actions name
nothing (a run finds no item). Close it when its target is gone from the new set. _Trigger:_ item actions in a story
whose data a person edits.

**The conformance kit's hit round trip, the rest** (2026-10-14, LC-17). An item hit must cover the position hit,
centered text is hit-tested and the reversed geometry runs the mirror; not yet: a mutant that drops every other item
(`dropHalf`) passes for ranges, points, variants, pairs, matrix and rows, and path hulls are no hit target (an arc's
vertices are never hit-tested), so a hit naming an item where nothing is drawn is not caught there. _Trigger:_ a new
renderer whose hit shape is not its drawn shape, or a hit bug on arcs.

**The composer proof's corpus, the rest** (2026-10-14, LC-41). `composer-throws.test.ts` covers the draw, dense band,
header and overlay catches; about a tenth of the composer's branch points are still reached by no corpus file (the
accent tint, the ring window across a cyclic origin, the ring labels' report and cached path, the linear and dot plot
drag rectangles, stacked `breaks`, lazy rows' other layers, the tree wedge's row labels). _Trigger:_ the next composer
change proved by the recording corpus (`viewer/_test/proof/`, `vitest.proof.config.ts`): add cases for its branches
first.

**A positive control that plants nothing** (2026-10-14, LC-23's remainder). The spec check reads `controlCaught(` in
code only, but a trivially true control (`controlCaught(() => 1 << 21)`) still passes. Refusing it means `controlCaught`
counting `leak` calls during its measure; the query guards' controls plant through a `(leak Q)` extension and the two
retention guards keep a host's reference instead, so each needs its own count. _Trigger:_ the guard check's port into
Mol\*'s tooling, or a guard found blind again.

**Reruns equal fresh opens, the rest** (2026-10-14, LB-25). `ref-families.test.ts` has a case per computed ref family; a
source's `coordinate_system` and a layer's `of` are not covered. _Trigger:_ a change to how either is resolved.

## Costs and cleanups

Known, measured, unowned: whoever next touches the code takes one, and its commit says so.

- **A manifest that lists only fine levels is not bounded by `Stream.SampleBudget`** (LB-22, VT-03; measured
  2026-10-14): the tiles server zoomed out over all rows with levels up to 2 listed asked 108 tiles and held 14 MB at
  400,000 rows (16 tiles, 2.5 MB with every level), growing with the rows; past `CacheBytes` at millions. The rule
  favors fidelity (sources.md, Levels and the zoom gate). If such a source appears: skip rows (ask every k-th tile row
  under the budget) or draw the coarsest held, and say so.

- **Nightingale's hover chain allocates per position change** (NG-13, by reading; not measured): `mount`'s `move` builds
  the event and its `extent` array, and a host's `extentLoci` sorts, copies and builds loci before `mounted.highlight`,
  which hover.heap.ts guards from there on. `move` dedupes within a position, so it is per position change, not per
  pointer event. The fix: a `highlightExtent(extent)` on `Mounted` through a reused `Extents` (no loci built), and the
  event's extent as reused `BlockRuns` (about 30 lines with a guard). With section highlights, each hover change also
  builds the highlight's text, a `CustomEvent` and a second `change` (`keptUnder`: a second attribute pass and
  `layLabels`). _Trigger:_ a host sweeping a long track in Nightingale, or a measured stutter on hover.
- **Placed ring labels are one per layout record** (LC-08, `draw/ring-labels.ts`'s `placedRingLabels`): the last
  composer to place a section's ring labels wins, an export's (`renderSvg`) included, so a hit on a label reads that
  placement until the screen's composer places again (its next frame). viewer.md's "a module cache keys on all its value
  depends on" does not hold for it. The fix: key by the layout record and the composer (a `WeakMap` per composer), or
  let an export place without publishing. _Trigger:_ an export drawn while the screen shows the rings at another size,
  or a label hit that names a label the screen does not show.
- **`forEachTo`'s scratch never shrinks** (CM-08, `clipped` in `model/coordinates.ts`): per nesting depth it keeps
  `2 × count` Int32s for the largest map ever walked (8 MB after a 1M-segment map). _Trigger:_ a host that walks one
  huge map once and stays long-lived on small ones.

- **The rest of a system's lookups on `Systems`** (`core/src/model/systems.ts`, VS-27): core's own copies in
  `format/resolve/section.ts` (`numberingOf`'s fallback and own lookup, the focus's numbering) read the lists by hand;
  the viewer's `refs.ts` resolves a ref's range args itself (`fromKey`, `narrow`) beside core's position resolution,
  which a `Loci.fromKey` in core would hold once.

- **The resolved view's tables past 16 columns are dictionary-mode objects** (`core/src/format/resolve/view.ts`,
  `pack`): its keyed stores (`out[k] = …`) leave V8's fast properties at the 17th, so since KS2 added `members`
  `v.layers` is a hash table and every `v.layers.x` read a lookup (the layout's scan per track ran three times slower, 7
  to 21 ms an arrival in the templated storm, before `firstLayer` became a table). Built with named fields, or read once
  into locals on hot paths.
- **A line under a continuous palette, a stroke per column** (`renderers/values.ts`): its color changes at nearly every
  position, so zoomed out it draws one subpath and stroke per device column (about 2–4k a frame per track at DPR 2, e.g.
  under Viridis); a banded line changes color only at band edges. Quantizing the color (or a gradient stroke) would
  bound it.
- **The nav and overlay windows on a wheel scroll** (`windows.ts`, `span`): about 340 bytes a step while the wheel
  scrolls a section (found by KS2's review round, which left the scrolled case out of its guard).
- **One batch per data column per frame** (`renderers/columns.ts`; KS1 added `row_stats` and `redundancy` as sources of
  data columns but no per-frame code).
- **`VerticalNav` folded into the composer** as passes.
- **The model's per-row parent array** at a million rows (`Alignment.parentRows` for a lazy alignment; the procedural
  tile server samples on demand and holds nothing).
- **A miniature of a rows track without a vertical overview** (the `seqv:map` miniature naming a rows track with no
  overview column).
- **`defaultColor` reaching every color channel,** not `color` alone.
- **The viewer commands' reference** beside the generated `_spec/reference/viewer.md`, once the `Search` part's commands
  are all registered.
- **Packaging the `.ts` bins under `node_modules`** (`seqv-render`, `svs-*`): before the first npm publish, with a test
  that a bin runs from a packed tarball.
- **Values through a path ask every region its exons reach** at the data level (a DMD-sized gene: 80 to 100 regions of
  4,096 a whole-protein view); items take a coarser level past 16 regions, values have no such bound (ST7).
- **The first zoom to a new letter size** still measures its font (95 `measureText`, two raster reads for its ink
  heights) and makes its glyph atlas cells, a page upload a text batch that made some: frames of 5 to 8 ms on case 15 at
  DPR 2 (they were 50 to 117 ms, the atlas uploading a page per new cell, before 2026-10-13). Measuring a family once
  and scaling its advances per px would drop the measuring (not the ink reads, which are hinted per size). Better for
  the animation (David, 2026-10-13): **a glyph mipmap for animating uncached state**: each glyph rasterized once at a
  few sizes (16, 32, 64 px) and drawn scaled while a zoom animates through sizes not yet in the atlas, then native
  glyphs at the settled size (crisp, hinted). No measuring or atlas page uploads mid-animation, and letters can grow
  continuously instead of stepping at whole px; softer letters only while moving. A WebGL backend would use a signed
  distance field atlas instead. _Trigger:_ a profile with a first zoom's frame over 8 ms in quiet time, or letters that
  should grow smoothly through a zoom.
- **Chrome's Canvas2D bindings box fractional arguments** (measured 2026-10-14, headless Chrome 155): 12 bytes a
  fractional double passed to `lineTo`, `arc`, `rect`, `drawImage` or `lineWidth` (none to `moveTo`, `fillRect` or
  `setTransform`), so a pairs, sector or path draw allocates per vertex in the browser where the recording backend's
  guards see nothing (the backend guard's sector and paths: up to 120 bytes a frame). And Chrome collects on its own
  while drawing, about every 2,500,000 `drawImage` calls and within long runs of `fillStyle` changes, with no JS heap
  growth (canvas/backend.heap.ts says how the guard stays clear of it). Placing by `setTransform` and drawing at whole
  numbers avoids the first, as text does. _Trigger:_ a profile of a pairs or path-heavy view with collections while
  panning.
- **The viewport's animation frame has no heap guard** (`viewport.ts`, `onFrame`): the manual scheduler allocates a
  frame's map in `advance`, so a guard needs a frame port that does not.
