# Pairs, grids and rows

## Intent

The figures milestone's second half ([version-1.md](version-1.md), [features.md](features.md): Pairs, Grid and dot
plots, Logo, Rows layer, Alignments, Trees). Four things a user reads off a family or a molecule:

- **pairs**: two locations that belong together (a tRNA's base pairs, nested and crossing; contacts across two chains),
  drawn as arcs over a track, chords on a ring, dots on a grid (EC-62, EC-63, EC-129);
- **grids**: a value per pair of positions (a base-pair probability matrix, a PAE), drawn as a dot plot with two axes,
  hovered as a crosshair and selected as a rectangle (EC-58, EC-59, EC-109, EC-145);
- **rows**: an alignment's rows over its columns, in tree or data order, clades collapsed to a summary, rows pinned,
  10,000 × 300 held and 10M lazy drawn at the same frame time (EC-38 to EC-44, EC-100, EC-103, EC-106, EC-154);
- **the logo**: a conservation logo over the rows' columns.

The author, 2026-11-18, on navigators: "we forgot to even mention the "vertical maps"? ... it should prolly be built with
the "navigation tracks"": the vertical overview reads what this design gives (below, What a navigator reads). Cost
follows pixels ([performance.md](performance.md)): "what is rendered depends on the pixels on the screen, not the amount
of data".

## Precedents

- **The previous implementation** (the precedent read for this design, never named in the package): a `rows` layer with
  a tree beside it, collapsed clades as wedge rows labeled "first +n", rows drawn from 256 × 256 tile images whose
  samples are a level pair (rows × columns, power-of-two each) of majority class and share, coarse under fine, built in
  place and kept in an LRU; a vertical overview of every row; a lazy alignment holding only bound rows; a `pairs`
  renderer of parabolic arcs in the band, ends clamped a quarter window past the edge; a `dot_plot` *layer* in a track
  whose y showed the section's window; a `logo` of stretched letters, one bar a device column below 2 px. Taken: the
  reductions, the arcs' rules, the crosshair and rectangle, the lazy rules, the wedges. Not taken: the dot plot as a
  layer (here a section with two axes, so its y can follow another window), reductions keyed to a theme's colors only
  (here the data's letters first, the theme's classes where cells are held), a level pair per tile protocol of its own
  (here the store's 1D tiles over the rows system, columns whole).
- **The previous implementation's clades**: collapsed clades are runtime UI state (a sorted set of tree-node indices per
  alignment, never in the document, gone with it); a single click on a tree node (a forgiving hitbox: the nearest
  branch within 6 px, a wedge's triangle included) toggles it, so does Space on a row header (a wedge's clade, else the
  leaf's parent) and the `clade.toggle` action; a collapsed clade is one display row, its wedge drawn in the tree and
  its label "first +n"; the keyboard cursor keeps its data row and stands on the wedge while folded; clades fold only in
  tree order. No collapse-all, no depth or threshold, no reveal of a row inside a folded clade, no summary cells (the
  wedge row is the tree's; a `summary` track showed under a collapsed *group*, EC-103). Taken: the gesture, the hitbox,
  the label, the cursor's rule, session state by default. Added: summary cells, bulk commands, reveal, the spec's
  opening state.
- **Jalview**: ids panel, tree panel, hide and group rows, block selection, annotation rows under the alignment, an
  overview window; in the tree panel, a click on a node selects its group, a threshold line drawn across the tree
  groups the rows by distance, and collapsing a node folds its rows (its state the session's). Taken: the header with
  tree, block selection, hidden rows that keep coordinates, the threshold as a command. Not: editing, per-panel state.
- **Tree viewers (iTOL, Dendroscope, FigTree)**: collapse a clade to a triangle, collapse by depth or by a support or
  distance cut, expand all; the collapsed set saved with the figure. Taken: the bulk commands, saving the set as the
  spec's opening state.
- **AliView**: very large alignments drawn to the pixel, keyboard first. Taken: zoom out to sub-pixel rows at frame
  time, the cursor. Not: editing.
- **MSAViewer (BioJS) and react-msa-viewer, NCBI MSA viewer**: a conservation bar and logo above, an overview bar, an
  anchor row kept on top, sort by label or field. Taken: the logo pinned above, pinned rows. Not: a row moved out of the
  tree to pin it (here it stays in place too).
- **Skylign, WebLogo**: letter heights by information content, the spec naming bits or probability. Taken. Not: a logo
  as a service apart from the view.
- **R2DT, forna, VARNA, R-chie**: base pairs as 2D drawings (R2DT, forna, VARNA) and as arc diagrams (R-chie, pseudoknots
  below the line). Taken: arcs, crossing pairs on the other side as an option. Not: 2D secondary-structure drawing (not
  an axis; a host's part).
- **Contact-map and PAE viewers** (AlphaFold DB's PAE plot, Juicebox, HiGlass): crosshair, a rectangle selecting two
  ranges, upper and lower triangles of different data, a 2D pyramid. Taken. Not: a rectangle that zooms (here it
  selects; Enter or a double-click focuses it).

### What a user can do there

| In the precedents | Here |
|---|---|
| click a row's name to select it; Shift, Ctrl to extend, toggle (Jalview, every viewer) | built |
| drag a block of rows × columns (Jalview, AliView) | built: a rect |
| click a column on the ruler or logo to select it across the rows | built: a columns range; a drag in the ruler's upper part selects columns, its lower 40% is the ruler rail's and pans (decided, [navigators.md](navigators.md)) |
| arrow keys over cells and rows, a page at a time (AliView) | built, over lazy rows too |
| sort by tree, by label, by a field (Jalview, MSAViewer) | built: tree, data, a field of the rows' table; joined tables wait for keyed tables |
| sort by pairwise identity, compute a tree (Jalview) | refused here: computed by a runtime the application binds ([maps.md](maps.md), computed sourcing); the order then reads its field |
| collapse or expand a clade by its tree node (the previous implementation, Jalview, tree viewers) | built |
| see a collapsed clade's summary (Jalview's group consensus) | built: summary cells |
| collapse all, expand all, collapse to a depth or a distance cut (tree viewers, Jalview's threshold) | built: commands |
| find a row inside a collapsed clade (search, a linked structure) | built: reveal opens its clades |
| save the collapsed set with a figure (tree viewers) | built: `clades.keep` writes the folds shown into the spec's opening state; a story step's view |
| hide a row, show hidden rows again (Jalview) | built: `rows.hide`, a count chip that shows them |
| hide columns (Jalview) | waits (decided): no command this slice; the axis's pieces could carry it later |
| pin an anchor row on top (NCBI) | built |
| a control beside each row (the previous implementation's row slot) | planned: the rows in view with their boxes, for a host to place controls; Seq\* draws none |
| an overview of the whole alignment (Jalview, MSAViewer, NCBI) | planned ([navigators.md](navigators.md)): the vertical overview, from what rows give |
| edit residues, slide gaps, realign (Jalview, AliView) | refused: data is not edited in the viewer ([ux-foundations.md](ux-foundations.md), the matrix) |
| hover a base, its partner lights (forna, VARNA) | built |
| crossing pairs drawn below (R-chie) | built, a config item |
| 2D secondary-structure drawing (R2DT, forna, VARNA) | refused: not an axis; a host links its own drawing through addresses |
| crosshair hover, rectangle selection lighting both ranges (PAE viewers) | built; in 3D through the structure link (the fourth slice, a library: [structure-syncing.md](structure-syncing.md)) |
| upper and lower triangles of different data, mirroring (RNAfold's dot plot) | built |
| a zoom box (HiGlass) | refused as a drag (a drag selects, as everywhere); a selected rectangle is focused by Enter or double-click |
| a color-scale slider (PAE viewers) | the theme's domain and legend ([themes.md](themes.md)) |
| a drag-select past the edge scrolls the view (Jalview, AliView, text editors) | planned (rows here, text in [layouts.md](layouts.md)): clock-timed, faster further past |

## Decision

### The author's decisions

The author, 2026-11-18: "grids-and-rows.md: ok as recommended"; in general: "a lot of this is "configurable behavior"
and where appropriate, the code should allow the configuration via spec ports etc."; on sizes: "don't overdo it on the
limits (they were sort of self-imposed and estimate), but obviously try to minimize the amount of code". On the review of the decisions: "ok" to its recommendations.

1. Clade folds are session state (the `ViewModel`'s, no undo step), the spec's `collapsed` their opening set, and
   `clades.keep` writes the set shown to the spec; pins, hidden rows and the order stay spec edits.
2. A click on a tree node toggles its clade's fold by default (a preset can make it select); a clade's rows are
   selected by its summary row or Cmd/Ctrl+click.
3. A reveal (a focus, a link, a search) opens folded clades; back folds them again.
4. Hiding alignment columns waits: no command this slice.
5. The selection's bookmark command (its name and key) is planned beside `refs.remove` and labels.
6. `mirror` is a plain `false`, set by the tRNA builder; a history step is restored by position, never dropped.
7. A drag on a ruler, the author: "i think drag on a rule should pan. or rather .. the original solution solved this by
   showing the rail as well ... maybe that's the way to go?" Adopted: on an alignment's ruler the lower 40% is the
   ruler rail's (pans) and the upper part selects columns; a grid's axes pan.
8. The author, on the remaining questions: "others as recommended": figure-meaning settings (the pairs layer's `side`,
   `reach`, `crossing`, the logo's `scale`, `correction`, `gaps`, `rows`, the rows layer's `summary`, `numbers` and the
   like) are spec params, not viewer config; `selection.copyLetters` copies a rows selection as FASTA lines.

### Through the layout interface: strips, a second axis, a band's scroll

The layouts design cuts a section into straight **strips** (one `AxisTransform`, its bands, a warp to the page); a
layer is laid out per strip (`ctx.axis`, `ctx.box` its band in the strip's frame), or once a section with
`scope: 'section'` and `ctx.place` for marks that cross strips. It owns `FrameInput.viewports` and the whole viewport
record: the section's `top`, **a band's own scroll and row zoom** (`bandTop`, `rowScale`, by band; `ctx.bandTop`,
`ctx.rowScale`), the grid's `y`; the y path in the emphasis passes; and the glyph runs' per-letter y, angle and
vertical scale. The rows' zoom limits are the host's (the author, 2026-11-21: "seems like this
needs to be configurable"; FX16): `rows.zoom` in `ViewerConfig` (`minPx` 1/64, `maxPx` 40, `enabled`), the floor the
larger of `minPx` and a lazy source's device pixel a row while it serves no row levels, a stored `rowScale` clamped
to them; `enabled: false` leaves a pan-only rail. Open: the source's coarse row levels; a band's own override in the
spec. This design builds on them and adds only:

- **The grid is a section layout** (`grid`, a `SectionLayoutType` beside tracks, text and rings): one strip, the plot,
  whose transform is x and which **carries a second transform** `strip.y` (an `AxisTransform` down the strip's frame,
  absent in every other layout), read by a layer as `ctx.axisY`. The `grid` layer and the `pairs` layer laid out in a
  grid section need it; on a strip without one they report `seq.render.layer.axes` and draw their fallback. `at` gives
  the strip-frame point, so a layer's hit inverts both transforms. Its viewport is the layouts design's record, its `y`
  set only when y has its own window (read by the grid's `nav` and `sameKey`); with `follow_y: x` (the default) the
  field is absent. Why a strip property: the layer stays a straight-axis layer; only the grid has two axes.
- **The rows layer owns its y inside its band** (the tracks layout knows bands, never rows): its rows transform is
  built from `ctx.box`, `ctx.bandTop` (the rows' scroll) and `ctx.rowScale` (px a row, row zoom), both the layouts
  design's band entry in the viewport, over its own order (prepared, tier 1); its hit inverts it.
- **An order under the rows' pieces**: a permutation from display offsets to rows, implicit for data order, an
  `Int32Array` and its inverse otherwise, so folded clades and the pinned band are pieces over permuted offsets. Why:
  tree order is no list of runs.
- **Hits stay flat** (0 bytes warm): `Hit` and the hover record gain `ySystem` and `yPosition` (a position on the
  second axis: a row of the rows system, a grid's y). A rows cell and a grid cell are the target `cell`, with `ySystem`
  and `yPosition` set; the viewer gains the targets `row` (a rows header), `tree-node` (a tree node) and `summary-row`
  (a folded clade's row). Why: a rows cell and a grid cell are a cell with a position on each axis, and a tree node and
  a summary row take different gestures (a node toggles its fold, a summary row selects its rows).
- **Emphasis along y**: the layouts design paints the section's emphasis per strip along x, and gives a layer with
  `axes: 2` the emphasis records' parts on its y system (its y path); the rows and grid layers paint them (a row's
  band, a block's rect in the rows, the crosshair's rows on a grid) through their own y transform, kept out of the x
  painting.
- **The rows' window** (history, sessions, navigators, the cursor) is derived from `bandTop` and `rowScale`; sessions
  and history steps store it as a row key (the row at the top) and `rowScale`. A step over a changed axis, layout or
  order (one rule, in every design that touches history) is restored by position: the top row by its key, clamped to
  what is drawn (a row now folded to its summary row, a row filtered out to the nearest shown); never silently
  dropped.
- **Per strip, by layout**: in `tracks` the rows band scrolls and folds; in `text` (the ESPript figure, the layouts
  design's "Alignment figures") the rows layer is a band of every shown row on every line: no band scroll, no tree,
  every clade unfolded (`seq.render.rows.text.unfolded`, info, when some are folded), a lazy alignment takes its
  fallback (a text figure is of rows held). Why: a figure's lines repeat the rows; a fold there would differ per line.

### Pairs: arcs, chords, dots

- **Data**: a `pairs` annotation, inline (a tRNA's pairs, Stockholm's `SS_cons`, a complex's contacts: read whole), each
  item two ends, each a range on its own system (one or two systems), an open `type` (the column's name: `kind` is reserved in the language; canonical, wobble, non-canonical,
  contact), a `strand` per end (a reversed side is valid, EC-206), `id` or `row`, fields. Self-pairs allowed; an end
  never wraps (EC-10); found from either end through an index on both (EC-62). Pairs read by region (EC-208):
  [region-reads](region-reads.md) (the fifth slice).
- **Scope**: `scope: 'section'` (an arc may join two blocks, two lines of text, or cross a ring's hole). Each end is
  placed with `ctx.place`, and the apex the same way as a point in its band's frame (the mid offset at the arc's height)
  mapped through the strip's warp: one curve is an arc over a block, an arc between two lines, a chord through a ring's
  hole, with no layout named in the layer. An end not placed in the view (a chain not shown) leaves a caret at the placed
  end, never a dangling arc. **An end on split copies** is placed by unit: on the block of its system or its row's unit (`a_unit`, `b_unit`,
  [molecules.md](molecules.md)), the block its copy is merged into, or a unit of its molecule (EC-48). A pair within
  one system is drawn on each copy (the resolver's groups: one system, one range). A pair between two systems shown as
  split copies whose data names no unit is never placed on a guessed copy: it is drawn on no copy, and
  `seq.render.pairs.copy` (warning, once, with the count) says its data names no units. Why (the author, 2026-11-20:
  "the whole point was to show what residue copy has contact with what chain. Imagine having 1tqn assembly. there are
  also few HEMs, but each has contact with just 1 copy of the chain"): a contact on the wrong copy states a fact the
  data does not hold. Contacts derived from coordinates are derived per pair of units. In a grid section the layer is laid out per strip instead, as dots over `ctx.axis` and
  `ctx.axisY`.
- **Arcs**, EC-129: in the band (the layer's `side`: `above`, `below`), never flat; an end past the window taken a quarter
  window past its edge (`reach`, 0.25), so an arc moves continuously on pan; height by span against the widest in
  view, so nested arcs nest and the widest fills the band (`height`: `span`, or a size theme); a label only under
  an apex in view whose ends are within reach. **Crossing pairs** (pseudoknots) are drawn crossing (`crossing`:
  `same`), or the crossing set on the other side (`split`: the largest nested subset above, the rest below; exact up to
  2,000 pairs, greedy past it, fixed: the exact subset is quadratic and no figure asks for another bound).
- **Kinds at one cell** (EC-63): the cell is cut into a slice a kind, in order, each slice hit as its own pair (built
  and tested; decided by the author, 2026-11-20: "Accept both"). A glyph per kind with "+" when they do not fit, the current Seq\*'s, is in
  `packages/seq/TODO.md`.
- **On a ring** the band is inside the innermost ring (`side: below` reaches the hole), so the apex placed toward the
  centre makes a chord, its depth by span. **In text** arcs between lines only with `acrossLines` (false: a line's
  own pairs, the others' ends marked). Fallback: `ranges` (each end a mark).
- **Bounded by the screen**: at most `pairs.maxArcs` (4,096); past it, arcs binned by their ends' device columns drawn
  once a bin (widest first), the rest a density strip of ends, said once (`seq.render.pairs.limit`, info). Contacts on
  copies fold by copy group (EC-62, EC-65): one arc a group, hit as every member.
- **Hit**: an index of arcs by hull; a point tests the candidates over its x by distance to the curve (`pairs.hitPx`,
  4). A hit is the pair's element, its extent both ends; hover and selection light both (a pair closes to its ends),
  projected on linked sections, and through the structure port later.

### The grid section and the dot plot

- **A section of its own**: `section` with `layout: grid` (stacked with track sections in a `tracks` view): `x` (`of`)
  and `y` (`of_y`, the same system or another), each may follow a section (`follow`, `follow_y`; `follow_y: x` by
  default, so y shows x's window and cells are square at aspect 1). Why a section: y needs a window, and a window is a
  section's. The Complex's PAE across two chains is one grid section whose x follows CDK2's section and y cyclin A's
  (`follow_y: <cyclin A's section>`); a single chain's PAE follows one section with `follow_y: x`. **A grid layer
  draws only on its own systems** (its data's x and y): data on another system is not projected onto a grid, so a
  generator computes a grid on the systems it is drawn on.
- **A grid over a complex** (the author, 2026-11-19: "can we have a PAE plot version for complex whose tracks are
  rendered next to each other?", then "ok as recommended"; AlphaFold's and PAE viewers' whole-complex plot): a grid's
  `x` and `y` may each follow a **section of blocks** (`follow`, `follow_y` naming it): the axis is the blocks'
  display axis, each block on its own system, the chains side by side in the section's order (copies split or merged
  as the section draws them). The data is one matrix over the blocks' positions in that order (an AF-style N×N
  `predicted_aligned_error`, read from a JSON or MessagePack array of arrays by path, [files.md](files.md)), its size
  checked against the blocks' total; each cell belongs to a (block, position) pair on each axis. Chain boundaries are
  the blocks' own gap, as in the tracks above, with no line in it and the cells clipped exactly at each block's edge
  (decided 2026-11-19, the author: the gap and an off-centre hairline together read as "very weird"); over two sections of blocks (a dot plot, whose ground is mostly empty) the gaps are tinted a step darker than the ground so the block-pair cells read as cells (the author, 2026-11-22: "Tint the gaps"); a hairline
  (`seq-break`) only where blocks touch with no gap, the blocks' labels on both axes; a hover names both residues with their
  chains and lights them in the tracks; a dragged rect is two ranges, each split at chain boundaries into a part a
  system. **A summary of a grid over a rect** (`gridSummary(view, rect)`: mean, min, max, count of the cells in view
  data, a library function) gives a host's readout ("Aligned A 360–821, Scored A 265–831, mean 10.1 Å"); the panel is
  the application's, a story draws one. Why blocks and not a new kind: the section of blocks already orders chains and
  copies on one axis, and the fourth slice's `complex`/`unit` kinds lay out as blocks.
- **Layers in it**: `grid` (a value per pair of positions, its color theme at its value) and `pairs` (a dot per pair of
  single positions at every zoom, a diagonal segment for a stretch; `mirror`, a spec param, `false` by default, draws
  (i, j) and (j, i) on a self grid: the tRNA builder sets it). Why `false`: a pair is drawn where its data puts it; a
  default read off whether the two axes are one system would change the picture when a writer moves y. A layer may keep
  to a triangle (`triangle`: `both`, `upper`, `lower`). Later layers over earlier; a missing value is transparent, so two
  grids on one pair of systems merge later over earlier (EC-59).
- **A plot's box is fixed** (decided, the author, 2026-11-19: "Fixed box"; found in use: zooming one axis of the
  Complex's PAE resized the plot, the page moved under the pointer and the wheel fell through to it): the box is set
  from the width and `aspect` below and never changes with a zoom or a pan; when x and y show different spans the cells
  are rectangles. Why: a box that follows its windows moves the page under the pointer, as Juicebox, HiGlass and maps
  avoid. **Double-click on a cell zooms a step** about it on both axes, keeping their ratio, again further, as a map
  does; Alt+Left goes back, `0` fits (the tracks' double-click still focuses an element). Why: focusing both positions
  jumped each axis to its residue, so a second double-click re-centred instead of zooming. On a grid whose y follows x (one window: tRNA, a complex's
  contacts) the step keeps the point under the pointer on x, y following; a double-click in a gap between blocks zooms
  about that point as on a cell, never fits (the orchestrator, 2026-11-20, PZ's question). **The wheel over a plot is
  the plot's**: plain wheel pans (dx x, dy y), Cmd/Ctrl+wheel and pinch zoom both axes about the pointer, and at an
  axis's end it stops: the page never scrolls from inside the plot (from outside it, as any page).
- **Size and aspect** (EC-109): height = the data box's width ÷ `aspect` (1), capped by the larger of `heights.grid` and nine
  tenths of the view's room (decided at PG's review: a view sized to its content otherwise feeds its own height back
  and the plot collapses), at least two rows; the width read without the scrollbar's room (the flicker rule).
  `aspect: none` leaves the height to the section (the language has no null here).
- **Direction**: an axis reversed (the axis design) flips its side of the plot; hits and the rectangle mirror (EC-145).
- **Grid data** (EC-58): an inline annotation, a-major `float32` (typed arrays by reference), `null` or NaN missing, a
  wrong length dropped with the count it expected; at most `maxBlockValues` (`DefaultMaxBlockValues`, 2^24) values a
  pair of systems, the bound the store already checks. The tRNA's 76² and a complex's PAE (a few hundred thousand to a
  few million values) fit. Grids read by region (Hi-C): [region-reads](region-reads.md) (2D tiles).
- **Drawn**: at `grid.rectMinPx` (12) a cell and above, rects with outlines; below, `image` primitives a tile (256 ×
  256 samples), never smoothed, from **the shared 2D reduction** (`render/pyramid.ts`, below: square levels, factor 2,
  a reducer by the layer's `aggregate`: `mean`, `max` or `min`; a third more than the grid), built once from the loaded entry
  and kept by its version, at the level whose sample is at least a device pixel (`k`, as for values). Images by
  (layer, level, tile, theme version) in **the shared image LRU** (below); a tile not built in the frame's budget
  (`grid.buildBudget`) draws the coarser level meanwhile, never empty.
- **Crosshair** (EC-145): a cell is two positions. The hover holds both positions and the cell's rect; the grid layer
  paints, for the hover, the selection and every highlight channel, a column at each lit x (the layouts design's x
  painting) and a row at each lit y (its own, the y path), so a hover on a tracks section shows as a cross on the grid
  and a range selected there as two bands; and the cell's outline. Linked sections light both positions; the hover
  card names both and the value.
- **Rectangle**: a drag selects `(rect :a <x range> :b <y range>)`, a click one cell's rect; Cmd/Ctrl+drag adds,
  **Alt+drag subtracts**: a rect minus a rect over the same two systems is at most two rects whose sides hold several
  ranges (decided at PG's review: exact, fewer parts than four), so the Location stays
  exact (the Location rule changes: rects subtract, never merge). Projected onto linked sections, a rect lights its two
  sides.

### One reduction and one image cache, shared

- **`render/pyramid.ts`**: one 2D reduction over a typed grid of samples, its levels by a factor and a pairing the
  caller gives (square levels for a grid; independent row and column levels for the rows, so a row stays a row zoomed
  out on columns) and a reducer (`mean`, `max`, `min`, or majority with its share); reused buffers, built in place. The
  grid, the rows layer and the navigators' rows image use it; the 1D max pyramid of the logo's bars is its one-axis
  case. Why one: three copies of the same level arithmetic would otherwise be written.
- **One image LRU** for every `image` tile (grid, rows, navigators), on core's `LRUCache` (`@molstar/core/util/lru-cache`),
  keyed by (layer, level, tile, theme version), its capacity `RenderConfig.imageBytes` (64 MB) over a tile's bytes
  (tiles are one size). Why one budget: two caches of 32 MB each could not lend each other room.

### The rows layer

**Data.** An `alignment`: a columns system and a rows system (rows are a system, EC-35: positions are rows in data
order, labels their names, `index` and `ordinal` numberings). Its cells are the **`cells` shape on 1D tiles of the
rows system, columns whole** (`codes` uint8, row-major, a tile's rows × every column; at a coarse row level, a row bin
per column as `class` and `share`). Why columns whole: an alignment's width is hundreds to a few thousand, so a tile of
1,024 rows stays within `maxBlockValues` up to 16,384 columns (wider is `too-large`, a MAF's case, Genomes), and the
store's keys and tiles stay one-dimensional. A rows table (`label`, the first residue's number from `NAME/27-137`,
EC-39, fields) comes on the same tiles. Rows bound to sequences (a row's `alignment` map) are held only when something
needs them whole: a path's endpoint, an annotated or reference row (EC-41); paths search bound rows only, a row bound
twice is ambiguous (EC-42). Gaps `-` and `.`, insert columns from Stockholm's `#=GC RF` or all-lowercase columns, A3M's
lowercase on no column (EC-38); readers' errors reported, never thrown (EC-40). A `tree` (Newick) matches leaves by
label, else by the label before its last `/`; a tree that fails keeps the alignment (EC-43).

**Column measures** registered for themes ([themes.md](themes.md) reads them): `share` (the most frequent letter's share in a column,
gaps counted by the logo's `gaps`), `identity` (the share equal to the reference row's letter), and per cell
`matches-reference`. Computed once from held cells, kept by the entry's version; over a lazy alignment, served or
unavailable (`seq.render.rows.measure`).

**Order** (EC-106): `order`: `data`, `tree`, or a field of the rows table (`{ field, direction }`, missing last, ties
in data order); `filter` by a field's values or range; `hidden` (row keys). The tree is drawn, and clades fold, only in
tree order (a clade is a band only there). A lazy alignment is ordered `data`, or `tree` when its source serves rows in
leaf order; any other order is `seq.render.rows.order` and draws in data order (why: ordering needs every row's field).

**The rows axis** (tracks): display rows over the expanded order, one display row a row, or one for a folded clade.
Rows in view are found by arithmetic, O(1) a row, never a scan.

- **Clades** (EC-100, EC-103) are named by the node's id where it is unique in the tree and not a number, else their
  first and last leaves in tree order (`(clade :of [a b])`), so a name survives a reorder of other clades and a tree
  read again. A support label is a number or numbers joined by `/` or `|` (IQ-TREE's `87.3/95`); a `collapsed` name that
  names no clade of the tree is reported once (`seq.render.rows.clade`, at render: the check never sees the tree, which
  is data), never dropped silently (decided at GD's review). Why the exception (decided, the author, 2026-11-19): inner labels are often support values that repeat
  (`((A,B)95,(C,D)95)`), and one name for two clades sends a hover, a fold or `clades.keep` to the wrong one.
- **Pinned rows**: `pinned` (row keys), drawn in a fixed band at the top of the rows band, and kept in their place in
  the order too, marked (`rows.pinInPlace`, true). Why both: the tree stays a tree and clades stay bands.
- **The band**: tracks before the rows track (a ruler, the logo) and after it (a summary) stay fixed; the rows scroll
  in their band (`ctx.bandTop`), which takes the section's room left (the layer's `height`: `fill`, or a number of rows). One
  rows layer a section (`seq.render.rows.second`): the cursor has one row order.
- **The row header**, in the view's one label column (ux-foundations section 3): the tree (`treeWidth`), the label
  (faded past its width), icons in fixed slots (pin, a bound row's link), data columns on the left (Row columns), the first and last residue numbers in gutters (`numbers`): of the residues **in view**, so they follow a
  pan and a zoom (decided, the author, 2026-11-19: "the offsets at the end of the sequences don't update"; Jalview's
  and ESPript's per-line numbers), equal to the row's own at the fit view; a residue is counted, a gap never; a row
  with no residue in view shows none. Labels,
  numbers and icons only at `rows.labelMinPx` (8) and more, the label column empty below it (confirmed, the author,
  2026-11-19: "Yes, leave it empty"; a dark slab was drawn there); a header hover lights the row band (drawing only, EC-146,
  EC-147). Hidden rows are a count chip ("3 hidden") at the header's foot; its click shows them (one patch).

**Row columns** (the author, 2026-11-19: "an extra column with some annotation next to each header… and being able to
position it say to the right of the sequence data"; the base's data columns, Jalview's annotation columns):

- The rows layer's `columns`: a list of `{ field, side, as, width, title }`. `side`: `left` (in the header, between the
  label and the numbers) or `right` (after the end numbers, before a vertical navigator); `as`: `text` (the default,
  and an unknown one's), `number` (right-aligned, its format), `bar` (a value's share of the column's domain) or `chip`
  (a category's color from a theme). A param of the layer, within `0.4.0`.
- **Fields** are the rows table's: the alignment's own per-row fields (Stockholm `#=GS` lines: `OS`, `DE`, `AC`), the
  measures computed per row (`length`, `gaps`, `identity` to the reference row), and a source's fields joined by row
  key (a records source whose `target` is the rows system, [files.md](files.md)). A field missing is an empty cell; a
  field named by no row is reported once (`seq.render.rows.column`).
- **Titles** in the band above the rows, beside the fixed tracks' labels; a click on one orders the rows by it
  (`order: { field }`, a second click the other direction). Hover names the value; a column is hit as its row's
  header. **A cut value shows whole only under the pointer** (decided, the author, 2026-11-19: "the full labels
  should only display when you actually hover on the species column ... otherwise it obscures the data"): a column's
  cell cut short is drawn whole (over the data, to the view's edge) only while the pointer is on that cell; a hover
  elsewhere on the row (its data, its label) tints the row and leaves its columns cut; the keyboard cursor on the
  row's header draws them whole, as the pointer on them would. Why: the expanded text covered the row's own residues
  while the user read them.
- **A narrow frame gives way in one order** (decided by the author, 2026-11-20, "As proposed"; a miniature before names,
  names before data; each step has a floor, never a column squeezed to nothing): (1) a vertical overview narrows from 96
  to 48 px; (2) below a plot of `navigator.overview.minPlotPx` (160) it is hidden, the rows rail staying; (3) the label
  column shrinks to its floor (80 px), then the tree collapses to its toggle glyph at the band's top; a click on it or
  Shift+T (`rows.tree`) opens the tree whole for the session (a band's `tree`, the frame scrolling) and folds it back;
  (4) the plot keeps its floor and the frame scrolls horizontally, never a hairline plot. The step is decided on the
  columns at their own widths, so a navigator never flips with the columns it caused to shrink.
- **Narrow**: columns yield before the label, right ones first, each down to a floor and then hidden (EC-107), a chip
  in the title band naming what it hides, on the side those columns live, never over the data ("Species, Identity
  hidden", "Species +1 hidden" when short of room; decided 2026-11-19, the author: "the 2 columns hidden is very
  misleading, because it doesn't show what columns … when you hover it": beside an alignment a count of "columns"
  read as alignment columns), its hover card listing each hidden column and why ("no room at this width"); the SVG and text figures draw the same columns.

**Opening and closing clades** (the author, 2026-11-18: "we should not forget about closing and opening clades too"):

- **Where the state lives** (decided): the folded set is session state, as the viewport is: no undo step. The
  `ViewModel` holds it beside the viewports (by the rows layer's address), and the `AddressModel` reads and sets it as
  addresses with the rest of the session (`(clade :of [a b])`); its opening value is the layer's `collapsed` (a spec
  param). `clades.keep` writes the set shown into `collapsed` (one patch). Why session: a reader folds as often as they
  scroll, and an undo stack of folds buries the edits. **A patch or an undo that changes `collapsed`** resets the
  session's set to it (an edit of the opening state is meant to show); other edits keep the set. **A story step**
  stores the folded set in its view (a step's view is data); its transition runs `enter`/`leave` channels on the rows.
  **An export** draws the set it is given (`renderSvg`'s `folds` option: the viewer passes its own, as it passes the
  selection), else `collapsed`.
- **Gestures and commands**, each with its inverse:

  | Action | Gesture, key | Command | Inverse |
  |---|---|---|---|
  | toggle a clade | click on its tree node (`tree-node`: 6 px to the nearest branch, a wedge included; decided, a preset may bind it to select instead); Space at the cursor on a node, a summary row, or a row's header (its parent clade) | `clades.toggle` | the same |
  | open a summary row, show its rows | double-click on it; Enter at the cursor on it | `address.focus` on the clade | double-click again, or Alt+Left |
  | collapse all, expand all | the header's menu | `clades.collapseAll`, `clades.expandAll` | each other; Alt+Left restores the set before |
  | collapse to a depth | the header's menu | `clades.collapseTo { depth }` | `clades.expandAll`; Alt+Left |
  | collapse by a distance cut | the header's menu (a drag of the cut line across the tree is planned, `packages/seq/TODO.md`) | `clades.collapseTo { height }` (every largest clade whose leaves lie within `height` of its node, by branch length) | `clades.expandAll`; Alt+Left |
  | expand the clades holding rows | the selection's context menu; done by a focus, a reveal, a link or a search step landing inside | `clades.reveal(rows)` | Alt+Left folds them again |
  | keep the folds shown as the opening set | the header's menu | `clades.keep` (one patch on `collapsed`) | undo |

  Bulk commands and reveal are jumps of the view's history, so back restores the set and the rows' window together; a
  single toggle is not (as a hand's move). A toggle by hand animates the rows' heights as a jump does (`view.motion`,
  `view.jumpMs`; a cut under reduced motion).
- **The tree's look** (decided, the author, 2026-11-19: "the (collapsed) clade visuals are quite ugly. and the different
  line thickness too"; "Light wedge"): every branch one width, snapped to device pixels (a 1 px line centred on a
  fractional coordinate smeared over two device pixels, so widths varied); a folded clade's wedge a triangle filled
  with a light tint of the branch ink and outlined in the branch ink at the branches' width, from the clade's node to
  its deepest leaf's x (by branch length, as its leaves would sit), its far side about half the row high with softened
  corners; hovered, its fill a step darker. Why: the solid triangle at the tree's full width read as a heavy blade.
- **The summary row**: the label "first +n" beside a wedge; its cells the clade's per-column consensus (the layer's `summary`:
  `consensus`, majority letter opaque by share; `logo`, a small logo; `bar`, the non-gap share), **computed exactly
  from the cells when the clade folds** and cached by (clade, entry version, theme): folded clades are disjoint, so
  folding all costs at most one pass over the cells (3M at 10,000 × 300, about 10 ms); a lazy clade summarizes from
  the served row bins (approximate, said once). Framed as one discrete mark, square where the view cuts it (EC-100).
  **Hit**: target `summary-row`, its element the clade (its extent its rows on the rows system, the column as
  `position`); a tree node hits as `tree-node`, its element the clade.
  **Hovered** (ux-foundations section 8: the element is hovered, here and on linked sections): a hover on a tree node
  or a summary row lights the clade: its branches in `seq-hover` ink and the band of its rows (one band, the clade's
  rows being adjacent in tree order; folded, the summary row's band), as a header hover lights a row's. Drawing only.
  **Selected**: a click selects the clade's rows; Cmd/Ctrl+click toggles them; a selection holding some of its rows
  draws a partial mark (the tint over the share, an outline when all). **Exported** as drawn, label and wedge included.
- **What a fold does to the rest**:
  - *the selection* is by data rows: a fold never changes it; folded rows show as the partial mark and on the overview;
  - *the cursor* keeps its data row: folded, it stands on the summary row (announced "in clade first +n"); unfolded, it
    is back on its row; ↑ ↓ step over a summary row as one stop;
  - *a focus or a reveal* of a row inside opens the clades holding it, one jump (fixed, decided by the author, 2026-11-20: "Fix all three"; no
    second use for focusing the summary row instead, which `packages/seq/TODO.md` keeps);
  - *a hover* projected onto a folded row (a structure, a sequence) lights the summary row's cells, partly;
  - *history*: windows stored as row keys, so a back after a fold lands on the same rows (a folded key on its summary
    row, the rule above); the window keeps the row at its top across a fold, so the rows under the pointer do not jump;
  - *an order leaving `tree`*: the folded set is kept but not applied (rows unfolded, no tree); back in tree order it
    applies again;
  - *filtered and hidden rows prune the tree*: a clade's rows are its shown leaves; a clade with none goes, one with
    one row folds to that row; "+n" counts shown rows (hidden ones are in the chip);
  - *a pinned row inside a folded clade* stays in the pinned band and counts in the clade's "+n".
- **Accessibility**: a tree node is announced with its state (`aria-expanded` on the node's accessible item, "clade,
  collapsed, 12 rows"), through the canvas's accessible layer (planned, with accessibility).

**Drawing by density** (ux-foundations section 4), per axis: columns step as any track's; rows step from letter rows
(20 px) to cells to sub-pixel rows. Below a device pixel on either axis the rows draw from **the shared 2D reduction**
(`render/pyramid.ts`), built by the layer in expanded order (only the layer knows the order):

- a sample of a level pair (`[rowLevel, columnLevel]`, factor 4 each, independent: zoomed out on columns a row stays a
  row) holds a majority letter and its share; **a majority is not mergeable**, so each level is the majority of the
  majorities below it (shares averaged): approximate, tested against the exact majority within the measured tolerance
  (decided, 2026-11-19: the share within 0.05 on average and 0.3 at worst; the letter exact where the exact winner
  leads by more than 0.25). Enough because a level's image is drawn below 3 device px a cell, where the shade reads and
  the letter does not; keeping the top two classes per sample was declined (memory per level);
- colored **by the theme's classes** (fixed, decided by the author, 2026-11-20: no `rows.reduce` item; a reduction by letter is in
  `packages/seq/TODO.md`, the tension kept: a theme's plurality against the letters' own): one pass mapping each cell to its theme class through a
  table by letter, or by (column, letter) when the theme declares column dependence (Clustal X's rules), so the
  per-letter shortcut holds both ways; the pyramid then reduces classes, not letters. A theme change rebuilds the
  classes and the pyramid (3M cells: under 50 ms, a target with its test). A lazy source's served row bins are by
  letter (`letter`), said once (`seq.render.rows.reduced`, info);
- built for held alignments at the first layout, up to `maxBlockValues` cells, in about 0.8 of the codes' bytes (2
  bytes a sample); larger and lazy alignments read the source's row levels (data or leaf order only);
- drawn as `image` primitives a tile (256 × 256 samples), never smoothed, coarse under fine while a tile arrives
  (EC-169: never drawn as empty); images in the shared image LRU.

**Lazy rows** (EC-41, EC-44): the layout wants the row tiles covering the rows' window at the row level that fits a
device pixel, at `view`; a tile above and below at `margin`; the last scroll's direction at `prefetch`
([data-store.md](data-store.md), want sets). A row not held answers none and never throws; its label shows `row 5,124`
until its tile arrives. **Bounded allocation**: a frame touches the tiles in view, the labels of the rows in view and
the image LRU; nothing sized by the row count is made per frame (the order's arrays are built at prepare, data order
makes none), so 10M lazy rows lay out in what 10,000 do.

**Selection**: rows are positions of the rows system, held as **a Location part on the rows system** (data rows), so a
reorder keeps it and it projects like any part; a block is a rect (`Location.rects`: a rows part × a columns part).

- a **row** (by its header): its part; Shift+click extends the run in display order from the anchor row (mapped to data
  rows through the order); Cmd/Ctrl+click toggles. `(row :key k)` addresses resolve to this part;
- a **cell**: a click selects the rect (that row × that column), projected onto the row's own sequence where it is
  bound (its residue, or the boundary at a gap: EC-127's hairline);
- a **block**: a drag selects a rect (the rows in display order between the ends, as data rows, × the columns; in
  text, the rows between the ends × the columns between them in reading order, across line blocks); Shift
  extends from the anchor, Cmd/Ctrl+drag adds one, Alt+drag subtracts (rects subtract, above);
- a **clade**: a click on a summary row selects its rows; Cmd/Ctrl+click on a node toggles its rows in the selection
  (a plain click on a node toggles the fold, decided);
- a **column**: a click on the ruler or the logo selects the column range on the columns system; a drag in the
  ruler's upper part (above `navigator.rail.grabShare`, 0.4 of its height, from the bottom) selects a run of columns,
  a drag in its lower part is the ruler rail's and pans ([navigators.md](navigators.md)); a grid's rulers are its axes and pan.

**The keyboard cursor** (EC-154; built, columns and rows (RK, 2026-11-19); text's line blocks not): `{ section, system, position, ySystem, yPosition, part }`, `part` `cell` or `header`.
← from the first column onto the header, → back; ← → step columns; ↑ ↓ step display rows (pinned band to scrolled
rows, a summary row one stop); Page Up and Page Down a window of rows; Cmd/Ctrl+↑ ↓ the first and last row;
Shift+arrows extend the block; Enter focuses the row or cell, or opens a summary row; Space (on key-up, and only when
no drag happened while it was held, so Space+drag pans without toggling) toggles the cell in the
selection, and on the header the row's clade (its node's, or the summary row's). Moves are arithmetic on display
offsets, so a million lazy rows cost what ten do; a cell whose tile is pending is outlined and announced pending. On a
relayout the cursor stays on its row; a folded row puts it on the summary row, a filtered one on the nearest shown.
Space and Enter stay the cursor's keys, acting by its part. ↓ from the pinned band's last row lands on the first scrolled row shown and the band keeps its
scroll (decided 2026-11-19); ↓ from a section's other tracks (its ruler, its logo) enters its rows. Cmd/Ctrl+Shift+
arrows start a new block, the others kept (`cursor.add-up`, `cursor.add-down` for ↑ ↓). The cursor stands on a pinned
row's copy it was moved or clicked to, drawn on both copies. Page Up and Page Down belong to `cursor.page-up` and
`cursor.page-down` (in `view.keys`), which move the cursor a page with the band, or scroll the band when no cursor is
shown; Enter on a header focuses the row's extent (its letters once its tile is held). On a grid ↑ ↓ step its y and
reveal it, as x is. **In text** (the ESPript figure) ↑ ↓ step rows within a
line block and, past its last row, cross to the first row of the next block (the previous block's last row going up).

**Drag past a scrolling box's edge** (the rows band, a text section): a drag-select whose pointer leaves the box
scrolls it toward the pointer, faster the further past (px a second from the distance), timed by the clock (the
`AnimationPort`'s frames carry their time), never by the frame count, so the speed is the same at 30 and 120 Hz; the
selection extends as the content moves under the held pointer. Planned (here for rows, [layouts.md](layouts.md) for text).

**Row zoom keys** (in `view.keys`, as every new command's): `rows.zoom-in` `]` and `rows.zoom-out` `[`, about the
cursor's row, else the band's centre; `view.fit` (`0`) fits the rows with the columns. Why these: a key's chord must
match only its own `e.key` (a printable chord matches whatever Shift types it, so Shift+= arrives as `+`, which is
`view.zoom-in`'s; Shift+− as `_`, Shift+0 as `)`); `[` and `]` are typed alone on the common layouts, never browser
zoom keys (EC-155). The bindings check, which reports conflicts when the table is built, covers `view.keys` too: two
commands whose chords one key event matches are a conflict. A layout that types `[` with AltGr rebinds them there.

**Projection**: a Location on a bound row's sequence (a structure residue, a sequence section's hover) projects onto the
rows as that row's cells; a column range onto every row's cells; a rows part onto its rows. A hover on a rows cell
lights its column and row bands and, through the bound row, its residue elsewhere.

**What a navigator reads** (the navigators design's asks): the rows in display order (count, uniform row height,
folded clades as one row with their summary cells, the pinned band outside the scroll, an order version); the band's
`bandTop`, visible height and `rowScale`; the rows' coarse tiles the rows draw, wanted at `view` priority within the
navigator's set (shared in the store), reduced by the shared pyramid and cached in the shared image LRU; and **one
value a row as the fallback** (`rowsValue`: identity to the reference row, else the non-gap share), computed with the
column measures. The emphasis on rows as data-row parts. The moves go through the section's `nav` (a pan by dy sets
`bandTop`, a zoom on rows sets `rowScale`). Idle when every row fits.

### The logo

- A `logo` layer over the columns, a fixed track above the rows: each column's letters stacked by `scale` (`bits`:
  `log2(k) − H`, less the small-sample correction with `correction`; or `probability`), lowest first; from the
  alignment's profile (the `matrix` shape on the columns: a frequency a letter a column), computed from held cells
  (`gaps`: `ignore` or `count`; `rows`: `all` or `shown`) or served by a lazy source, else its fallback with
  `seq.render.logo.profile`.
- Each letter placed and stretched with the glyph runs' per-letter y and vertical scale (the layouts design's form
  change). Below 2 device px a column: one bar a device column, the tallest stack there, in its top letter's color,
  from the shared reduction's one-axis max. Fallback: `values` (the column's information as bars).
- **A letter's color is its class's, made to read as a glyph** (the author, 2026-11-19: "this letter coloring doesn't
  really work for the logo plot"): a scheme's colors are cell fills, pale by design (letters drawn on them); on the
  track's surface a stacked letter takes its color moved in oklab lightness until it reaches 3:1 against the surface,
  in each theme, the `ink` rule of [themes.md](themes.md) with the surface as the fill. Hue kept, so the logo and the
  cells below it still read as one scheme (WebLogo's and Skylign's letters are saturated for this reason). The bars
  below 2 device px take the same color.
- **And saturated as a glyph** (the author, 2026-11-19, on LG's light logo, its pale fills darkened to greyish tones:
  "Yes, schedule it"): a pale fill has little chroma, so lightness alone gives a muddy letter. The glyph color also
  raises its oklab chroma, hue kept, to a floor (`seq-logo-chroma`, a token per theme) clipped to the sRGB gamut at its
  lightness, then takes the lightness the 3:1 rule needs; a color already above the floor keeps its chroma. A grey
  (a gap, `missing`, a scheme's neutral) stays grey. Why a floor and not the scheme's own strong colors: a scheme
  names one color a class; a second, ink color a class would double every scheme's table.
- Hit: the column, and the letter under the point as the hit's `state`; a click selects the column.

### Data: what this slice adds

The protocol deferred `pairs`, `matrix`, `grid` and `rows` to the Genomes milestone ([protocol.md](protocol.md)).
This slice takes the least of them:

- **`cells`** (with the rows table) on 1D tiles of the rows system, columns whole; coarse row levels served by a
  source as `class` and `share` per row bin and column. The store's keys, tiles and want sets are unchanged.
- **`pairs`** and **`grid`** as inline annotations, read whole through the `seq` format as inline data is; their
  pyramids and images are the layers'.
- **`matrix`** for the profile, on the columns system, computed or served.
- **Deferred to the Genomes milestone**: a request's second axis (`RegionRequest.y`), 2D tiles, the store's grid
  summaries, pairs and grids read by region (EC-208, EC-223).

### Performance and size (each target with its path and test)

- **The 10,000 × 300 alignment**: first frame under 300 ms (`attach`, the measures and profile, the theme classes, the
  pyramid, the first layout: about 10 ms each pass over 3M cells, so well inside; heavy, `itWithin`); a pan or zoom
  frame under 8 ms at the fit, at letters and between two level pairs; a cursor step and a block drag step under
  4 ms; a theme change's rebuild under 50 ms; folding every clade under 16 ms. Quick representative: one fit
  frame and one letters frame bounded.
- **Scaling**: frame time at 10,000 and 100,000 held rows and 10M lazy (synthetic) within 1.2× (heavy, ratio); layout
  allocation equal at 10,000 and 100,000 rows (the opt-in tier); draw and hover 0 bytes warm.
- **A dot plot's frame** (its quick representative): a 2,000 × 2,000 PAE's pan frame under 8 ms once its tiles are built, a
  tile build within `grid.buildBudget` a frame. **Arcs**: 4,096 arcs a frame under 8 ms; a hover over nested arcs under
  2 ms.
- **Size, minified, estimated** (the author: "don't overdo it on the limits ... but obviously try to minimize the
  amount of code"; no size cuts a feature): pairs and grids about 12 kB in headless rendering (the grid layout 3, the grid layer 2,
  the pairs layer with split and hit 6, emphasis along y 1); rows about 17 kB (`cells` and the rows table 2, the rows
  layer and axis 8, the tree and folds 3, the logo 2, the column measures 2); the shared pyramid and image cache about
  2 kB once; the readers about 6 kB in their own entry (aligned FASTA, Stockholm, Newick), and about 6 kB in the viewer
  (commands, cursor, bindings). Measured when they land.

## The interaction matrix's cells

Cells this design fills (ux-foundations section 8); the sub-rows name the new things acted on.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element: a row** (header) | click; `selection.set` | Shift+click: the run in display order | Cmd/Ctrl+click toggles | double-click: the row shown, x its residues; `address.focus` | Alt+Left | — | — | `rows.hide` (patch); shown again by the chip | pin: `rows.pin` (patch, the header's icon); no: residues are not edited |
| **an element: a pair** (arc, chord, dot) | click: both ends | Shift+click: the range from the anchor end | Cmd/Ctrl+click | double-click: both ends' hull | Alt+Left | — | — | no: hide the layer or filter by kind (a data filter, [themes.md](themes.md)) | no: data |
| **an element: a cell** (grid or rows) | click: a one-cell rect | Shift+click: the rect from the anchor | Cmd/Ctrl+click | double-click: a grid's cell zooms the plot a step about it (decided 2026-11-19); a rows cell, both positions | Alt+Left | — | — | no | no |
| **a range: a block or rectangle** | drag on cells | Shift+drag | Cmd/Ctrl+drag adds, Alt+drag subtracts (exact: rects whose sides hold several ranges) | Enter: x and y to the rect | Alt+Left | — | — | — | the selection's bookmark command (planned) |
| **a range: columns** | click or drag on the ruler or logo | Shift+click | Cmd/Ctrl+click, Alt+drag removes | Enter | Alt+Left | — | — | no: hiding columns waits (decided) | — |
| **the selection** (rows, rects) | — | Shift+arrows: the block | Cmd/Ctrl+Shift+arrows start a block | Enter | Alt+Left | — | — | — | the selection's bookmark command (planned): a ref of its Location |
| **a leading section: rows** | — | — | — | — | — | rows: Cmd/Ctrl+wheel over the header, `]` and `[` (`rows.zoom-in`, `rows.zoom-out`), `0` fits rows and columns | the wheel's vertical delta over the band, drag on the header, Page keys | `rows.filter`, `rows.hide` | order, filter, pins by patch; `collapsed` by patch, `clades.keep` |
| **a leading section: grid** | — | — | — | — | Alt+Left | Cmd/Ctrl+wheel, pinch: both axes about the pointer (y alone on y's ruler when y has its own window); -, = | wheel: dx x, dy y; drag on a ruler pans that axis (decided: a grid's axes pan); Space+drag both (Space's toggle acts on key-up only when no drag happened) | no: hide its layers | `aspect`, `triangle` by patch |
| **a follower: a grid's y** | as a leading one | | | | | moves its leader (`lead`) | moves its leader | | |
| **the view** | — | — | — | 0 fits, a grid on both axes | the history holds the rows' windows and the folds a jump changed | — | — | — | — |
| **a track or layer: the rows layer** | click a header: a row | — | — | — | — | — | — | `layer.hidden` | `order`, `filter` by patch, previews while chosen |
| **a clade** (a `tree-node`, a `summary-row`) | click a summary row: its rows | Shift+click: the run of display rows to it | Cmd/Ctrl+click a node or summary row toggles its rows | double-click or Enter on a summary row: opened and shown; `clades.reveal` opens the clades holding a row | double-click again, Alt+Left | — | — | fold (session state): click on its node, Space on key-up with no drag while held (`clades.toggle`), inverse the same; `clades.collapseAll` / `clades.expandAll`; `clades.collapseTo` (depth, cut), inverse `clades.expandAll` or Alt+Left | the opening set: `clades.keep` (one patch) |
| **a bookmark or link** | — | — | a rect or a rows part as a ref | runs it: x and the rows (opening folded clades) | Alt+Left | — | — | — | its label (planned) |
| **a navigator: the vertical overview** ([navigators.md](navigators.md)) | — | — | — | a click centers the rows there | Alt+Left | drag the `window` target's `start` or `end` part: rows | drag its `body` | the host hides it | — |

No cell is the hover's: hover is a method, not a verb ([viewer-models.md](viewer-models.md)); on a grid it is the
crosshair, on rows the row and column bands, on a pair both ends.

## Configuration

Whose, by the rule the axis design states (spec param, `RenderConfig` preset, `SpecBounds`, `ViewerConfig`, a port, a
registry entry, an option of a call, or fixed with its reason).

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `side`, `reach`, `height`, `crossing` (pairs layer) | spec param | `above`, 0.25, `span`, `same` | arc diagrams' convention; an off-window end moves continuously (EC-129); nested arcs nest; crossing shown as crossing | `below`; any fraction; a size theme; `split` (R-chie) |
| `acrossLines` (pairs layer, text) | spec param | false | arcs across many lines cross the text | true |
| `mirror` (pairs layer, grid) | spec param | false | a pair drawn where its data puts it; the tRNA builder sets it | true |
| `pairs.maxArcs`, `pairs.hitPx` | `RenderConfig` preset | 4,096; 4 | bounded by the screen; a curve is thin | any |
| `aspect`, `follow_y` (grid section) | spec param | 1; `x` | square cells; y shows x's window | a ratio, `none`; a section, none |
| `triangle` (grid layers) | spec param | `both` | the whole matrix | `upper`, `lower` |
| `aggregate` (grid layer) | spec param | `mean` | a sample reads as its cells | `max`, `min` (PAE: confident pairs) |
| `grid.rectMinPx`, `grid.buildBudget` | `RenderConfig` preset | 12; 65,536 samples a frame | rects only where few; frames under 8 ms | — |
| `imageBytes` (grid, rows and navigator tiles) | `RenderConfig` preset | 64 MB | one bounded cache | — |
| `height`, `summary`, `numbers`, `treeWidth` (rows layer) | spec param | `fill`, `consensus`, `gutters`, 120 px | the rows take the room; one row reads as its clade; numbers per row (EC-39); a readable tree | a number of rows; `logo`, `bar`; `none`; 0 hides the tree |
| `collapsed`, `pinned`, `hidden`, `order`, `filter` (rows layer) | spec param | none; `data` | the writer's opening state and order | — |
| `rows.pinInPlace` | `RenderConfig` preset | true | the tree stays a tree | false: moved only |
| `rows.labelMinPx` | `RenderConfig` preset | 8 | a label needs a row | — |
| `scale`, `correction`, `gaps`, `rows` (logo layer) | spec param | `bits`, `none`, `ignore`, `all` | the logo's usual reading; exact for given frequencies; gaps are no letter; the alignment's logo, not the view's | `probability`, `small-sample`, `count`, `shown` |
| a click on a tree node | `ViewerConfig` (the bindings' preset) | toggles the fold (decided) | a node is chrome, as a group's chevron | a preset binds it to select |
| keys of `rows.zoom-in`, `rows.zoom-out` | `ViewerConfig` (`view.keys`) | `]`, `[` | each matches only its own `e.key` (above) | any key the host's layout types alone |
| the folds an export draws | option of `renderSvg` (`folds`) | the layer's `collapsed` | the viewer passes its own | — |
| folds as session state; a fold's motion (`view.motion`, `view.jumpMs`); split exact to 2,000 pairs | fixed | — | decided; the jump's motion reused; no second use | — |

## In the specs

| Decision | Where |
|---|---|
| Clade folds session state, `collapsed` the opening set, `clades.keep` a patch; an edit of `collapsed` resets them (Decided 1) | viewer › The rows' commands |
| A click on a tree node toggles its fold; a summary row selects (Decided 2) | viewer › A tree's nodes |
| A reveal opens folded clades, back folds them (Decided 3) | viewer › The rows' commands |
| Hiding alignment columns waits (Decided 4) | design-only: not built (decided to wait) |
| The selection's bookmark command (Decided 5) | viewer › Configuration |
| `mirror` a plain `false` (Decided 6) | core › The pairs layer |
| A drag on a ruler: the lower part the rail's, the upper selects columns; a grid's axes pan (Decided 7) | viewer › A ruler rail |
| Figure-meaning settings are spec params; `selection.copyLetters` copies FASTA (Decided 8) | viewer › Configuration |
| The grid a section layout with a second transform (`strip.y`, `ctx.axisY`) | core › One strip, the plot |
| The rows' zoom limits the host's (`rows.zoom`) | viewer › Configuration |
| The rows layer owns its y inside its band (`bandTop`, `rowScale`) | core › The transform |
| An order under the rows' pieces | core › Display rows |
| Hits stay flat; `ySystem`, `yPosition`; the targets `row`, `tree-node`, `summary-row` | core › The rows layer |
| Emphasis along y | core › Emphasis along y |
| The rows' window kept by a row key, restored by position | viewer › Bands |
| Per strip by layout: rows in a text section, every clade unfolded | core › In a text section |
| Pairs: data, either end, self-pairs, an end never wrapping | core › Kinds |
| Pairs laid out by section through `ctx.place`; arcs, `reach`, `height`; on a ring a chord; in text `across_lines` | core › Arcs |
| An end on split copies placed by unit; never on a guessed copy (`seq.render.pairs.copy`) | core › Ends on blocks |
| Crossing pairs (`same`, `split` exact to 2,000) | core › Crossing pairs |
| Kinds at one cell sliced (EC-63) | core › The pairs layer |
| A glyph per kind with "+" | design-only: not built (the package's TODO) |
| Pairs bounded by the screen (`pairs.maxArcs`); copies fold by group; the hit by hull and distance | core › The pairs layer |
| The grid section: `of_y`, `follow_y` | core › The grid layout |
| A grid layer draws only on its own systems | core › Only on its own systems |
| A grid over a complex (blocks on both axes), gaps, tinted gaps in a dot plot | core › The grid layout |
| A summary of a grid over a rect | core › A summary over a rect |
| Layers in it: `grid`, `pairs` as dots; `mirror`, `triangle`; later over earlier | core › Triangles |
| A plot's box fixed | core › A plot's box is fixed |
| A double-click on a cell zooms a step; the wheel over a plot is the plot's | viewer › Wheel |
| Size and aspect (EC-109) | core › The grid layout |
| Grid data inline, bounded by `maxBlockValues` | core › Only on its own systems |
| Drawn: rects at `grid.rectMinPx`, image tiles from the shared reduction, `grid.buildBudget` | core › The grid layer |
| The crosshair | core › Emphasis along y |
| Rectangles subtract exactly | core › rects subtract |
| One 2D reduction and one image LRU (`imageBytes`) | core › Render |
| The rows layer's data: `cells` on 1D tiles, columns whole; readers; the tree | core › `cells` |
| Column measures for themes (`share`, `identity`, `matches-reference`) | core › Themes |
| Order (`data`, `tree`, a field); filter; hidden | core › The rows layer |
| A lazy alignment in tree order, a lazy clade's summary from row bins | design-only: not built as designed (drift note `grids-and-rows`) |
| The rows axis by arithmetic | core › Display rows |
| Clades named by a unique non-support label, else first and last | core › Clades |
| Pinned rows in a band and in place | core › Display rows |
| The band: the rows take the section's room left; one rows layer a section | core › One rows layer a section |
| The row header: tree, label, numbers of the residues in view, `rows.labelMinPx`, the hidden chip | core › The header |
| Row columns (`columns`), their fields, titles ordering, a cut value whole under the pointer | core › Row columns |
| A narrow frame gives way in one order | core › A narrow frame |
| Narrow: columns yield, a chip naming them | core › Narrow |
| Gestures and commands of the clades | viewer › The rows' commands |
| The tree's look | core › The tree's look |
| The summary row, its hit, hover and partial selection mark | core › The rows layer |
| What a fold does to the selection, cursor, focus, history | viewer › The rows' commands |
| Clades announced with their state (`aria-expanded`) | design-only: planned with accessibility |
| Drawing by density: the shared reduction by theme classes, majority of majorities within tolerance | core › The rows layer |
| Lazy rows: the rows' window wanted, bounded allocation | core › The rows layer |
| Selection: rows as a Location part on the rows system; blocks as rects | viewer › A cell of two axes |
| The keyboard cursor on rows | viewer › The cursor |
| Drag past a scrolling box's edge | viewer › Pointer |
| Row zoom keys `]` and `[`; the keys check | viewer › The keys check |
| Projection onto the rows | viewer › Projection |
| What a navigator reads: the rows track itself | core › The picture |
| One value a row as a navigator's fallback (`rowsValue`) | design-only: not built as designed (drift note `grids-and-rows`) |
| The logo: scale, correction, profile, bars below 2 device px, glyph colors with the chroma floor; hit | core › The logo |
| Data this slice added (`cells`, inline `pairs` and `grid`, the profile) | core › Cells |
| Performance and size targets | design-only: measures, held by the tests and the performance tiers |
| Precedents, the interaction matrix's cells | design-only: intent |

## Tensions

- **Folds as session state against spec edits**: a figure needs the folds saved; an undo step per fold buries the
  edits; `track.hidden`, a toggle of the same family, is a patch (ux-foundations section 8). Decided: session state,
  `collapsed` the opening state, `clades.keep` writes it; pins, hidden rows, order and filter stay spec edits. No
  `spec` mode (each fold a patch): no second use yet. Standing: a fold is session state while hiding a track is an edit.
- **A click on a tree node toggles, not selects**: the previous implementation's gesture, against "click selects"
  (Jalview selects the group). Decided: the toggle (a node is chrome, as a group's chevron); a preset may bind it to
  select.
- **Reductions by the theme's colors or by letter**: colors are what is drawn, letters what a lazy source serves.
  Colors where cells are held, letters otherwise, said once; open: a served class histogram per bin would close it.
- **A pyramid of majorities is approximate**; counts per letter a sample would be exact at about twenty times the
  bytes. Majorities, tested within a tolerance; clade summaries exact.
- **Images in SVG**: an `image` is runs of rects, so a zoomed-out 10,000-row export is large; a PNG inside SVG needs an
  encoder in core. Rects in this slice; measured with the PNG export, which decides.
- **Lazy order**: tree order only when a source serves leaf order; a field order on a lazy alignment is refused.
- **No navigator for a grid's own y window**: a grid whose y has its own window (not following a section) has no
  navigator over y; answered by [genome-coordinates](genome-coordinates.md)'s square (two axes over blocks × blocks).
- **Arcs between sections** (ribbons, EC-209): closed by [across-sections](across-sections.md) (a `bridge`, the `pairs`
  layer drawing ribbons in it; a grid's `follow_y` a second section of blocks).
- **Emphasis along y** is the two-axis layers' own, through the layouts design's y path (answered there: `axes: 2`).

## Self-check

The claim: **this design covers what a user does with pairs, grids and rows.**

1. **Assume it holds**: a user reads arcs and a dot plot, scrolls and zooms an alignment's rows, folds clades, pins and
   selects rows and blocks, steps with the keyboard, reads a logo, in tracks and in a text figure; each has a gesture,
   a command, a bounded draw and a test.
2. **What it rests on, negated**:
   - *rows are y and columns x* → the symmetric case: positions down and records across;
   - *one alignment a section* → a plural: two alignments of one family, scrolled together;
   - *the user acts in the rows* → the inverse: acting in a structure or a sequence and expecting the rows to answer;
   - *hide, pin, fold, filter* → their inverses: showing what was hidden, finding what was filtered;
   - *a clade can be opened or closed at will*, negated both ways: **always open** (no tree, a lazy source in data
     order, a reader who never folds) and **always closed** (a figure saved folded, a depth cut, then a search or a
     structure pointing inside);
   - *folds are the reader's* → a fold that must be kept (a figure, a story);
   - *a section is one straight strip* → the strip world: an alignment wrapped in text, pairs between lines, a chord on
     a ring;
   - *pairs on one system or section* → pairs between sections and grid axes that are not positions.
3. **The negated world**: a residue selected in the structure shows nothing in the alignment; three hidden rows are
   lost; two alignments scroll apart; a PAE's rectangle must show in the tracks and the structure; always open, the
   tree's groups are invisible in 10,000 rows; always closed, a row asked for is unreachable, a selection inside
   invisible, the cursor homeless, a folded set lost when the figure is saved; in text, a fold would differ per line
   and an arc from line 3 to line 9 would cross six lines of letters; samples × samples is asked of a grid.
4. **Both held**:
   - *Integrated*: projection onto the rows from a bound row's sequence; the hidden-rows chip; a rect lighting both
     sides on linked sections and the crosshair drawn for every channel, so a selection in the tracks shows on the grid;
     rect subtraction; one flat hit pair for rows and grid cells; reveal opening folds as a jump that back undoes; the
     partial mark; the cursor on the summary row and on the header; the bulk commands; `clades.keep` and the
     story step's folds; a history step restored by row key across a fold, a filter or a reorder, never dropped; folds kept but not applied off tree order; in text every row unfolded on every line and
     `acrossLines` off by default; arcs and chords by `ctx.place`, so the strip world needs no layout-specific
     layer code.
   - *Standing*: rows always run down; ribbons and categorical grid axes closed by [across-sections](across-sections.md) (EC-209,
     EC-223); a lazy alignment in data order has no clades to fold; a fold is session state while hiding a track is
     an edit; two alignments share no row window (a scroll group over rows axes waits for a dataset that needs it).

## Consequences

**Order of work**: the rows land the shared changes first, pairs and grids build on them.

- **Rows** (first): `Hit`'s and the hover record's `ySystem`, `yPosition` and the targets (`cell` with them, `row`,
  `tree-node`, `summary-row`); the shared reduction (`render/pyramid.ts`) and image LRU on core's `LRUCache`, reused by
  the grid and the navigators; the rows layer reading the layouts design's `bandTop` and `rowScale`; the `cells` shape and the rows table on 1D row tiles; the column
  measures `share`, `identity`, `matches-reference` registered for themes; the aligned FASTA, Stockholm (`#=GC RF`,
  `SS_cons` as pairs) and Newick readers, with fuzz tests (malformed input reported, never thrown); then the
  `alignment` and `tree` kinds, the rows layer (order, folds, pins, header, reduction, lazy rows), the logo and the
  profile, the clade commands, the cursor's `part` and `y`, the selection's rows part.
- **Pairs and grids** (on the rows): the `grid` section layout with `strip.y` and `ctx.axisY`; the `pairs` annotation and layer (section
  scope, arcs, chords, dots, split, kinds, hit); the `grid` annotation and layer (the shared reduction and images,
  triangles); the
  crosshair and the rectangle, rect subtraction in the Location.
- **`seqx-core`**: language (`alignment`, `tree`, `pairs`, `grid`, layers `rows`, `logo`, `pairs`, `grid`, a section's
  `layout: grid` with `of_y`, `follow_y`, `aspect`; the version moves, pre-release); data (`cells`, the rows table,
  `matrix`); resolve (rows as a system, bound rows, the tree matched); render (as above); query (`(row :key k)`,
  `(clade :of [a b])`, rects over rows × columns, rects subtract).
- **`seqx-viewer`**: `rowScale` and `bandTop` moved through the section's `nav`; **the history's record** gains, per
  rows band, the top row's key and `rowScale`, and the folded set when a jump changed it; the folded set per rows layer
  in the `ViewModel` (session state, read and set as addresses by the `AddressModel`), passed to `renderSvg`'s `folds`;
  the hit's targets; the selection's rows part and rects; the cursor's
  `part`, `ySystem`, `yPosition`; the commands `rows.pin` (toggles: a pinned row is unpinned), `rows.hide`,
  `rows.show` (the hidden-rows chip's), `rows.filter`, `rows.order`, `clades.toggle`,
  `clades.collapseAll`, `clades.expandAll`, `clades.collapseTo`, `clades.reveal`, `clades.keep`, `rows.zoom-in`, `rows.zoom-out` (keys in `view.keys`); the bindings.
  Tests for the folds: selection, cursor, hover and history across a fold and an unfold, reveal and its back, the bulk
  commands and their inverses, an edit of `collapsed` resetting the set, a story step's folds, an export's set.
- **`seqx-examples`**: tRNA (arcs nested and crossing, the base-pair dot plot), Globins (family with logo and tree,
  folded clades, the ESPript figure's rows in text, the 10,000-row view), the Complex's PAE and contacts; drawing
  checked against the data at the fit and across a tile edge (every pair in view has its arc and every arc is hit;
  every row in the window has its band; every cell at letters is hit as its row and column).
- **Storybook**: those views, canvas and SVG side by side, each with a headless twin on a fake clock (the lazy 10M-row
  story on the synthetic async source); a Design page for this design.
- **Tests**: property tests over random alignments, orders and folds (every row placed once, the hit the inverse of
  placement, a block's data rows those drawn between its ends, a reorder or a fold keeping the selection), the pyramid
  against the exact majority within its tolerance, over random pairs (nesting, every arc in its band, a hit on the
  curve), over random grids and rects (the crosshair at the hovered cell, a rect minus a rect covering exactly the
  difference); reused buffers against fresh ones for the shared reduction and image LRU; the bindings check reporting
  a `view.keys` chord another command's matches (Shift+= against `view.zoom-in`'s `+`); records for the tRNA arcs and the Globins family;
  the targets above.
- **`_spec` sections**: core Data (`cells`, the rows table, `matrix`), Render (`strip.y`, the rows, logo, pairs and
  grid layers, the grid section, emphasis along y), Queries (row, clade, rect subtraction); viewer ViewModel (the
  history's record, `rowScale`, the band's scroll), HoverModel and SelectionModel (the flat y fields, rows parts,
  rects), Bindings, Commands, Configuration.
- **Edge cases with an outcome here**: 10, 35, 38–44, 58, 59, 62, 63, 100, 103, 106, 107, 109, 129, 145, 146, 147,
  154, 155, 206 kept; EC-208, EC-209 and EC-223 still deferred.
