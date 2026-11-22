# Layers across sections: bridges, ribbons, grids over two sections, categorical axes

Lane LX's design (the fifth slice, stage 0; agreed by the author 2026-11-21 and promoted from its draft in `execution/fifth-slice/drafts/`). It amends
[grids-and-rows](grids-and-rows.md) (pairs in a bridge, a grid over two sections of blocks, categorical
axes, its Tensions and its self-check's standing residue), [layouts](layouts.md) (the view stacks
sections and bridges; a third layer scope), [maps](maps.md) (a map's runs as a layer's data),
[viewer-models](viewer-models.md) (a lock beside follows and scroll groups),
[ux-foundations](ux-foundations.md) (the matrix's bridge rows) and
[extension-openings](extension-openings.md) (the extension entry) (below, Consequences). Built on Gate
0's answers as recommended (the fifth slice's plan): the cut order puts contact maps before the liftover (1), so the
liftover figure is the part that goes last; public API changes stay in `packages/seq` (no `@molstar/core` or
`@molstar/io` change). Built on the three designs before it and taking what they decided:
[genome-coordinates](genome-coordinates.md) (GA: assemblies, blocks, `placesOf`, runs of small blocks, **a ribbon's
twist**, **hit tests per drawn run**, follow between unequal block lists refused, stubs where no block shows a
position, the `block` field), [region-reads](region-reads.md) (RR: pairs tiled in 1D by either end, the pairs level of
bin-pair counts with a strand, items gated by counts with hysteresis, **contact-map tiles** in 2D, the whole-assembly
level, the joint small-block rule, byte bounds, sided pairs sources) and [genome-formats](genome-formats.md) (GF: the
`seq` runtime's sided pairs found from either end, a map's runs through `map.source` read by either side, side `a` the
from side, the language additions riding GC1's bump, the lane order RG → FR → LA). The synteny, dot plot and liftover
views read their alignment from the `seq` format (PAF and chain runtimes are designed only, GF).

## Intent

The JBrowse reply's fourth ask, "ribbons joining two sections, a dot plot with its own y section" (EC-209), deferred
until now: [features](features.md) 151 ("Layers across sections … trigger: the Genomes milestone") and
152 ("Contact maps … trigger: layers across sections"); [grids-and-rows](grids-and-rows.md) Tensions
("Arcs between sections (ribbons, EC-209) stay deferred") and its self-check's standing residue ("ribbons and
categorical grid axes wait for the Genomes milestone"); [version-1](version-1.md), Deferred. The
milestone's views that need it ([examples](examples.md), Genomes): synteny ribbons with an inversion, a
translocation and a gap (EC-206); a dot plot of one genome against the other; a liftover through an alignment's map (EC-212); a
contact map by region (EC-223). What a user does with them: compare two genomes (where does this chromosome go, what is
inverted, what is missing), carry a locus from one assembly to the other, read contacts between two chromosomes, and
read a value between categories (samples × samples) on a grid.

What the agreed designs give already: one section's pairs as arcs or chords (`scope: 'section'`), a grid section whose
y follows another section (EC-209's second half, built), projection of hover and selection along the view's links.
What nothing gives: a place for a mark that belongs to two sections, a layout and hit test for it, how the two sides
move together, a grid whose two axes are two different sections of blocks, and a grid axis that is not positions.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| JBrowse 2's linear synteny view: two (or more) linear views stacked, a "level" between each consecutive pair holding synteny tracks; each feature a filled band (a quadrilateral with curved sides) from its query extent on the upper view to its target extent on the lower; drawn at the true x even off screen; colour and alpha by strand or a field; a minimum width; "link views" (pan and zoom both by the same delta); "square" (equal zoom); click a feature for its card | the level as a `bridge` between neighbours in the stack; a ribbon a pair, its sides cubics with vertical tangents; ends at their true x, off-screen ends leaving sideways; a minimum width; `lock: px` (link views); a ribbon's card; levels chained for three genomes | a synteny view as its own view type (here a bridge in any tracks view); per-level renderer settings beyond the layer's params |
| JBrowse 2's dot plot view: two assemblies on x and y, each a list of displayed regions, PAF lines as segments, a diagonal by strand | a grid section whose x and y follow two sections of blocks (GA), the built `pairs` dots and diagonal segments, RR's bins as dots | the dot plot's own region editor (here the sections' blocks) |
| Circos: links and ribbons between ideogram positions inside one circle, twisted ribbons, colour by chromosome, a ribbon's thickness from its ends | twist (GA's rule), colour by the a end's block (EC-218), ribbons wide by their ends | ribbons between two circles: Circos draws every link inside one circle; so does every ring precedent (Rondo, pyCirclize): a ring's links stay the section's (`pairs` chords; RN's `rondo:ribbons`) |
| UCSC chain and net tracks: a chain drawn in one assembly's coordinates as boxes (aligned blocks) joined by lines (gaps), coloured by the other assembly's chromosome; nets as nested fills; liftOver with `minMatch` (0.95) refusing a location that maps less | the chain track as the `pairs` layer's `ranges` fallback on one side, coloured by `b_system` (RR's column) over GA's block-ordered domain; liftover as projection through the chain's map (GA, GF), a mapped share for a host's readout deferred (Liftover) | nets (levels of fills and gaps); `minMatch` as a refusal |
| Juicebox: one map, x and y each a chromosome chosen from a list ("chr1 vs chr4"), "All" for the genome, 1D tracks along the top and down the left, a control map in the other triangle | x and y following two sections (a chromosome each, or an assembly each), "All" as RR's whole-assembly level, the control in the other triangle (EC-59, built) | tracks down the left of a map (a section beside another: [layouts](layouts.md) "side by side later"; a Gap) |
| HiGlass: views side by side with zoom and location locks, a viewport projection (one view's region drawn as a rectangle on another), 2D annotations | the lock (`lock: px`, for two genomes; a follow for one); a projection is a navigator's window (navigators.md) | views side by side (stacked here) |
| This codebase: sections stacked, `scope: 'section'` layers placed by `ctx.place`, a grid section's `follow` and `follow_y`, a grid over a section of blocks (one leader for both axes), the rows layer's labels and `order`, projection along links, `via` | a bridge laid out after its two sides from their placements; `follow_y` naming a second section of blocks; a categorical axis drawing the rows layer's labels in the rows layer's `order` | a fourth section layout: a bridge has no axis of its own, so it is not a section |

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| hover a ribbon: both ends lit, a card (JBrowse) | built: the pair's element, its ends lit on both sides, the card (Decision, Hover) |
| click a ribbon: select it, open its card (JBrowse) | built: the pair selected (both ends); its card is the hover card pinned by the host |
| double-click a ribbon: go there | built: each side focuses its end, one jump |
| pan or zoom one genome alone (JBrowse, default) | built: each side is a leading section of its own |
| link views: pan and zoom both by one delta (JBrowse, HiGlass's locks) | built: the bridge's `lock: px` (patch), or a drag or wheel over the bridge (moves both, whatever the lock) |
| bring the other genome to what this one shows ("square", "center on") | built: `bridge.align` (Shift+A at a side, the bridge's menu), a jump |
| hover a position on one genome: its partner lit on the other | built: through the map's path where the data is a map (projection, as built); ribbons with an end at the hover lit, their far ends with them, where it is only pairs |
| select a region on one genome, see and select its partners (JBrowse's "select in other view") | built: the selection's ribbons drawn in the accent; `bridge.selectAcross` adds their far ends |
| colour by strand, by identity, by chromosome (JBrowse, Circos) | built: the layer's colour channel over `strand`, a field (`identity`, `mapq`), `block`, `b_system` |
| filter short or low-quality alignments (JBrowse's min length, MAPQ) | built: the data filter's numeric form (EC-225, the query sublanguage) |
| three genomes in a chain (JBrowse's levels) | built: A, bridge, B, bridge, C |
| A against C across B (non-neighbours) | refused: a ribbon would cross B's tracks; reorder the sections |
| a dot plot, zoom into a block pair (JBrowse, Juicebox) | built: the grid section's gestures, the square navigator (GA) |
| choose a chromosome for x and one for y (Juicebox) | built: two sections, each a chromosome; a host's lists patch their `of` |
| tracks down the left of a contact map (Juicebox, HiGlass) | refused this slice: sections stack; side by side is layouts' later (Decided, 13) |
| liftOver a region (UCSC) | built: projection through the alignment's map, every image lit; a readout of the share mapped deferred (Decided, 15) |

## Decision

### Whose each choice is

By [axis](axis.md)'s rule. No new port: a bridge is laid out from its sides' placements in core, its
data read through the store as RR's pairs are. What the figure says (which sections a bridge joins, its height, whether
its sides move together) is a spec param; how ribbons are drawn and bounded is a `RenderConfig` preset; what a gesture
over a bridge does is the viewer's.

### A bridge: what a layer between two sections is in the language

- **A new kind, `bridge`, a child of `view`** beside its sections and navigators: `a` and `b` (refs to two sections of
  the same view, required, distinct: `seq.lang.bridge.sides`), `label`, `height` (CSS px; absent: `bridges.height`,
  96), `lock` (`none`, `px`; absent: `none`), `hidden` (false). Its children are `layer`s (no tracks: a bridge has one
  band, the room between its sides, and its layers composite in it, later over earlier). Why a kind and not a section:
  every section has an axis, a viewport, a follow, a history entry and a cursor; a bridge has none of them (its two
  axes are its sides'), so a section with `layout: bridge` would make every reader of sections skip one. Why under
  `view`: the navigators are there already, view-level marks bound to sections by ref.
- **Neighbours in the stack**: the view stacks its children in order ([layouts](layouts.md), "The view
  stacks its sections"); a bridge's `a` and `b` are the section just above it and the one just below, in either order.
  Otherwise `seq.view.bridge.neighbours` (error at the bridge) and it is not drawn, its room not taken. Why: a ribbon
  to a section beyond a neighbour crosses that neighbour's tracks; three genomes are two bridges (A, A–B, B, B–C, C).
- **Sides are rows**: each side is a section in the `tracks` layout, its blocks (if any) `arrange: row`. A ring, a text
  section, a grid, or blocks in a stack as a side is `seq.view.bridge.layout` (error), not drawn. Why: a ribbon leaves
  a side at its edge facing the bridge; a ring's rim faces the bridge on one half only, a text section's line holding a
  position is anywhere in its box, a stack's block anywhere in its rows, so the ribbon would cross the side's own
  tracks. Links inside one ring stay that ring's (`pairs` chords, RN's ribbons, GA's two assemblies on one ring,
  EC-210): every ring precedent draws them so.
- **Its layers** are `pairs` layers (below). Another type is `seq.view.bridge.layer` (warning) and drawn as its
  fallback chain's first `pairs`, else not drawn. Why: a bridge has no axis, so a type that draws along one has nowhere
  to draw. The layer contract's `scope: 'bridge'` is internal this slice, not in the extension entry: no extension
  draws in a bridge yet (RN's ribbons are a ring's), and the entry is what a review must freeze (Decided, 16).
- **Its data**: a `pairs` node (inline, or `pairs.source` read by region, RR), or **a `map` with `method: segments`**
  (inline `segments`, or `map.source`, GF) read as pairs: a run is a pair of equal length, `a` its `from` range,
  `b` its `to` range, `b_strand` its strand. So one alignment (a `seq`-format pairs source of runs, a host's chain or
  PAF adapter alike) is drawn and projected through from one node: the ribbons and the liftover. Amends [maps](maps.md): a segments
  map may be a `pairs` layer's data, anywhere (in tracks its runs draw as arcs). A map of another method is
  `seq.view.layer.shape`.
- **Which systems**: each end is placed on the side whose section shows its system (a system of the side, a block of
  it, GA's `placesOf`); when both sides show it (one genome against itself, an unsided BEDPE over one assembly), the
  data's `a` end on the bridge's `a` side, `b` on `b` (with `mirror`, built, also the swap). A sided source (RR's sides,
  served by the `seq` runtime: `chr1` on each side, distinct) is placed by its side's assembly, so a bridge whose `a` is grape's section and whose
  data's `a` is peach draws the same ribbons as the reverse. **No path**: pairs are not projected
  ([grids-and-rows](grids-and-rows.md), Pairs), so an end on a system neither side shows is not drawn,
  counted in `seq.render.bridge.unplaced` (info, once a layer, with the first unknown systems). Why: a ribbon through a
  composed path would draw a claim the data does not make.
- **A reader without bridges** (an older minor, or another renderer, EC-235, EC-238 to EC-244) leaves the bridge out:
  its sides stack as sections, each still drawing its own layers. Said once (`seq.view.bridge`, info, from the kind's
  unknown-kind rule). Why: the sides are whole figures; the bridge adds a relation, never their content.

### Layout: the bridge's band and its ends

- **Placed after its sides**: `layoutView` places each section (as built), then each bridge in the box between its
  sides (its `height`, plus the sections' gap above and below). A bridge's **placement** (`BridgePlacement`) holds its
  box and its two sides (`BridgeSide`: the side's `Placement`, the y of its edge facing the bridge (the side's box's
  bottom when above, top when below), and `placesOf` over its blocks).
- **Lazy, all but the language and the resolved view** (the plan's risk: the default viewer at its bound): the bridge
  layout, the ribbons, and the viewer's bridge parts (the band's pointer, the lock's carrying, `bridge.align`,
  `bridge.selectAcross`, the cursor across a band, the resize grip) load as one chunk when a resolved view has a
  bridge, the viewer's parts as a `BehaviorProvider` registered then. Eager: the `bridge` kind in the language,
  `ResolvedView.bridges` and `.stack` (resolve reads them; a view without a bridge has an empty list) and `Hit.bridge`
  (one number). Measured at LA's gate against the default viewer's bound; past it, the resolved parts move into the
  chunk with the view module's bridge pass, and the author decides.
- **Height by drag** (JBrowse's levels): the band's lower edge is a grip (4 CSS px, `row-resize`); a drag previews and
  ends in one patch of `height` (one undo step), Esc cancels. Why: synteny readers resize the band to unflatten long
  ribbons.
- **The key**: a bridge's key is its own params and both sides' keys (viewport, the key's reads); a pan of one side
  relays that side and the bridge only. Why: a ribbon's ends are wherever its sides draw them.
- **A layer of scope `bridge`** (a third value of `LayerType.scope`) is laid out once a bridge with a `BridgeContext`:
  `ctx.box` the bridge's band, `ctx.sides` the two `BridgeSide`s (each end's page x from `placesOf` then the strip's
  transform and warp), the store and theme as any layer's. Internal this slice (above). Why a scope and not `section`:
  a section layer reads one placement; a bridge layer reads two.
- **An end's x**: its extent on each block that shows it (GA: every place of an overlapping block, each counted against
  `pairs.maxArcs`); a reversed block gives the end's first position the larger x. **Off the window** (a side's window
  shows none of the range): the end is placed toward its true x, clamped to one data-box width past the side's edge
  (JBrowse draws the true x; the clamp keeps a curve from a far chromosome from going flat), so the ribbon leaves the
  band sideways and moves continuously on pan until the clamp (EC-129's continuity).
  **In a hidden stretch** (a side with `pieces`, a collapsed intron, EC-31, EC-32): the part of an end in a hidden
  stretch is drawn at its break's hairline, at `ribbons.minPx`, and the card names it ("chr2:…, in a hidden stretch");
  an end wholly hidden keeps the ribbon (its foot on the break), never drops it. Why: a hidden stretch is a layout map,
  not data (EC-31); the ribbon's partner is still true. An end outside a block's `range` has no place (below).
  **No place** (its system shown, its position on no block: a sequence `include` left out, outside a block's
  `range`, GA): a **stub** (core's, in a bridge; a ring's outward stub stays RN's), a
  ribbon from the placed end tapering into the bridge for `ribbons.stubShare` (a quarter) of its height, hovered as the
  pair ("→ chr7:12,000,001–12,250,000, not shown"). **On an unread block inside a run** (GA's and RR's joint rule): at
  the block's place in its run, read from its other end's tile.
- **A ribbon's shape**: its two sides are cubics with control points at a third and two thirds of the band's height,
  vertical at both edges, so at height `u` (0 at side `a`'s edge, 1 at `b`'s) a side's x is
  `x_a + (x_b − x_a) · S(u)`, `S(u) = 3u² − 2u³`. A side joins **corresponding positions**: the a end's first position
  to the b position the pair maps it to (b's first when the strand product is 1, b's last when −1). So a ribbon is
  drawn crossed exactly when GA's rule says (`run direction × display direction of A × display direction of B = −1`),
  derived by geometry, not by a flag: reversing either block untwists it, as GA requires. A bin of strand 0 (RR: mixed)
  is never crossed.
- **Width at least `ribbons.minPx`** (1 device px) at each end, centred on the end; a ribbon whose ends both fall under
  that width draws as one curve (`shape: line`'s look). `shape` (the pairs layer's param, read in a bridge only:
  `ribbon`, absent; `line`, a curve between the ends' middles) forces lines. Why: a 10 kb run on a 3 Gb side is a
  hundredth of a pixel and must still be seen; a line is what an arc diagram draws.
- **An alignment's blocks merge under a pixel** (whole-genome scale): consecutive blocks of one alignment (rows sharing
  `row`, RR) of one strand whose gaps are under one device pixel on both sides are drawn as **one ribbon**, the hull of
  those blocks on each side; where a gap reaches a pixel, the ribbon parts there. Recomputed on a zoom, never on a pan
  (as GA's runs). Why: a chromosome-scale alignment is thousands of gapless blocks; at the fit they are one band, and
  thousands of slivers would spend the cull on one alignment.
- **Fill and order**: filled at `ribbons.fillOpacity` (0.35) in the layer's colour, no outline; drawn widest first
  (by the sum of its ends' widths), so narrow ones lie on top and stay hittable. **Colour**: by default `field block`,
  the block of the ribbon's end on the **bridge's** `a` side (whichever of the data's ends that is) over GA's
  block-ordered domain (EC-218, chromosome painting: each chromosome of the upper genome one hue, carried to the
  other); a side with no blocks has one block, so one colour. Twist already shows the strand.

### What is read, and at which level

- **Each side reads its own view** (RR's pairs, 1D tiles by either end, a pairs tile holding every pair with an end in
  it), each side at the level its own px pick (RR's rule: bins at the coarsest level no wider than `pairs.binPx`,
  items at a tile size; items count as the finest level). Why not one level for both: a level fine enough for the
  zoomed side is thousands of tiles on the other, and a coarse one is a single bin across the zoomed window.
- **Each pair drawn exactly once, by ownership**: the **near side** is the finer of the two (the fewer positions a
  device px over its window; items finer than any bins), a tie the bridge's `a`; the other is the **far side**. A pair
  or bin belongs to a side by its end's **anchor** there (the end's first position, the one RR bins it by). The near
  side **owns the region `R`**: its window on each of its systems widened outward to the far side's bin boundaries (the
  far level's grid on the near side's systems; with the far side on items, `R` is the window itself).
  - The near side reads `R`: its window at its own level, the widening (at most one far bin each side) at the coarsest
    level whose tiles fit `maxTilesPerWant`, never coarser than the far level; it draws every pair or bin anchored in
    `R`, wherever its far end is (off the far window: leaving sideways).
  - The far side reads its window at its own level and draws every pair or bin anchored in its window whose near
    anchor is **outside `R`**.
  - Exact: levels are a factor apart (RR's `summaries.factor`, a source's listed levels nested), so `R` is a union of
    whole far bins and every coarser bin a union of finer ones; no far bin straddles `R`'s edge, and a pair anchored
    in `R` is the near side's alone. A pair anchored in both windows is drawn by the near side; anchored in neither, by
    none (it is in neither view). A source whose levels are not nested (a format's own) widens `R` to the far bins
    that touch it, the same rule.
  - Why the finer side: it holds the detail the reader zoomed for; the far side's coarse bins fill only what the near
    side does not show.
- **The gate, per side**: each side asks items when the bins it draws (after ownership) sum to at most `pairs.maxArcs`,
  kept to `pairs.itemsHysteresis` × it (RR's rule, per side), so a side at the fit with thousands of bins never keeps
  the zoomed side off items. Both sides' wants together within `maxTilesPerLayer`. Items dedupe by RR's identity rule.
- **Bins in a bridge** are ribbons of their bins' widths, crossed by the bin's strand, **opacity by aligned length**:
  the bridge asks RR's `sum` of the source's aligned length (a run's length for a map's runs or an alignment's blocks,
  a `block_length` field where the source has one; else the count), its share of the bin's width on the narrower side, on a square-root scale over the bins in
  view (`ribbons.binOpacity`, 0.15 to 0.8). Why aligned length, not count: one chromosome-scale alignment is thousands
  of blocks, many short ones are few bases; the bases aligned are what synteny shows.
- **One cull over both sides** past `pairs.maxArcs`: every ribbon drawn (items, merged blocks, bins, stubs) ranked by
  aligned length in positions (an item its length or `matches`, a merged ribbon its blocks' sum, a bin its `sum`),
  ties by span, the first `maxArcs` drawn, the rest a density strip of ends along each edge of the band (the built
  strip, laid along the edge). Said once (`seq.render.pairs.limit`). Why one unit: bins and items then compare as the
  bases they carry.
- **Pending** (EC-169): a side's tiles not loaded draw the held coarser bins; with none held, the band's pending hatch
  over the stretch of its edge in view on that side, never empty. **Failed**: the placeholder with its issue.

### Hit testing between two sections (GA: per drawn run)

- **The bridge's `at`**: a page point in its band is a bridge point (`kind: 'data'`), else the sections' as built. A hit
  carries `bridge` (its index in `ResolvedView.bridges`, −1 elsewhere), one flat field; hits stay flat and allocate
  nothing warm.
- **A ribbon's hit, closed form**: at the point's `u`, `S(u)` is one value for every ribbon; each ribbon's two sides are
  `x_a + (x_b − x_a)·S` at that height, and the point is in it when its x lies between them (widened to `pairs.hitPx`,
  4, for a thin one). A linear scan over the drawn ribbons (bounded by `pairs.maxArcs`), last drawn first (the
  narrowest on top): about 4,096 multiply-adds, inside the hover's 2 ms. A stub is hit by its taper the same way. Why
  not an index: the ribbons' x at a height changes with `u`, and the scan is already bounded by the screen.
- **What a hit is: the drawn run, one rule** (GA: "a ribbon's hit is its run (its two blocks and ranges)"; GF: "click:
  the block"; RR: "what is drawn is what is selected"): a hit names the ribbon as drawn, its two ranges on its two
  blocks: an item (a pair, a map's run, one block of an alignment), a merged ribbon (the blocks of one alignment drawn
  as one, above: their hull on each side), a bin pair, a stub (its placed range and its unshown end). Its element is
  the rect form `(rect :a <range on a's system> :b <range on b's system>)` (built for grids, RR's for bins), selected
  as a Location of two parts; the card names what it belongs to (the alignment by its `row` or `id`, the map's run by
  its row, the bin's level). An item's own address (`(item :row k)` under its node) is the card's link, not the hit.
  Why one rule: a click must select what the reader sees under the pointer, whatever the zoom merged or binned.
- **Checked as the inverse of drawing**: property tests over random ribbons, both directions, reversed blocks, off-window
  ends and stubs: a point inside the drawn polygon (sampled finely) hits exactly the topmost ribbon covering it.

### Hover, selection and the keyboard across a bridge

- **Hover a ribbon**: the pair hovered, both ends lit on their sides (a pair closes to its ends, built), the ribbon
  drawn at full opacity with the hover's ink (`seq-hover`) as its outline. The card: "peach chr2:1,200,001–1,450,000 ↔
  grape chr5:880,001–1,130,000 · reversed · 250 kb · identity 0.97", a bin's as RR's ("37 pairs · 1 Mb bins · zoom in to
  see them"), a stub's naming its unshown end.
- **Hover a position on a side**: as built, projected along the view's links (a map between the two sides' systems,
  GA's one map between assemblies: so a liftover's map lights the image on the other side); and **the ribbons with an
  end at the hovered position** are lit with their far ends (as arcs light the pairs ending at a hovered base). Why both:
  a map projects positions; pairs that are not a map (BEDPE, Hi-C links) reach the other side only through their
  ribbons.
- **Select a ribbon**: click selects its drawn run (above: both ranges); Shift+click extends from the anchor end on its
  system; Cmd/Ctrl+click toggles. **Double-click** focuses each side on the hull of what the ribbon belongs to (GF: an
  alignment's hull on each side; a lone pair, its ends; a bin, its two bins), one jump, both viewports in one history
  step; again, back.
- **A selection on a side**: projected onto the other along the links (as built); the ribbons with an end inside it
  drawn in the selection's accent. **`bridge.selectAcross`** (the bridge's menu: across that bridge; the context menu
  over a side: across every bridge touching it, so a middle genome with a bridge above and below gains both; no key by
  default): adds the far ends of those ribbons to the selection (one selection change). Why a command and not
  automatic: a selection is what the user chose; partners are added when asked.
- **The keyboard cursor** crosses sections with ↑ and ↓ (built: tracks have one line); **across a bridge** it lands on
  the image of its position: through the topmost ribbon covering the cursor's position (its corresponding position on
  the far end), else through the view's link when a map joins the systems, else the far side's window centre; announced
  with the side and the way ("grape chr5:880,123, through a reversed ribbon"). ↓ crosses the bridge below, ↑ the one
  above, so a middle section is never ambiguous. Why: following a ribbon must have a keyboard way, and ↓ already moves
  to the next section.
- **ARIA**: a bridge is a group labelled by its `label` (else "a to b"), its description the count in view ("1,204
  ribbons between peach and grape in view"), set when the frame settles; ribbons are not focusable one by one (there are
  thousands); the cursor's announcement reads them.

### Moving two sections joined by a bridge (follow, lock, align)

- **No follow across a bridge**: GA refuses a follow between unequal block lists and through a map read by region
  (EC-211); taken. Each side is a leading section, moved alone by default (JBrowse's default). Two sides joined by a
  follow (GA lifts it for equal block lists: one genome against itself) or one scroll group already move together.
- **`lock: px`** (a bridge param, a patch): a pan of either side by dx CSS px pans the other by dx; a zoom by a factor
  about a page x zooms the other by that factor about the same page x. Why px and not positions: blocks of two
  assemblies, scales and gaps differ; the same px keeps the ribbons' slant as the reader moved them (JBrowse's link
  views, HiGlass's locks). A lock is a relation of its own beside follows and scroll groups
  ([viewer-models](viewer-models.md)): it never derives one viewport from the other, so it never reads
  the data.
- **One move a section a gesture**: a gesture's move is carried through follows (to a leader, `lead`, built), scroll
  groups and locks as one closure over the view's sections, each section moved **once**, by the first relation that
  reaches it in that order (follow, group, lock); a section reached again is skipped. So a lock between a leader and
  its follower moves nothing twice (the follow carries it; the lock is then redundant, `seq.view.bridge.lock`, info),
  and a cycle through follows, groups and locks ends where it began. Why one closure: three relations each moving a
  section would move it two or three times for one gesture.
- **A drag or a wheel over the bridge's band** moves both sides (the px rule, whatever `lock`, through the same
  closure), Cmd/Ctrl+wheel and pinch zoom both about the pointer's x. Why: the band belongs to both, and a drag there
  selects nothing (no axis). Bound by `bridge.pointer` (`ViewerConfig`: `both`, the default; `none`). Its lower edge's
  grip resizes it (Layout, above).
- **`bridge.align`** (a jump, one history step; Shift+A with the cursor or pointer on a side: every neighbour across a
  bridge touching that side is aligned to it, a middle genome's both in one step; the bridge's menu has each direction
  alone): the other side reveals the image of this side's window, the hull on the far system holding the most aligned
  length among the ribbons with their near end in the window (through the view's link when a map joins them and no
  ribbon is read). Nothing in view: no move, said in the readout ("nothing aligned in view"). Why a jump and not a
  follow: it reads the ribbons held now, once, so a map arriving by region never moves a viewport under the reader
  (EC-211).

### A liftover (EC-212; cut last, Gate 0, 1)

- **The figure**, from the nodes above with nothing of its own: a section of the old assembly (or one chromosome), a
  bridge whose `pairs` layer's data is the alignment's map (GF's `map.source` over the `seq` runtime's runs; a host's
  chain adapter alike), a section of the new assembly; the old assembly's annotations drawn on the new section through the same map (a layer's data along a
  path, as built; by region, GA, EC-211); hover and selection projected through it. Unaligned stretches are gaps in
  the ribbons; a part of an item that maps nowhere is not drawn on the new side (EC-21).
- **Runs read by either side**: the new section's layers and its projections need the runs whose `to` range falls in
  its window, where GF reads `map.source` by the `from` side. A map's runs are a pairs source's rows, tiled by either
  end (RR), so a run is asked by region of `from` or of `to` (`RegionRequest.side`, RR): forward projection reads by
  `from`, a layer drawn on the `to` section by `to`. The `seq` runtime reads its pairs whole and indexes both ends
  (GF), so this costs no second read.
- **Many to one, one to many** (runs overlap on `to` after a duplication; a deletion maps nowhere): a projection is
  the union of every image, as built (EC-20); the card of a projected hover lists the images ("2 places").
- **Deferred**: `mapSummary` (a mapped share and the unmapped ranges for a host's readout) and the builder's
  `liftover(old, map, new)`. Why: the liftover is cut last (Gate 0, 1), the figure needs neither, and "the share
  mapped" is ambiguous where positions have several images (a share of the location with at least one image, or with
  exactly one?), a question for when a host asks (Decided, 15).

### Contact maps and dot plots across sections (EC-209, EC-223; RR's tiles)

- **A grid over two sections of blocks**: a grid section with no `of` and no blocks that `follow`s a section of blocks
  takes its blocks for x (built), and with **`follow_y` naming another section of blocks** takes that section's blocks
  for y (new; today y is x's blocks). Each block pair is a strip pair: its cells clipped at both blocks' edges, the
  blocks' labels on both axes, GA's gaps and runs on each axis, a reversed block flipping its side (EC-145). Still
  refused (`seq.view.section.system`): a grid following blocks in a stack or on a ring. Amends
  [grids-and-rows](grids-and-rows.md) (A grid over a complex) and the resolved view's Grids.
- **Reads**: RR's 2D tiles by block pair, the upper triangle only when x and y are one assembly and the coverage says
  `symmetric`, the whole-assembly level at the fit; a block pair is read at the finer of its two axes' levels (RR).
- **A dot plot of two genomes**: the same grid section with a `pairs` layer (the alignment's sided pairs, the `seq`
  runtime): items as dots and
  diagonal segments with or against the diagonal by strand (built), bins as dots of their bin's size, opacity by count
  (as in a bridge). So an inversion is an anti-diagonal run, a translocation a run in another block pair's cell.
- **Between two chromosomes** (Juicebox's "chr1 vs chr4"): x follows one chromosome's section, y another's (built for
  systems: `follow_y` a section); a host's two lists patch their `of`.
- **The square navigator** over the grid is GA's (blocks × blocks, `navigator.square.maxCells`); nothing new here.

### Categorical grid axes

- **A grid axis may be a rows system** (a system whose positions are named members: an alignment's rows, RR's matrix
  states, the records of RR's `rows` shape): `of` or `of_y` names it, and the grid's data is over it (`x` or `y` a rows
  system; a sample × sample kinship, a gene × sample expression grid with positions on x). Why rows systems: they are
  systems already (EC-35), addressed, selected and ordered as rows, so a member selected on the grid lights its row in a
  matrix layer (the same system) with no map.
- **Its axis is one unit a member** in an order: the grid section's `member_order` (x) and `member_order_y` (y; named apart from the rows layer's `order` and an assembly's `order`, which the grid section also reads through its sides), the rows layer's
  values (`data`, `tree`, `{ field, direction }`, [grids-and-rows](grids-and-rows.md), Order), absent
  `data`. The axis transform holds the order's permutation (as the rows layer's: an `Int32Array` display offset to
  member and its inverse), so the grid's cells, its shared 2D reduction and its hits read members through it; a reorder
  rebuilds the grid's pyramid (its key holds the order's version). Why the rows layer's order: one grammar for ordering
  members, and a tree beside the samples' rows orders the kinship grid alike.
- **Chrome**: no ruler; the members' labels (the rows header's, `rows.labelMinPx` 8 px a label, thinned past it), down
  the label column for y, along the top at 45° for x (the plot's box grows by the labels' room, fixed: the box never
  changes with a zoom). Hover names both members and the value ("NA12878 × NA12891 · 0.42 kinship"); a click on a label
  selects the member (a `row` target, built), on a cell its rect (built), Alt+drag subtracts (built).
- **Mixed and moved**: a categorical axis follows nothing of positions (`follow` or `follow_y` from a positions section
  onto it is `seq.view.follow.unplaced`, built); zoom and pan act on members as on rows (`rowScale`-like, the grid's own
  `y`), `0` fits. A categorical axis larger than `maxBlockValues` cells over the grid is the grid's bound (built).

### Ribbons in other places (for RN and the rest)

- **Where ribbons and stubs live, one rule**: core draws ribbons and stubs **in a bridge only** (`pairs` there); in a
  section core's `pairs` draws arcs and chords (a bin pair one arc stroked by count, RR); ribbons and outward stubs on a
  ring are RN's (`rondo:ribbons`, fallback `pairs`). An extension type in a bridge falls back to `pairs` ribbons this
  slice (`scope: 'bridge'` internal, Decided, 16).
- **Ribbons inside one row** (an arc with width between two blocks of one section): not core's; an extension may draw
  them (`scope: 'section'`). Why: no Genomes view asks for it (the synteny views are across sections, the ring is RN's).

## Edge cases this design takes or assigns

- **Taken here**: EC-209 (both halves: ribbons between sections, and a grid whose x and y follow two sections of
  blocks), EC-206 (the twist by corresponding positions, GA's rule), EC-212 (the liftover figure; the
  runs' reading GF's `map.source`, its projection GA's), EC-223 (contact maps across two sections of blocks; tiles RR's), EC-129
  (an off-window end moves continuously: ribbons leave sideways at their true x), EC-169 (pending in a bridge: held
  coarser bins, else the hatch along the edge), EC-235 and EC-238 to EC-244 (a reader without bridges drops them, the
  sides whole), EC-218 (a ribbon's `block`, the a end's), EC-35 (rows systems as categorical grid axes), EC-145 (the
  crosshair on categorical axes, built).
- **Taken from the other designs, not redesigned**: EC-207 and EC-210 (GA: one map between assemblies; two assemblies on one
  section), EC-211 (GA: no follow through a map read by region; here `align` and `lock` instead), EC-208 and EC-62
  (RR: pairs read by region, found from either end; here per side), EC-171 (RR: no "load anyway").
- **Assigned**: EC-213 (MAF, cut); EC-225 (a minimum length or MAPQ: the query sublanguage's numeric
  filter).

## Defaults, with a picture

What a reader sees with no setting but the bridge's `a`, `b` and one `pairs` layer over an alignment's sided pairs (the builder's
`synteny(a, b, data)` writes the three nodes; nothing inferred).

1. **Two assemblies, the fit view**: peach above, grape below, a 96 px band between; a ribbon a run, filled at 0.35 in
   its peach chromosome's hue; an inversion crossed; a translocation a ribbon to another chromosome; a gap no ribbon.
   ```
   peach  |1        |  |2          |  |3      |
          ░░░░  ▒▒▒▒   ▓▓▓  ▓▓▓▓▓▓     ████
           ╲ ╲   ╲╲     ╲╳╱  │││        ╲
            ╲ ╲   ╲╲    ╱╳╲  │││  (gap)  ╲___
   grape  |1        |  |2          |  |3     |  |4 …   (chr3's tail went to 4)
   ```
2. **Zoomed out, bins**: ribbons a bin wide on each side, opacity by count; zoomed in past the gate, items replace them
   with no gap.
3. **One side zoomed in**: the zoomed side's runs fan to their true x on the other; ends past a side's edge leave the
   band sideways.
   ```
   peach   |   chr2:1.2–1.5 Mb   |
            ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
             ╲               ╱
              ╲─────────────╱
   grape  |1 | |2 ▏| |3  | |4 |   (the ribbon narrows to grape chr2's 250 kb)
   ```
4. **A stub**: an end whose chromosome the side's `include` left out: a quarter-height taper from the placed end.
5. **The liftover**: the old assembly above, the map's runs as ribbons (straight, a few narrow gaps, a duplicated run
   fanning to two places), the new assembly below with the old one's genes drawn through the map; a hover on a gene
   lights its image below.
6. **A dot plot of two genomes**: x peach's blocks, y grape's, each block pair a cell; runs as diagonal segments, an
   inversion against the diagonal.
   ```
           1      2       3
   1   ╲      ▏       ▏
        ╲     ▏       ▏
   2    ▏ ╲   ▏   ╱   ▏      (an inversion in 2 × 2)
        ▏     ▏ ╱     ▏
   ```
7. **A contact map between two chromosomes or two sections**: RR's picture 5, its axes following two sections.
8. **Categorical axes**: members' labels down the left and at 45° along the top, cells at least `grid.rectMinPx` as
   rects with outlines, the crosshair lighting a member's row and column.
   ```
              NA1 NA2 NA3
            ╱   ╱   ╱
   NA12878 ▓▓  ▒▒  ░░
   NA12891 ▒▒  ▓▓  ░░
   ```
9. **Pending, failed, nothing**: the pending hatch along the edge in view of the side not loaded; the placeholder with
   its issue; nothing in view: an empty band, its ARIA description "no ribbons in view".
10. **A reader without bridges**: the two sections stacked with the sections' gap, the band gone.
11. **Dark**: tokens only (the categorical palette's dark variant at the same opacity, `seq-hover`, the selection's
    accent, the pending hatch); nothing hard-coded.

## The interaction matrix's cells

Cells this design fills; the others stand as GA, RR, grids-and-rows and ux-foundations give them.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element: a ribbon** (a pair, a run, merged blocks, a bin, a stub) | click: the drawn run, `(rect :a … :b …)`, both ends lit | Shift+click: from the anchor end | Cmd/Ctrl+click | double-click: each side to the hull of what it belongs to (an alignment, a pair, a bin), one jump | double-click again, Alt+Left | — | — | no: a data filter (length, MAPQ, `type`) | no: data |
| **a bridge** (its band) | no: a bridge is not selected | — | — | Shift+A at a side: `bridge.align` (a jump); the menu has both ways | Alt+Left | Cmd/Ctrl+wheel, pinch over it: both sides about the pointer's x | drag or wheel over it: both sides (`bridge.pointer`) | `layer.hidden`; `hidden` (patch) | `lock`, `height`, `label` by patch; its menu toggles `lock` (one undo step) |
| **a section beside a bridge** | as built | as built | as built | as built | as built | alone; with `lock: px` both | alone; with `lock: px` both | as built | as built |
| **the selection** across a bridge | — | — | `bridge.selectAcross`: the far ends of its ribbons (menu, no key) | Enter: as built | Alt+Left | — | — | — | as built |
| **the keyboard cursor** at a side's last track | — | Shift+↓: no (a range is on one system) | — | ↓ or ↑ across a bridge: the image of its position, announced | Alt+Left | — | — | — | — |
| **a grid over two sections of blocks** | as built (a cell, a rect) | as built | as built | double-click: a zoom step (built); on a block's label: that axis to the block (GA) | Alt+Left | as built (both axes) | as built | as built | `follow`, `follow_y` by patch |
| **an element: a categorical member** (a label) | click: the member (`row`) | Shift+click: the run in display order | Cmd/Ctrl+click | double-click: the member's row and column shown | Alt+Left | — | — | no (a filter of the rows) | no |
| **a categorical axis** | — | — | — | — | Alt+Left | as the grid's: members | as the grid's: members | — | `member_order`, `member_order_y` by patch |
| **a bookmark or link** to a ribbon | — | — | the drawn run's rect as a ref (positions, never the level or the merge) | runs it: each side to its range | Alt+Left | — | — | — | its label (planned) |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| `bridge.a`, `bridge.b` | spec params | — (required) | the two sections joined; their neighbours in the stack |
| `bridge.label` | spec param | "a to b" (the sides' labels) | names the band in the label column and to ARIA; any text |
| `bridge.height` | spec param | `bridges.height` | the figure's; any px |
| `bridges.height` | `RenderConfig` preset | 96 CSS px | room for a crossed ribbon to read (JBrowse's level: about 100); any |
| `bridge.lock` | spec param | `none` | each genome moved alone, as JBrowse; `px` |
| `pairs.shape` (in a bridge) | spec param | `ribbon` | a pair's extents are what synteny shows; `line` |
| `mirror` (pairs layer, built) | spec param | false | a pair drawn where its data puts it; true adds the swap on one genome against itself |
| `ribbons.fillOpacity` | `RenderConfig` preset | 0.35 | overlaps read through; any |
| `ribbons.binOpacity` | `RenderConfig` preset | 0.15 to 0.8, square root of the count in view | strong bins read strong without hiding weak ones; any range |
| `ribbons.minPx` | `RenderConfig` preset | 1 device px | a short run on a long axis is still seen; any |
| `ribbons.stubShare` | `RenderConfig` preset | 0.25 of the band | a stub says "elsewhere" without reaching the far side; any |
| `pairs.maxArcs`, `pairs.hitPx`, `pairs.binPx`, `pairs.itemsHysteresis` | `RenderConfig` preset (built, RR's) | 4,096; 4; 2; 1.25 | now read by bridges too; any |
| `bridge.pointer` | `ViewerConfig` | `both` | the band belongs to both sides; `none` |
| key of `bridge.align` | `ViewerConfig` (`view.keys`) | Shift+A | a keyboard way to bring the other genome over; any (checked against the bindings at LA) |
| `follow_y` a second section of blocks | spec param (grid section) | `x` (built) | a section of blocks: the grid's y over its blocks |
| `member_order`, `member_order_y` (grid section, categorical axes) | spec params | `data` | the members' own order; `tree`, `{ field, direction }` (the rows layer's values) |
| `of`, `of_y` a rows system (grid section) | spec params | — | a categorical axis; a positions system as built |
| `bridge.hidden` | spec param | false | kept in the spec, its room collapsed; true |
| keys of `bridge.selectAcross`, the lock toggle | `ViewerConfig` (`view.keys`) | none (menus) | rare acts, and keys are scarce; any |
| `bridge.resize` grip | `RenderConfig` preset (`bridges.gripPx`) | 4 CSS px | JBrowse's resizable levels; 0 for none |
| categorical x labels at 45° | fixed | — | upright labels overlap at a member a few px wide; vertical ones read sideways |
| `aspect` with `follow_y` a second section | spec param (built) | 1 | the box fixed from the width, as built; the cells rectangles when the two windows differ in span (built) |
| the colour of a ribbon | the layer's channel | `field block` (the a end's block) | chromosome painting (EC-218); `strand`, a field, `b_system` |
| ends placed by system, `a` to `a` when both sides show it; no path for pairs | fixed | — | the reasons in the Decision |
| a bridge's sides neighbours and rows | fixed | — | a ribbon must not cross a side's tracks |
| the twist by corresponding positions; ends at their true x, clamped one box width out | fixed | — | GA's rule by geometry; continuity on pan |
| each side read at its own level, a pair drawn once by ownership (the finer side's region widened to the far bins) | fixed | — | one ribbon drawn once at the level its side needs |

Not built: every row (LA; `follow_y`'s second section and categorical axes LA too).

## Tensions

- **No follow across a bridge** (GA, EC-211): two genomes drift apart as the reader moves one. Options: a follow
  through a map read whole at open (a chain), `lock: px`, `bridge.align`. Recommendation: lock and align now; a follow
  only if the use gate asks (Decided, 7).
- **Sides are rows only**: a ring's links stay inside it, a text section takes no bridge. Standing: two rings joined by
  ribbons has no precedent and no view asks for it.
- **Neighbours only**: A against C across B needs a reorder. Standing.
- **Long translocations are nearly horizontal** in a 96 px band and overdraw others. Accepted: `height` is the figure's,
  colour by block tells them apart, the cull keeps the widest.
- **Two levels in one band**: the near side's items or fine bins beside the far side's coarse bins read as two
  densities, and the near side's region is widened by up to a far bin each side (read coarser there).
  Accepted: each is true at its scale; the card says which is a bin.
- **A bridge reads two want sets**, up to twice a layer's tiles. Bounded by `maxTilesPerLayer` (RR); measured in LA on
  the whole-genome fixture.
- **Pairs that are not a map do not project**: a BEDPE's partners light through ribbons and `selectAcross` only.
  Standing: a projection through pairs would be a composed claim.
- **A categorical reorder rebuilds the grid's pyramid**. Accepted: a reorder is a patch, not a frame.
- **Tracks beside a grid** (Juicebox's, HiGlass's left tracks) wait for sections side by side (layouts). Open.

## Self-check

The claim: **this design covers what a user does with layers across sections.**

1. **Assume it holds.** A bridge joins two stacked sections; its `pairs` layer draws each pair as a ribbon between
   them, twisted where GA's rule says, hit by closed form; each side moves alone, locked or aligned; hover and selection
   cross through ribbons and maps; a grid's two axes follow two sections of blocks; a grid axis can be members.
2. **What it rests on, negated.** (a) Two sections; plural: three genomes, A against C, one genome against itself. (b)
   The sides are rows; negated: a ring, a text, a stack, a grid. (c) Both ends are in view; negated: one off the window,
   on a chromosome not shown, on an unread block, on two overlapping blocks, nowhere. (d) Both sides at one scale;
   negated: one zoomed in a thousandfold. (e) The data is a map; negated: pairs only (BEDPE, Hi-C links), bins. (f) The
   user moves one side; inverse: the other must come along; symmetric: the user acts in the band itself. (g) The user
   acts with a pointer; negated: the keyboard, a screen reader. (h) Axes are positions; negated: members. (i) The reader
   has bridges; negated: one without, an older minor. (j) Data arrives; negated: pending, failed, past the screen. (k)
   What is drawn is true; negated: a ribbon drawn at the wrong end, twice, or uncrossed.
3. **The negated world.** (a) chains of bridges, a refused non-neighbour, `mirror` for self-synteny; (b) a ribbon off a
   ring's far rim, a text line in the middle of its block; (c) sideways ends, stubs, the joint run rule, a ribbon a
   place; (d) a level per side and a rule against drawing twice; (e) ribbons lit through their ends, `selectAcross`,
   opacity by count; (f) `lock`, `align`, the band dragged; (g) the cursor across the band, ARIA counts; (h) the
   categorical axis, its order, its labels; (i) the bridge dropped, the sides whole; (j) held coarser bins, the hatch
   along an edge, the cull; (k) placement by system, corresponding positions, the hit as the drawing's inverse.
4. **Both held.** Integrated: every item of (a), (c) to (k) in the Decision. Standing: (b) rings and text as sides
   refused (a ring's links stay inside it); a follow across a bridge (f, Decided 7); A against C (a); tracks beside a grid
   (h's neighbour, layouts' side by side).

## Consequences

- **`seqx-core`** (LA, stage 3, after RG and FR; its language part declared by GC1's bump, as GF's: Decided, 12):
  - language: the `bridge` kind (`a`, `b`, `label`, `height`, `lock`, `hidden`; children `layer`); `pairs.shape`; a
    `pairs` layer's `data` a segments `map`; a grid section's `follow_y` a second section of blocks, `of`/`of_y` a rows
    system, `member_order`, `member_order_y`; the builder's `synteny(a, b, data)`; issues `seq.view.bridge.lock`,
    `seq.lang.bridge.sides`, `seq.view.bridge.neighbours`, `seq.view.bridge.layout`, `seq.view.bridge.layer`,
    `seq.view.bridge`, `seq.render.bridge.unplaced`.
  - resolved view: `ResolvedView.bridges` (`ResolvedBridge`: address, `a`, `b`, label, height, lock, layers, issues) and
    `ResolvedView.stack` (sections and bridges in order); a grid's `y` over a second section's blocks; `GridAxis`
    gains `categorical` and its order.
  - render: `layout/bridge.ts` (`BridgePlacement`, `BridgeSide`, its `at`, its key); `layoutView` placing bridges;
    `LayerType.scope: 'bridge'` and `BridgeContext`; `pairs-ribbons.ts` (ribbons, stubs, bins, the per-side reads and
    the drawn-once rule, the closed-form hit, emphasis); `Hit.bridge`; the categorical axis (permutation transform,
    labels, the pyramid's key); export (a bridge a group, `data-bridge`, `data-a`, `data-b`; ribbons as cubic paths);
    nothing in the extension entry (`scope: 'bridge'`, `BridgeContext` internal). Lazy as the Layout says.
- **`seqx-viewer`**: the bridge's pointer (drag and wheel both), `lock` carried in the `ViewModel` beside follows and
  groups, `bridge.align` (`navigate`), `bridge.selectAcross` (`select`), the lock toggle (`patch`), the cursor across a
  bridge, the cards, ARIA.
- **`seqx-examples`** (GD's): peach and grape with a bridge (inversion, translocation, gap); one genome against itself
  with `mirror`; three genomes chained; the liftover (generated, its truths computed in code); the dot plot
  of two genomes; a Hi-C between two chromosome sections; a kinship grid on categorical axes beside the genotype matrix.
- **Storybook**: a Design page "Across sections"; a story per default above, light and dark, canvas and SVG, each with
  its headless twin on a fake clock (tiles arriving per side).
- **Tests**: property tests (the closed-form hit equals the drawn polygon's topmost cover; the twist by geometry equals
  GA's three-direction rule over every combination of block and section directions; each pair drawn exactly once
  (items and bins, every pair of levels on two sides, windows overlapping or not, a tie in scale)
  across two sides at two levels; a bridge's reads equal the whole read clipped to either window; `lock` moves round
  trip; `align` lands on the hull it names); drawing against the data (every ribbon's ends at its data's positions on
  each side, at the fit, across a tile edge, across a block gap, inside a run, reversed blocks); EC-129, EC-169,
  EC-206, EC-209, EC-212, EC-218, EC-223, EC-235 named by tests. Performance, each with its path: a bridge frame of
  4,096 ribbons under 8 ms (quick representative at 1,024), its hover under 2 ms at 0 bytes warm (quick), a pan of one
  side relaying that side and the bridge only (quick, counted), the whole-genome synteny at the fit within
  `maxTilesPerLayer` (heavy, generated); sizes: the chunk about 6 kB, measured at LA's gate.
- **`_spec` sections that change**: core Language (bridge, pairs, grid section), Resolved view (Bridges, Grids:
  categorical axes, a second section's blocks), Render (Section layouts: the stack and bridges; Layer types: scopes;
  The pairs layer: ribbons; The grid layer: categorical axes; Lanes and Hits: `bridge`; Export; Configuration); viewer
  Models (the lock), Input, Commands, Bindings, Configuration.
- **Other designs amended**: [grids-and-rows](grids-and-rows.md) (Tensions "Arcs between sections"
  closed; the self-check's standing "ribbons and categorical grid axes" closed; Pairs: in a bridge, data a map; A grid
  over a complex: `follow_y` a second section); [layouts](layouts.md) (the view stacks sections and
  bridges; `scope: 'bridge'` in the layer contract); [maps](maps.md) (a segments map as pairs data);
  [viewer-models](viewer-models.md) (the lock); [ux-foundations](ux-foundations.md)
  (the matrix's ribbon and bridge rows); [extension-openings](extension-openings.md) (the entry exports
  the bridge context; `scope: 'bridge'`); [features](features.md) (151, 152 planned with LA);
  [edge-cases](edge-cases.md) (EC-209's outcome).

## Decided by the author

2026-11-21: every gap as recommended (the author on the drafts: "generally reasonable and we can fold in the reviews";
on the whole-designs review's twelve recommendations: "yes with the recommendations"; on the fold-in: "yes"). An item
the author decided by name carries its own date.

1. **A bridge as a new kind under `view`** (layers only, no tracks), not a section with `layout: bridge`: **decided yes**.
2. **Its name `bridge`** (JBrowse says "level", Circos "link"; `link` is the view's resolved links already): **decided yes**.
3. **Sides must be the bridge's neighbours in the stack**; A against C needs a reorder: **decided yes**.
4. **Sides in the `tracks` layout only, blocks in a row**; a ring's links stay inside it (chords, RN's ribbons), text
   and stacks refused: **decided yes**.
5. **The `pairs` layer draws ribbons in a bridge** (one type, `shape: ribbon | line`), no core `ribbons` type: **decided yes**.
6. **A segments map may be a `pairs` layer's data** (an alignment's runs: one node for ribbons and projection): **decided yes**.
7. **No follow across a bridge**; `lock: px` (a spec param), a drag or wheel over the band moving both, and
   `bridge.align` (a jump) instead: **decided yes**. (A follow through a map read whole at open: no, until the use gate asks.)
8. **Each side reads its own window at its own level**, a pair drawn once by ownership (RR's level rule): **decided yes**.
9. **↓ and ↑ across a bridge land on the image of the cursor's position**: **decided yes**.
10. **`bridge.selectAcross`** (the far ends added on asking, never automatically): **decided yes**.
11. **Ribbons coloured by the a end's block by default** (chromosome painting), the twist showing the strand: **decided yes**.
12. **LX's language additions ride GC1's bump** (`bridge`, `pairs.shape`, a map as data, `follow_y` a second section,
    a rows system on a grid axis, `member_order`, `member_order_y`), so the language changes once: **decided yes**.
13. **Tracks beside a grid** (Juicebox's left tracks): not this slice, with sections side by side: **decided yes**, defer.
14. **Categorical grid axes as rows systems** with the rows layer's order values (`member_order`, `member_order_y`) and labels: **decided yes**.
15. **`mapSummary` and the `liftover()` builder deferred** (the liftover is cut last; "the share mapped" is ambiguous
    with several images): **decided yes**, defer.
16. **`scope: 'bridge'` and `BridgeContext` internal**, not in the extension entry, until an extension draws in a
    bridge: **decided yes**.
17. **Ownership by the finer side, `a` on a tie** (each pair drawn once; the far side's bins cut at the near side's
    region widened to the far grid): **decided yes**.
18. **A hit is the drawn run** (a block, blocks merged under a pixel, a bin), selected as its rect; the alignment named
    in the card and focused by double-click: **decided yes**.
19. **One closure moves each section once a gesture** (follows, groups, locks): **decided yes**.
