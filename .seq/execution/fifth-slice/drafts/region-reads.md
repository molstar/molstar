# Pairs and grids by region: the answer shapes, 2D tiles, the pairs level

Draft of lane RR (the fifth slice, stage 0) for Gate 0's agreement; on agreement it moves to `design/` and amends
[protocol](../../../design/protocol.md), [grids-and-rows](../../../design/grids-and-rows.md) (Data, the pairs and grid
data rules, Tensions), [data-store](../../../design/data-store.md) (summaries), [layouts](../../../design/layouts.md)
(an image on a ring), [performance](../../../design/performance.md) (structures, scale) and
[files](../../../design/files.md) (held pairs from files) (below, Consequences). Built on Gate 0's answers as
recommended ([plan](../plan.md)): bundling is **core's pairs level**, which the ring extension reads (4); MAF is cut
first (1), so wide alignments are here only as what falls out; public API changes stay in `packages/seq` (no
`@molstar/core` or `@molstar/io` change). Built on GA's draft ([genome-coordinates](genome-coordinates.md)): assemblies,
blocks, the places index (`placesOf`), runs of sub-pixel blocks, the twist rule, the square navigator's coordinates.

## Intent

The protocol deferred four answer shapes to Genomes ([protocol](../../../design/protocol.md): "`pairs`, `matrix`,
`grid`, `rows`: with the Genomes milestone"); grids-and-rows deferred "a request's second axis (`RegionRequest.y`), 2D
tiles, the store's grid summaries, pairs and grids read by region (EC-208, EC-223)", and the ring refuses images
("dense matrices on rings wait for the Genomes milestone", [layouts](../../../design/layouts.md), Rings). The milestone
needs them for: synteny pairs from a PAF or BEDPE over a whole genome (EC-208), a Hi-C map by region and level
(EC-223), samples × positions (genotypes, methylation, expression: EC-220 to EC-222), results keyed by id (EC-224's
data, not its join), and the ring's 12,000 links bundled by count (research-genome-ring, Q2: "a `pairs` level read by
region"). Cost follows pixels ([performance](../../../design/performance.md)): each read below is bounded by the screen,
never by the data.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| HiGlass and Clodius: 2D tiles of 256 × 256 bins (`uuid.z.x.y`), zoom levels a factor 2 apart, a tileset-info answer (bounds, resolutions), the upper triangle served and mirrored | 2D tiles of 256 bins a side, square levels, a store tile one image tile, symmetric grids served once | a preprocessing server before anything shows; tiles by zoom index (here bin widths in positions, the protocol's levels) |
| cooler (`.cool`, `.mcool`): a bin table and a sparse pixel table (`bin1_id`, `bin2_id`, `count`, upper triangle), resolutions coarsened by summing, balancing weights per bin | the sparse form (`xBin`, `yBin`, `value`) with its `fill`; coarse levels by `sum` for counts; a balanced value as a channel (`field`) | bins of uneven width (restriction fragments): bins here are a level's width |
| Juicebox and `.hic` (straw): a matrix per chromosome pair, resolutions in bp (2.5 Mb … 5 kb), normalizations (`NONE`, `KR`, `VC`, `SCALE`) and observed/expected chosen per read | a grid over blocks read by block pair; the normalization and `oe` as `RegionRequest.field`; the source's own levels as listed | a resolution picker as the reader's control (here the level fits the pixels; a host may still pin one through the source) |
| pairix and 4DN `.pairs`, BEDPE with tabix: pairs indexed by both ends (pairix) or by the first only (tabix) | "a pairs tile holds every pair with either end in it": the runtime's contract; one indexed by one end reads whole at open | a 2D index for pairs in the protocol (pairs are tiled in 1D here) |
| Rondo: connections clustered hierarchically, precomputed at several resolutions, a cluster's own links shown only zoomed in | the pairs level: bin-pair counts at each level, items only where they fit | the clustering itself (RN's picture, over the bins) |
| JBrowse: "zoom in to see features" by a byte estimate from the index; EC-171's zoom gate | items gated by the counts of the level above, measured, not estimated | "force load" (refused: below) |
| Circos heatmap tracks: a sector a bin a row on a ring | a matrix on a ring as sector runs, one run a device column a row at most | a polar bitmap (a second drawing path; a Tension) |
| This codebase: 1D tiles of 1,024 bins, value summaries (factor 16, read whole up to 2^24 positions, else binned per tile), `cells` (rows tiled, columns whole), the shared 2D reduction and image LRU, `readTiles`, `boundedLevel`, `RegionRequest.field` | every one, extended: `y` on a request and a want, the 2D key, summaries for pairs and grids, `matrix` as `cells`' mirror (positions tiled, states whole) | a second store, a second cache, a level rule a layer invents |

### Interactions in the precedents

| What a user does there | Here |
|---|---|
| pan and zoom a Hi-C map, tiles arriving coarse to fine (HiGlass, Juicebox) | built: the grid section's gestures (built) over 2D tiles; a pending tile draws the coarser held one |
| choose a resolution (Juicebox) | refused as a gesture: the level fits the pixels; a host or writer lists the levels (the source) |
| choose a normalization, observed/expected (Juicebox, HiGlass) | built: `grid.field` (by patch) |
| hover a bin: its two ranges and its value (HiGlass, Juicebox) | built: the hover card (below, A hover over a bin) |
| a contact map over the whole genome, chromosome pairs (Juicebox "All") | built: a grid section over assemblies (GA's blocks), read by block pair |
| zoom into a cluster of links to see its links (Rondo) | built: a zoom reads items where the counts fit |
| filter links by a minimum count (Rondo) | built: the layer's data filter on `count` (EC-225's numeric filter, the query sublanguage) |
| export the links in view as BEDPE, the matrix as a table (Rondo) | refused here: data export is a host's; the SVG and PNG are built |
| a heatmap ring of samples (Circos) | built: a `matrix` layer on a ring |
| load anyway past the zoom gate (EC-171) | refused: bins are drawn instead, and say how many; reading items past `pairs.maxArcs` draws no more than the screen holds |

## Decision

### Whose each choice is

By [axis](../../../design/axis.md)'s rule. No new port: the shapes are the protocol's, the levels the store's and the
layers'; formats are GF's runtimes. What changes the meaning of a figure (which source, which channel) is a spec param;
how a read is planned and bounded is a store setting or a `RenderConfig` preset.

### The answer shapes

`Shape` gains `pairs`, `grid`, `matrix`, `rows`. Every table keeps a `form` (protocol, Tables), so one shape has an
encoding a level:

```ts
// pairs: 1D tiles over an end's system; a tile holds every pair with either end in it
interface PairsList {                 // the finest level the source lists
  form: 'list';
  aSystem?: CategoryColumn; bSystem?: CategoryColumn;  // absent: both ends on the requested system
  aStart: Int32Array; aEnd: Int32Array; bStart: Int32Array; bEnd: Int32Array;
  aStrand?: Int8Array; bStrand?: Int8Array;            // −1, 0, 1 (EC-206)
  type?: CategoryColumn; id?: CategoryColumn; row?: Int32Array; fields?: RangesField[];
}
interface PairsBins {                 // a coarser level: bin-pair counts (the pairs level)
  form: 'bins';
  aSystem?: CategoryColumn; bSystem?: CategoryColumn;
  aBin: Int32Array; bBin: Int32Array;                  // bins of the level's width on each end's system
  strand?: Int8Array;                                  // a × b strand: 1 straight, −1 reversed (the twist), 0 none
  category?: CategoryColumn;                           // the `type` binned by, when the coverage lists it
  count: Float32Array;                                 // pairs in the bin pair
  sum?: Float32Array;                                  // the requested `field` summed, when asked
}
// grid: 2D tiles; one value a cell at every level
interface GridDense  { form: 'dense'; x: number; y: number; nx: number; ny: number; value: Float32Array } // x-major, NaN missing
interface GridSparse { form: 'sparse'; x: number; y: number; xBin: Int32Array; yBin: Int32Array; value: Float32Array;
                       fill: number }                  // cells not listed: 0 (counts) or NaN (missing)
// matrix: 1D tiles over positions, every state (sample) whole: cells' mirror
interface MatrixItems { form: 'items'; start: Int32Array; end: Int32Array;   // a column an item: a variant, a CpG, a gene
                        states: number; value?: Float32Array; codes?: Uint8Array; // item-major, states × items values
                        id?: CategoryColumn; row?: Int32Array }
interface MatrixBins  { form: 'bins'; start: number; states: number;        // a bin per column, a level's width
                        mean?: Float32Array; min?: Float32Array; max?: Float32Array; count?: Float32Array; // values
                        classes?: Uint8Array; shares?: Uint8Array }          // codes: majority and share (as cells)
// rows: keyed records, on the source's records system (a position a record)
interface RowsKeyed  { form: 'keyed'; start: number; key: CategoryColumn; fields: RangesField[] }
```

- **Pairs are tiled in 1D, by either end.** A request `{ system, start, end, level }` answers every pair with an end
  overlapping the region (EC-62: found from either end), the other end anywhere; a pair with both ends in view is in two
  tiles and readers dedupe it by **its row's identity and both ends** (`row`, else `id`, with both ends' systems and
  starts: BEDPE writes `.` for every unnamed pair, so an `id` alone would merge distinct rows; an `id` of `.` or empty is
  never an identity); bins by their whole key (both systems, both bins, strand, category).
- **Items are read at a tile size, not a bin.** The `list` form may be asked at **any listed level**: the level is then
  only the tile's width (`level × tileBins` positions), and the answer is every item with an end in that tile, uncut.
  The request says so (`RegionRequest.items: true`, the key's level part `i<level>`, so an items tile and a bins tile at
  one level never share a key). The reader asks items at the **finest level whose tiles over the view fit
  `maxTilesPerWant`** (`boundedLevel` with no pixel floor): a 10 Mb view reads items in 10 tiles of 1,024 kb, not
  9,800 tiles of 1,024 positions; a genome ring at the fit reads them a chromosome a tile or two. A source with one
  level serves an items tile at any level as a region read of that tile (the summaries' wrapper maps it). Why: with the
  level as a bin only, the gate opened over a wide view asked past the tile bound and stayed pending for good. Why 1D: a ring, an arc track and a ribbon read pairs by the region of
  one end, and a whole genome's pairs over 25 chromosomes are 25 tiles, not 625; a dot plot of pairs in a grid section
  filters the x tiles by its y window (pairs are sparse: a tile's rows are few where a grid's cells are not). Why
  both ends: a pair is seen from either end (an arc from chr4 must show the link to chr1). **A runtime indexed by one
  end** (tabix over BEDPE) serves the other by reading whole at open, within the open's bounds (`limits.maxRows`, `limits.maxBytes`, built), and indexing both (GF's
  runtimes; else `too-large` with the hint).
- **The pairs level is bin-pair counts.** At a level coarser than the finest a source lists, `pairs` answers `bins`:
  each bin pair with at least one pair, its count, by relative strand (so a ribbon's twist, GA's rule, is a bin's own:
  straight and reversed runs of one bin pair never merge) and by `type` when the coverage lists it (`ShapeInfo.binBy`),
  with a field summed when `field` is asked (a PAF's matches, a link's score). Why counts and not a sample of items: a
  count is exact at every level and sums across levels (the pyramid's check, below); a sample would show what it
  dropped as absent. Why in the protocol and not a layer: one level serves core's arcs, RN's ribbons, the square
  navigator and a dot plot of pairs (Gate 0, 4).
- **Grids are tiled in 2D**: `RegionRequest.y?: { system, start, end }`, the grid's second axis, at the **same level**
  (bin width) as x: square levels, as HiGlass, cooler and `.hic` have them. Why square: cells stay square at aspect 1
  and one bin width names a tile. A tile is `gridTileBins` (256) bins a side, clipped to each system's length (a
  chromosome shorter than a tile is one small tile, never a padded 256²). Why 256, not the 1D tiles' 1,024: 1,024² float
  cells are 4 MB a tile (the eviction's grain, and 16 tiles fill the cache); 256² is 256 kB, and one store tile is one
  image tile (the image LRU's 256 × 256 samples), so a region grid needs no second reduction.
- **Symmetric grids are served once**: the coverage marks a grid on one system, or on a pair of systems of one
  assembly, `symmetric` (Hi-C); the store asks only the tile at or above the diagonal (the x tile not after the y tile,
  by the coverage's system order across two systems) and the reader transposes the other. Why: Hi-C formats store the
  upper triangle; asking both doubles the bytes for the same cells. EC-59's two grids on one pair (upper and lower of
  different data) stay two grids, each with its own `triangle`.
- **One value a grid cell, aggregated as the source says**: the coverage's `aggregate` for a grid (`sum` for counts,
  `mean` for scores; `mean` when absent) is how its coarse levels were made; `field` picks the channel (`count`,
  `balanced`, `oe`: the normalizations as fields, `RegionRequest.field`, built). Why one value: four floats a cell
  (min, max, mean, count) make a dense tile 1 MB and the view's tiles past the cache; the grid layer's `aggregate`
  (`mean`, `max`, `min`) still reduces samples finer than a device pixel (held grids as built; region grids only where a
  source's levels jump). Standing: a region grid has no `max` across levels (a Tension).
- **A matrix is cells' mirror**: positions tiled, states whole. A tile at the finest level is the items over the region
  (a variant, a CpG, a gene: each a column with its own `start` and `end`, so a long deletion is one item across the
  SNPs it covers, EC-222) × every state, `value` (float32, NaN missing: methylation, expression, EC-220, EC-221) or
  `codes` (uint8: genotypes `0`, `1`, `2`, 255 missing); at a coarser level a bin per column, `mean`/`min`/`max`/`count`
  for values, majority `classes` and `shares` for codes (cells' coarse form). The states are a system (as an
  alignment's rows, EC-35), its names from the coverage (`ShapeInfo.states`), so a sample is addressed, selected and
  ordered as a row is. Bounded **in bytes**: a tile within `maxTileBytes` (16 MB: about 4,000 states of float values
  at 1,024 bins, 16,000 of codes); more is `too-large` (a Tension: a states window by `y` would lift it). The long form (a row a sample × position) is a
  records transform (files.md's), converted once at open (EC-220's "wide and long forms").
- **`rows` are keyed records** on the source's records system (a position a record, in the source's order): `key` and
  fields, read by region of that system, held whole within `summaries.maxBytes` when something joins by key. Why a system: the
  store reads and evicts it as any 1D shape, and the key index is built once from the held tiles. The join itself
  (EC-224: results keyed by id colored onto transcripts) is not here: a theme or a query reading a key (Gaps, 8).
- **Wide alignments fall out, built only if MAF stays** (Gate 0, 1: cut first): `cells` may tile its columns instead
  of its rows (`ShapeInfo.tile: 'columns'`), rows whole, exactly the matrix's tiling; a MAF block over a genome is then
  a few dozen species × 1,024 columns a tile. Nothing else is designed for it.

### Bounds are bytes, for every shape

A count of cells is not a size: 2^24 float cells are 64 MB, the whole cache. Every bound below is in bytes, each a
`StoreSettings` item, and each answer's bytes are counted before it is kept (`checkAnswer`, as built for its columns):

| Shape | One answer (a tile) | A whole read (summaries, a one-end index) |
|---|---|---|
| `pairs` | `maxTileBytes` (16 MB: about 600,000 items at about 28 bytes, or 1M bins at 16) | `summaries.maxBytes` (16 MB: about 600,000 items) |
| `grid` | `maxTileBytes`; a dense 256² tile is 256 kB, a sparse cell 12 bytes (so sparse past a third of a tile is refused as dense's job) | `summaries.maxBytes` (about 4M float cells: a 2,000 × 2,000 grid) |
| `matrix` | `maxTileBytes` (about 4,000 float states or 16,000 code states at 1,024 bins) | — (never read whole) |
| `rows` | `maxTileBytes` | `summaries.maxBytes` when joined by key |
| `cells` (built) | `maxBlockValues` stays its bound, now also `maxTileBytes` | as built |

`maxTileBytes` (16 MB) and `summaries.maxBytes` (16 MB) are a quarter of the cache each, so one answer or one pyramid
never empties it. `maxRows` and `maxBlockValues` stay as the parse's guards before allocation (security.md).

### Coverage

`ShapeInfo` gains, all optional: `y` (the grid's second systems: a list, or `'*'`, every system of the coverage: a
genome's 25² pairs are not listed one by one); `symmetric`; `aggregate` (`sum`, `mean`); `binBy` (a category the pairs
bins keep: `type`); `states` (a matrix's state names, in order); `tile` (`cells` only: `'rows'` or `'columns'`);
`fields` (as built: a grid's or matrix's channels, a pairs field that can be summed); `whole` (a grid's or pairs'
**whole-assembly system**: `{ system, sequences }`, a system whose positions are the listed sequences end to end in
that order, served at the coarse levels only: Juicebox's "All", below). Why on the shape: one source (a
`.hic`) serves one grid over every pair; a VCF serves ranges and a matrix over the same systems. `Coverage.levels` stays
the source's; a format's own (`.hic`'s nine resolutions) are listed as they are.

### The store: keys, tiles, wants

- **The 2D key** is the 1D key with the second axis: `<source>|<x system>|<level>|<x tile>|<y system>|<y tile>`, each
  part escaped as built. Six parts never equal four, so 1D and 2D entries never collide. `TileRegion` and `Want` gain
  `y?` (system, start, end; for the region, the y tile's positions). Why one key family: eviction, settle, priority,
  grace and retry are the store's and stay one code path.
- **A 2D want's tiles** are its x tiles × its y tiles, symmetric ones folded; `maxTilesPerWant` (256) bounds the
  product, so a grid want past it is `seq.store.too_many_tiles` (as built) and the layer goes coarser first
  (`boundedLevel` with the product). **Margin** on a 2D want is one ring of tiles round the view
  (`grid.marginTiles`, 1), at `margin`, evicted first; no prefetch.
- **A bound on a layer's wants together** (`maxTilesPerLayer`, 1,024): a layer over blocks makes a want a block or a
  block pair, each within `maxTilesPerWant`, and their sum is bounded too; past it the layer reads coarser, and past
  the coarsest it reads the whole-assembly level (below), else `seq.store.too_many_tiles`. Why: a whole-genome Hi-C at
  the fit is about 325 block-pair wants (4,096 under the square), thousands of range requests on one `.hic` file.
- **The whole-assembly level** (Juicebox's "All"): a source over an assembly may list a `whole` system (Coverage); a
  grid or pairs layer over blocks reads it, one want for the section, whenever the level that fits is at or coarser
  than the whole system's finest level, or the block wants pass `maxTilesPerLayer`. The layer cuts each block pair's
  rect from the whole tiles by the blocks' offsets in the whole system (so GA's `order`, `include` and directions
  need not match the source's: a reversed block flips its cut). The store's summaries build it for a one-level source
  read whole within `summaries.maxBytes`. Why: at the fit a genome is one image of a few million cells, not 325 files'
  worth of requests.
- **Bytes**: an entry weighs its tables (a dense tile `nx × ny × 4`, a sparse one 12 bytes a cell, pairs and matrices
  their columns); the cache's bound (64 MB), its least-recently-used order and never evicting the view are as built.
  The reference view (a 1440 × 900 plot at DPR 2, every cell a device pixel) is at most 13 × 9 dense tiles, about
  30 MB, and about 42 MB with its margin ring: inside the bound. **A second grid on the same pair** (EC-59: upper and
  lower of different data) doubles it, past the bound: `seq.store.over_budget`, and the layer with the most bytes reads
  a level coarser (a quarter of the bytes), as built. EC-59's "size cap per block pair" is `maxTileBytes` a tile and
  `summaries.maxBytes` a whole read, per block pair.
- **Pending in 2D** (EC-169): a tile not loaded draws the held coarser tile under it; with none held, a **pending
  rectangle** (the pending hatch over its rect, clipped to the plot), hit as pending ("loading"), never as empty or
  missing values.
- **Pairs and matrices use the 1D key** on the system read; **rows** on their records system.

### Summaries: the store's pyramids for one-level sources

Beside the value summaries (built, `summaries.ts`), for a source listing one level only, owned by the attached source,
counted against `maxBytes`, gone with it:

- **Built in time slices** on the main thread (performance.md: "on the main thread in time slices until the compute
  port exists"), `summaries.sliceMs` (4 ms) a frame, the source's coarse levels pending meanwhile (bins drawn as they
  land); a 600,000-pair pyramid or a 4M-cell grid is a few dozen slices, never one long task.
- **Pairs**: a source whose items read whole fit `summaries.maxBytes` (16 MB, about 600,000 pairs) is read whole once
  at its level (every system's region), each coarser level's bins (`summaries.factor`, 16, as values) built from the
  items in one pass a level, then the items let go (the source serves them again by region, as items tiles). **Only
  the levels that halve the bins are kept**: a level whose bin pairs are more than half the level below's (or, for the
  first, half the items) is dropped from the coverage the wrapper gives, since fine levels of sparse pairs barely
  compress and would cost bytes and a read for nothing; the reader then takes the next kept level. More pairs: a coarse
  tile is binned on the fly from an items read of that tile alone within `maxTileBytes`, else `too-large` (a format with
  levels, said once: `seq.store.pairs.levels`). Why the values' factor: the pairs tiles are 1D, their density per device
  column the same question.
- **Grids**: a source whose finest cells over every listed pair fit `summaries.maxBytes` (about 4M float cells: a
  2,000 × 2,000 grid) is read whole once and summarized at factor 2 a side (square levels, a third more bytes), `sum` or
  `mean` as its `aggregate`, with its whole-assembly level when it is over an assembly. More: a coarse tile is reduced
  on the fly from the finer tiles under it while they fit `maxTileBytes` (two levels up), else `too-large`
  (`seq.store.grid.levels`: a Hi-C genome at 5 kb is a `.hic` or `.mcool` with its own levels, never summarized in a
  page). Why factor 2 for grids: a factor-16 step in cells is a factor 4 a side, and the view would hold up to 16 times
  its pixels in cells between levels.
- **Held pairs get the same levels** (inline pairs, a tRNA's, a complex's links, a small generated genome's):
  `binPairs(data, level)`, the function the store's summaries call, computed once per level from the held
  `PairsData` and kept by its version. One reader interface over both backings (`PairsReader`: items in a region, bins
  at a level), so the pairs layer and RN's ribbons read held and region pairs alike. Why: the ring's 12,000 links are
  as likely inline (generated) as a file, and the bundled picture must not depend on where they came from.

### How a level is chosen at a zoom

One rule a shape, each in `render/layer-data.ts` beside `layerLevel`, each bounded by `maxTilesPerWant` as built:

- **Values, ranges, letters, cells** (as built): the coarsest level no wider than `k` device pixels' positions.
- **A matrix**: the same along positions (`k` = `values.coarser`, 2); states whole.
- **A grid**: the **finest level whose bin is at least `grid.binPx` (1) device pixel** on both axes (the larger of the
  two axes' positions a pixel), then coarser while the tiles exceed the bound. Why the other way from values: a cell
  under a pixel is reduced again for nothing, and the tiles in view are what bounds the bytes (above); this is the
  rule the held grid's image pyramid already follows ("the level whose sample is at least a device pixel").
- **Pairs**, two reads: **bins** at the coarsest level no wider than `pairs.binPx` (2) device pixels' positions, always
  (bounded by the screen: at most a bin pair a pair of device columns, in practice far fewer); then **items**, read as
  items tiles at the finest level whose tiles fit the bound (above), **with hysteresis**: they are asked when the bin
  pairs with an end in view (deduplicated across tiles and blocks, their counts summed: pairs, not whole tiles) sum to
  at most `pairs.maxArcs` (4,096), and kept while the loaded items with an end in view stay at most
  `pairs.itemsHysteresis` × `pairs.maxArcs` (1.25: 5,120); past it, back to bins. Bins draw until the items arrive;
  items replace them with no gap (EC-169: pending is never empty). Why gated by counts: the count is known exactly from
  the bins already read, where EC-171's gate guessed from tiles; and a zoom into a cluster shows its links as soon as
  they fit the screen, Rondo's behaviour. Why hysteresis: a pan at the threshold would otherwise flip bins and arcs
  every frame. **EC-171's outcome**: its catalog's "load anyway forces" is refused; bins show what is there, and items
  past the screen's bound would be culled again.
- **Past `pairs.maxArcs`, the built cull orders by strength**: the pairs layer's device-column binning (built: past
  `maxArcs`, one arc a bin of end columns, the rest a density strip) stays for items between `maxArcs` and the
  hysteresis bound and for held pairs; for bins it keeps the strongest bin pairs (by count, ties by span) up to
  `maxArcs`, the rest in the density strip of ends. One cull, ordered by count for bins, widest first for items.
- **Copies** (EC-62's "contacts fold by copy group"): bins on blocks that are copies fold by copy group as items do
  (built): one arc a group, hit as every member; a bin between two systems on split copies whose data names no unit is
  drawn on no copy (`seq.render.pairs.copy`, as built).
- **On blocks** (GA): each block reads at its own scale's density (a scaled block finer); a grid over blocks reads a
  block pair at the finer of its two axes. **Small blocks**: the joint rule with GA, below. Why: 3,000 contigs ×
  layers are more wants than pixels.
- **On a ring**: the device columns are the outer circumference's (EC-108), as built.

### Small blocks: drawn in runs, read while a device pixel wide (the joint rule, for GA's draft word for word)

> **Small blocks run together, and are read while they have a pixel.** Consecutive blocks each narrower than
> `blocks.runPx` (2 CSS px, `RenderConfig`) at the current scale are drawn as one run: one band over their joined
> extent, no gaps or labels inside, a hairline at each end, labelled when it has room. Inside a run a block's layers
> read and draw its data at its place while the block is at least one device pixel wide (`blocks.readPx`, 1 device
> px); a block under one device pixel reads nothing at that zoom, and its stretch of the run shows the run's band
> alone, its hover naming the block and "not read at this zoom". Pairs with an end on an unread block are still read
> from their other end's tile and drawn at the block's place (`placesOf`); a pair with both ends on unread blocks is
> not drawn, and the run's hover counts it. Runs and the read set are recomputed on a zoom, never on a pan. Why two
> thresholds: a gap and a label need two CSS pixels to read, data needs one device pixel to be seen; under a device
> pixel a block's data is more requests than pixels (3,000 contigs × layers), and its pairs, the data a genome view
> reads across blocks, survive through their other end.

At DPR 2 the run starts at 4 device px and reading stops at 1, so a run's blocks between are drawn with their data;
at DPR 1 at 2 px and 1 px.

### Host adapters (amends seq-data-ports.md)

A host's adapter is a runtime ([seq-data-ports](../../../design/seq-data-ports.md)) and may serve any of the new
shapes: `y` on a grid request, pairs `bins` and items tiles at any level, `matrix`, `rows`. The conformance suite gains
the new shapes and applies to every adapter as to a runtime: each listed shape read by region equals the whole read
clipped, at every level it lists; pairs found from either end; a symmetric grid's transposed tile equal to the tile
asked; answers within `maxTileBytes`. An adapter that serves a shape at one level only gets the store's summaries
like any source. Why: a host holding a contact map in memory (a structure viewer's predicted contacts) is code Seq\*
does not own, read through the same bounds.

### What the pairs layer and the grid layer draw from them

- **The core `pairs` layer at a bins level** draws each bin pair as one arc (a chord on a ring, a dot of the bin's size
  in a grid) between the two bins' middles, its stroke by count (`pairs.binStroke`, 1 to 6 CSS px on a square-root
  scale over the counts in view), twisted per GA's rule where `strand` is −1, coloured by its `category` or the layer's
  theme (`count` is a field the theme reads). Each end is placed at every block showing its bin (GA's places index;
  a bin cut by a block's edge is clipped to each block). `pairs.maxArcs` still bounds what is drawn. Why one arc and
  not a ribbon: the ribbon, its trumpets and its bundling are RN's picture over the same bins; core draws the honest
  minimum a viewer without the extension shows.
- **The grid layer over a region grid** draws each store tile as one image tile (its cells coloured by the value
  theme, `fill` honoured, NaN transparent), cached in the shared image LRU by (layer, the entry's key and version, the
  colours' version, the triangle); a pending tile draws the held coarser tile under it, never empty; a failed one its
  placeholder; with no coarser tile held, a pending rectangle (above). Rects at `grid.rectMinPx` and above as built. A
  region grid on one system draws its mirror from the served triangle.
- **A `matrix` layer** (new type, `axes: 2` in its band, fallback `values` of the states' mean): the states as rows in
  the band, with the rows layer's band machinery (its scroll `bandTop`, its zoom `rowScale`, labels from the states
  system, `order` by a field of a states table when held, `hidden`); items at their own extents, a long one solid
  across the columns it covers; coarse bins as image tiles of states × device columns (the shared reduction for rows
  finer than a device pixel, majority for codes). Why its own type and not `rows`: an alignment's rows draw letters on
  one columns system; a matrix's columns are items or bins of a genome. Its emphasis along y is the rows layer's.

### Dense matrices on a ring (`seq.render.ring.image`, lifted for rows)

- **An image on a polar strip is drawn as runs**: each of its sample rows becomes runs of equal colour along the angle,
  each run a rect the polar warp turns into an annular sector (built), at most one run a device column (at the outer
  radius) a row; SVG export draws an image as runs of rects already, so a ring gets one path for both. So a matrix's
  states are thin concentric rings (Circos's heatmap tracks), an alignment's rows on a ring the same. Bounded:
  `rings.imageRuns` (65,536 a frame); past it rows are taken from the next coarser row level of the shared reduction,
  said once (`seq.render.ring.image`, now info: "rows reduced on the ring"). A `grid` layer is never on a ring (a ring
  has one axis; it is laid out in a grid section only), so the issue stays for it. Why runs and not a polar bitmap: the
  bitmap is a second drawing path, canvas only, its SVG an embedded PNG; runs reuse the warp and the export (a
  Tension: measured in RG).

### A hover over a bin, and what a bin is when acted on

- **A pairs bin** hovers as its two bin ranges, its count and its level, never as a pair: "chr1:12,000,001–13,000,000
  ↔ chr4:50,000,001–51,000,000 · 37 pairs (contact 30, other 7) · 1 Mb bins · zoom in to see them" (the by-category
  breakdown when `category` is held, at most `hover.maxCategories`, the rest "+3 more"); a strand −1 bin says
  "reversed". Its element is `(rect :a <bin on a> :b <bin on b>)` (the built rect form): a click selects the rect,
  projected as its two ranges on every section (both ends light, as a pair's); a double-click focuses both ranges'
  hull (as a pair's), which reads finer bins and, past the gate, the items.
- **A grid cell read at a coarse level** hovers as its bin on each axis and its value with its channel and aggregate:
  "chr1:12.000–12.005 Mb × chr4:50.000–50.005 Mb · 1,204 (count; balanced 0.83) · 5 kb bins" (a `sum` level says
  "sum of 4 bins" only when the source's level is finer than the cell drawn); a sparse tile's unlisted cell reads its
  `fill` ("0"), a NaN "no data", a pending tile "loading". A click selects the bin's rect (what is drawn is what is
  selected); a double-click zooms a step about it (built). `Hit` gains one flat field, `bin` (the level read, in
  positions; absent at 1), so the card computes the ranges and a hit still allocates nothing.
- **A matrix cell** hovers as its item, its state and its value or genotype ("rs123 (chr7:117,559,590) · sample
  NA12878 · 0/1"), at a coarse level the bin and its mean or majority with its share.

## Edge cases this design takes or assigns

- **Taken here**: EC-208 (pairs read by region, either end, the pairs level), EC-223 (a binned contact map by region and
  level), EC-220 to EC-222 (the matrix shape and layer: wide and long forms, missing values, a long deletion solid),
  EC-62 (found from either end, by region too; bins fold by copy group as items), EC-206 (a bin's strand), EC-169
  (pending bins and tiles drawn coarse, a pending rectangle in 2D with none held, never empty), EC-170 (the level rules
  by shape), EC-171 (outcome: its "load anyway forces" refused; the gate is by counts, bins drawn instead), EC-172 and
  EC-173 (2D wants ordered and cancelled as 1D), EC-125 (an export waits for 2D tiles as for any), EC-108 (a ring's
  columns and an image on it), EC-59 (two region grids on one pair merge later over earlier, as held ones; its size cap
  per block pair is `maxTileBytes` a tile and `summaries.maxBytes` a whole read).
- **Its data part here, the rest assigned**: EC-224 (the `rows` shape here; the join a theme's or a query's: Gaps, 8),
  EC-209 (a contact map's tiles here; the layer across sections LX's), EC-213 (the column tiling falls out; MAF cut
  first, Gate 0, 1), EC-225 (a minimum count as the data filter's numeric form, the query sublanguage's).
- **Assigned**: the bundled picture of the pairs level (RN's `rondo:ribbons`), each format's reads and bounds (GF's
  runtimes: `.hic` and `.mcool` are not in the plan's formats list; a Gap, 10).

## Defaults, with a picture

What a reader sees with no setting but the data's `source` and a layer.

1. **Pairs zoomed in** (a few Mb of a synteny PAF, under 4,096 pairs in view): items as arcs, as built.
   ```
     ╭───────────╮      ╭──╮
   ──┴──────┬────┴──────┴──┴───  chr4
            ╰─ (to chr1)  caret
   ```
2. **Pairs zoomed out** (the whole chromosome): bin pairs, each an arc stroked by count; the hover names the bin pair.
   ```
     ╭━━━━━━━━━━━━━━━╮   ╭─────╮
   ──┸───────────────┸───┸─────┸──  chr4   (thick: 120 pairs; thin: 3)
   ```
3. **A genome ring at the fit** (core's pairs layer, no extension): bin-pair chords stroked by count between the
   chromosomes' arcs, the strongest drawn first and kept by the cull; with RN's extension, ribbons over the same bins.
   ```
              chr1
         chrX ╭────╮ chr2
            ╭─╯ ┃╲ ╰─╮
       chr22│   ┃ ╲  │chr3     (thick chord: chr1 ↔ chr17, 140 pairs;
            ╰─╮ ┃  ╲╭╯          thin: 3 pairs)
              ╰─┸───╯ chr4 …
   ```
4. **A Hi-C map by region**: tiles of 256 bins at a device pixel or more a bin; a pending tile drawn from the held
   coarser one (blocky, then sharp), the plot's box fixed (built).
   ```
   ┌────────┬────────┐
   │▓▓▒▒░░  │░░      │
   │▒▓▓▒▒░  │ ░░     │   (coarse tile, sharpening)
   │░▒▓▓▒▒░ │  ░     │
   ├────────┼────────┤
   ```
5. **A Hi-C map over the genome** (a grid section over an assembly, GA's blocks): at the fit, one read of the
   whole-assembly level cut into block pairs; zoomed in, block-pair tiles; the upper triangle served, the lower
   mirrored, chromosome gaps as GA's, labels on both axes.
   ```
          chr1      chr2     chr3 …
   chr1 ▓▓▓▒▒░░▏ ░░ ░  ▏ ░  ░ ▏
        ▒▓▓▓▒░░▏ ░ ░░  ▏  ░   ▏
   chr2 ░ ░░   ▏▓▓▓▒░  ▏ ░░   ▏
        ░░  ░  ▏▒▓▓▓▒  ▏  ░ ░ ▏
   chr3 ░  ░   ▏ ░ ░   ▏▓▓▓▒  ▏
   ```
6. **Samples × positions** (a VCF's genotypes, 40 samples): a row a sample, an item a variant at its extent, a long
   deletion solid; zoomed out, a bin's majority genotype, its share as the colour's strength.
   ```
   NA12878 ▏▌ ▌▌   ███████   ▌ ▏
   NA12891 ▏ ▌ ▌             ▌▌▏
   ```
7. **A matrix on a ring**: each state a thin ring inside the chromosomes' arcs, at most a run a device column a ring.
   ```
           ╭────────╮   chromosomes (bands)
         ╭─┼▒▒░░▓▓▒─┼─╮ sample 1
        ╭┼─┼░▓▓▒░░▒─┼─┼╮ sample 2
        ││ │  hole  │ ││
   ```
8. **Bins with nothing in view**: no arc, no placeholder; pending: the band's pending hatch (built); failed: the
   placeholder with its issue (built).
9. **Dark**: tokens only (the value and categorical themes, `seq-hover`, the pending hatch); nothing hard-coded.

## The interaction matrix's cells

Cells this design fills or changes; the others stand as grids-and-rows and GA give them.

| | select | extend | add / remove | focus | back | zoom | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element: a bin pair** | click: its rect (both bins, both ends light) | Shift+click: the hull from the anchor | Cmd/Ctrl+click | double-click: both bins' hull (finer bins, then items) | Alt+Left | — | — | no: a data filter on `count` or `type` | no: data |
| **an element: a grid cell at a coarse level** | click: the bin's rect | Shift+click: the rect from the anchor | Cmd/Ctrl+click, Alt+click subtracts | double-click: a zoom step about it (built) | Alt+Left | — | — | no | no |
| **an element: a matrix cell** | click: the item on its system and the state | Shift+click: the block from the anchor | Cmd/Ctrl+click | double-click: the item | Alt+Left | — | — | no | no |
| **an element: a state (a matrix row's header)** | click: the state | Shift+click: the run in display order | Cmd/Ctrl+click | double-click: the state shown | Alt+Left | — | — | `rows.hide` (patch, as rows) | no |
| **the selection over a matrix** (items × states) | — | Shift+arrows: the block of items and states | Cmd/Ctrl+Shift+arrows start a block | Enter: the items' hull, the states shown | Alt+Left | — | — | — | the selection's bookmark command (planned, grids-and-rows) |
| **a bookmark or link** to a bin or matrix block | — | — | a bin's rect or a matrix block as a ref (positions, never the level) | runs it: both ranges (a bin), the items and states (a matrix) | Alt+Left | — | — | — | its label (planned) |
| **a navigator over a matrix's states** (the vertical overview, navigators.md) | — | — | — | a click centres the states there | Alt+Left | the window's grips: states | drag its body | the host hides it | — |
| **a track or layer: a region grid** | — | — | — | — | — | the grid section's (built) | the grid section's (built) | `layer.hidden` | `field` by patch (a normalization) |
| **a track or layer: pairs by region** | — | — | — | — | — | — | — | `layer.hidden` | `field`, the data filter by patch |
| **a leading section: a matrix's band** | — | — | — | — | — | rows: as the rows layer (`]`, `[`, Cmd/Ctrl+wheel over the header) | the wheel's vertical delta over the band | — | `order` by patch |
| **the view: the level read** | no: the level fits the pixels | — | — | — | — | — | — | — | no: a gesture forcing items past `pairs.maxArcs` is refused (the screen holds no more) |

## Configuration

| Item | Whose | Default | Why the default; the alternatives |
|---|---|---|---|
| `pairs.source`, `grid.source`, `matrix.source`, `field` on each | spec params | — (instead of inline `rows`, `values`) | what the figure reads (EC-208, EC-223); `field` a normalization or a summed column |
| `pairs.a`, `pairs.b`, `grid.x`, `grid.y`, `matrix.target` an assembly | spec params | — | a file over every chromosome (GA's `target` an assembly) |
| `gridTileBins` | `StoreSettings` | 256 | 256 kB a dense tile, one image tile; any power of two |
| `maxTileBytes` | `StoreSettings` | 16 MB | one answer never more than a quarter of the cache; any |
| `summaries.maxBytes` | `StoreSettings` (`summaries`) | 16 MB | a whole read for a pyramid (about 600,000 pairs, 4M grid cells), a quarter of the cache; more: per tile |
| `summaries.sliceMs` | `StoreSettings` (`summaries`) | 4 ms a frame | a pyramid built in slices on the main thread; any |
| the one-end index's whole read (GF's runtimes) | `OpenLimits` (`maxRows`, `maxBytes`, built) | 1,000,000 rows; 64 MB | the open's bounds as built; past them `too-large` |
| `maxTilesPerLayer` | `StoreSettings` | 1,024 | a layer's wants over blocks together; past it coarser, then the whole-assembly level |
| `grid.marginTiles` | `RenderConfig` preset | 1 (a ring of tiles) | a pan's next tiles, evicted first; 0 |
| `pairs.itemsHysteresis` | `RenderConfig` preset | 1.25 | items kept to 1.25 × `maxArcs`: no flip at the threshold; ≥ 1 |
| `blocks.readPx` (joint with GA) | `RenderConfig` preset | 1 device px | a block's data is read while it can be seen; 0 reads every block |
| `navigator.square.maxCells` (GA's) | `RenderConfig` preset | 4,096 | bounds the square's block-pair wants and its lattice; any |
| the grid summaries' factor | fixed | 2 a side | square levels; the view's bytes stay a few times its pixels (above) |
| the pairs summaries' factor | `StoreSettings` (`summaries.factor`, built) | 16 | the values' step, one setting; any |
| `grid.binPx` | `RenderConfig` preset | 1 device px | a cell under a pixel is reduced again for nothing; 2 halves the tiles |
| `pairs.binPx` | `RenderConfig` preset | 2 device px | a bin pair a couple of pixels wide; any |
| `pairs.maxArcs` | `RenderConfig` preset (built) | 4,096 | now also the items' gate; any |
| `pairs.binStroke` | `RenderConfig` preset | 1 to 6 CSS px, square root of the count | thick reads as many without hiding the rest; any range |
| `rings.imageRuns` | `RenderConfig` preset | 65,536 a frame | a ring's sector runs bounded; past it the next row level |
| `hover.maxCategories` | `RenderConfig` preset | 4 | a card stays short; any |
| the matrix's states bound | `StoreSettings` (`maxTileBytes`) | about 4,000 float states (16,000 codes) at 1,024 bins | a tile within its bytes; a states window later (Tensions) |
| pairs tiled in 1D by either end; grids in 2D at square levels; one value a grid cell; symmetric grids served once | fixed | — | the reasons in the Decision |
| the level rules by shape | fixed | — | one rule a shape, each with its reason; their bounds are the presets above |

Not built: every row above but `pairs.maxArcs`, `summaries.factor` and the `OpenLimits` (RG; `blocks.readPx` GC2,
`navigator.square.maxCells` GC3).

## Tensions

- **No `max` across a region grid's levels**: one value a cell keeps the bytes within the cache; a PAE-like `min` or a
  Hi-C `max` zoomed out reads the source's aggregate. Recommendation: accepted until a region grid asks (a source may
  serve a second field, `max`, as a channel).
- **Items behind bins cost a round trip at zoom-in**: the counts must arrive before items are asked. Accepted: bins
  draw meanwhile, and the margin's bins are read ahead.
- **A matrix's states are whole**: past about 4,000 float states (16,000 codes) a tile is `too-large`. Open: a states window (`y` on a matrix
  request, the 2D key) when a dataset asks (a biobank).
- **One-end indexes**: a BEDPE under tabix is read whole at open for its other end; a large one is `too-large` and needs
  a pairix-like or a levelled format. Standing; GF's.
- **Images on rings as runs** cost more than a polar bitmap at full density. Recommendation: runs (one path, the
  export's), measured in RG against `rings.imageRuns`; the bitmap if the measure says so.
- **Images in SVG** (grids-and-rows' Tension) grow with region grids: an export of a Hi-C view is runs of rects a tile.
  Standing until the PNG-in-SVG decision.
- **The pairs level is not a clustering**: bins are a fixed lattice, so a cluster across a bin edge splits. Accepted: a
  lattice is exact, cheap and sums across levels; clustering over bins is RN's picture.
- **Reads inside runs** (GA's Tension): closed by the joint rule (read while a device pixel wide); its cost is that a
  run's values under a pixel are not drawn at the fit. Open until GC2's measure on the 3,000-contig fixture.
- **The whole-assembly level is a source's or the summaries'**: a `.hic` without it, over many contigs, still makes a
  request a block pair up to `maxTilesPerLayer`, then refuses. Standing until a real file is measured (GF, M5).
- **Summaries on the main thread**: time-sliced, so a 600,000-pair pyramid lands over a few dozen frames. Standing until the
  compute port.

## Self-check

The claim: **this design covers what a user does with pairs and grids read by region.**

1. **Assume it holds.** Pairs, grids, matrices and keyed rows come from a source by region and level through the
   store; each layer picks its level by one rule and draws within the screen; a reader hovers, selects and zooms into
   bins as into items.
2. **What it rests on, negated.** (a) The data fits in memory; negated: 100M Hi-C contacts, a genome at 5 kb. (b) A
   source has levels; negated: one level only (a PAF, a `.cool` of one resolution). (c) Pairs are on one system;
   plural: across chromosomes, two assemblies. (d) The reader zooms in; inverse: the whole genome, a ring. (e) Cells
   are square, positions × positions; symmetric case: samples × positions. (f) A grid is symmetric; negated: two
   genomes, upper and lower of different data (EC-59). (g) A value is a count; negated: a score, missing against zero.
   (h) Data arrives; negated: pending, failed, too large, an index on one end only. (i) One consumer; plural: the view,
   the square navigator and an export on one store. (j) A bin is a fact; negated: a bin is an aggregate a reader may
   mistake for one pair. (k) A track is straight; negated: a ring.
3. **The negated world.** (a) the store's summaries refusing whole reads past a bound and a format with levels named;
   (b) summaries for pairs and grids; (c) bins keyed by both systems, ends placed through GA's index, the twist kept by
   a bin's strand; (d) the pairs level, items gated by counts, blocks under a pixel unread; (e) the matrix, states
   whole, cells' mirror; (f) `symmetric` declared, never inferred, two grids by `triangle`; (g) `aggregate` and `fill`;
   (h) coarse tiles under pending ones, `too-large` with a hint, a one-end index read whole; (i) one key family and
   want sets shared (built), the square reading by block; (j) a hover that says the bin, its count and its level, a
   selection of the bin's rect; (k) images on rings as runs.
4. **Both held.** Integrated: every item above in the Decision. Standing: no `max` across region grid levels (g);
   states whole (e); one-end indexes read whole (h); images on rings as runs, unmeasured (k); a lattice, not a
   clustering (j); reads inside runs, until measured (d).

## Consequences

- **`seqx-core`** (RG, stage 3, alone in the store):
  - data: `protocol.ts` (the shapes and forms above, `RegionRequest.y`, `ShapeInfo`'s new members), `keys.ts` (the 2D
    key, `tileSpan` for 2D, symmetric folding), `store.ts` (`Want.y`, `keysOf` and `tooMany` over the product, the
    margin ring, bytes of the new forms), `summaries.ts` (pairs bins and grid levels, `binPairs`), `tables.ts`
    (`checkTables`, `checkAnswer`, `sliceTable`, `joinTables` for each new form; bins and cells within the region and
    the level), `seq-format.ts` (the new shapes round trip), the conformance suite (a region read equal to the whole
    read clipped, at every level, for each shape; either end found).
  - language (a minor version): `pairs.source`, `pairs.field`, `grid.source`, `grid.field`, `rows` and `values` then
    optional; the `matrix` kind (`target`, `source`, inline `values` or `codes` with `states`, `field`); the `matrix`
    layer type; `target`s an assembly (GA). Issues `seq.store.pairs.levels`, `seq.store.grid.levels`,
    `seq.render.matrix.states`; `seq.render.ring.image` moved to info for rows.
  - resolve: access plans for `pairs`, `grid`, `matrix` by source (held ones unchanged); the states system from the
    coverage (as GA's derived sequences enter resolve, `ResolveOptions`).
  - render: `layer-data.ts` (`gridLevel`, `pairsLevels`, `readTiles` with `y`, the bins-then-items read), `pairs.ts`
    (`PairsReader` over held and region pairs, bins drawn, the gate), `grid.ts` (store tiles as image tiles, `fill`,
    mirrored tiles, `Hit.bin`), the polar warp (images as runs, `rings.imageRuns`), the `matrix` layer (on the rows
    layer's band machinery), the hover texts.
- **`seqx-viewer`**: the hover card's bin and matrix texts; nothing else (gestures as built).
- **`seqx-examples`** (GD's): a generated Hi-C over the toy assemblies (symmetric, two levels and a one-level twin for
  the summaries), a synteny PAF by region, 12,000 generated links (held and as a source: one picture), a genotype
  matrix with a long deletion, methylation with missing values.
- **Storybook**: a Design page "Region reads"; stories per default above, light and dark, canvas and SVG, each with
  its headless twin on a fake clock (tiles arriving coarse to fine).
- **Tests**: property tests (a region read equals the whole read clipped at every level, every shape; the pairs
  pyramid's counts sum to the items' at every level and per bin pair to the items inside it; a pair found from either
  end, deduped once; a symmetric tile transposed equals the tile asked directly; held and region pairs give equal
  bins); drawing against the data (every bin pair in view has its arc and every arc is hit; every cell in view hit as
  its bin; across a 2D tile edge; a long deletion solid across its SNPs); EC-62, EC-169 to EC-173, EC-206, EC-208,
  EC-220 to EC-223 named by tests. Performance, each with its path: the pairs pyramid of 12,000 pairs built under
  10 ms (quick representative) and of 600,000 (the `summaries.maxBytes` bound) within 300 ms of slices (heavy, `itWithin`); items read over a 10 Mb view in at most `maxTilesPerWant` tiles (quick); the hysteresis holding across a pan at the threshold (quick); a region grid's pan frame under 8 ms
  with its tiles built, a tile's image within `grid.buildBudget` (quick: 4 × 4 tiles); the reference view's 2D bytes
  with its margin ring under 44 MB (quick, counted); a whole-genome Hi-C at the fit read through the whole-assembly level in one want (quick, generated); a ring of 40 matrix states within `rings.imageRuns` under 8 ms (quick); a whole genome
  Hi-C by region (heavy, by block pair). Sizes: the shapes, keys and summaries about 6 kB in headless rendering, the
  matrix layer about 5 kB; measured at RG's gate.
- **`_spec` sections that change**: core Data (Protocol, Tables, Summaries, Store), Language (pairs, grid, matrix),
  Resolve (access plans by source), Render (Data, The pairs layer, The grid layer, a new The matrix layer, Rings
  layout's images, Lanes and Hits' `bin`, Configuration).
- **Other designs amended**: [protocol](../../../design/protocol.md) (the shapes, `y`, the 2D key, the level rules by
  shape); [grids-and-rows](../../../design/grids-and-rows.md) (Data: the deferred list built; Pairs: "read by region
  … wait" replaced; Grid data: region grids; Tensions: "No navigator for a grid's own y window" answered with GA's
  square); [data-store](../../../design/data-store.md) (pairs and grid summaries); [layouts](../../../design/layouts.md)
  (Rings: an image as runs, `seq.render.ring.image` for grids only); [files](../../../design/files.md) (held pairs from
  files: `pairs.source` read by region, so its option (a) is not needed for pairs);
  [seq-data-ports](../../../design/seq-data-ports.md) (a host's adapter serving the new shapes; the conformance
  suite's new cases for adapters); [performance](../../../design/performance.md) (the structures table: pairs gain a bin-pair pyramid, a matrix row
  as cells' mirror, grids 2D tiles; Scale: a Hi-C over a genome); [features](../../../design/features.md) (Contact
  maps, Annotations' matrix, Keyed tables' data part).

## Where this differs from GA's draft

- **The square navigator's picture** (GA: "reads RR's 2D summary of the whole at its coarsest level"): for pairs there
  is no 2D summary; the square reads the pairs level's coarsest tiles **by block** (one want a block: a pairs tile holds
  every pair with an end in it), so its cost is the blocks, not their square, and sums the bins into one image of
  block × block cells, reduced to its pixels by the shared reduction. For a grid source it reads the coarsest 2D tiles
  by block pair, the upper triangle when symmetric, at most `navigator.square.maxCells` (4,096) block or run pairs:
  agreed with GA's bound, which then bounds the wants, not only the lattice's lines.
- **Small blocks**: GA's run "keeps every block's data drawn"; the first draft here read nothing under one device
  pixel. The joint rule above (runs at `blocks.runPx`, 2 CSS px; data read while at least `blocks.readPx`, 1 device
  px; pairs through their other end) replaces both, for GA's draft word for word. The measure stays GC2's.
- **Pairs across blocks**: agreed, no new owner; this draft adds that a bin cut by a block's edge is clipped to each
  block, and that a bin's strand carries GA's twist.

## Gaps for the author

Each a yes or no; the recommendation first.

1. **Pairs tiled in 1D by either end** (not 2D tiles; a runtime indexed on one end reads whole at open): yes.
2. **The pairs level is bin-pair counts** (by strand, by `type` when listed, a field summed), in the protocol and the
   store's summaries: yes.
3. **Items gated by the counts of the bins in view** (`pairs.maxArcs`), no "load anyway": yes.
4. **Grids in 2D tiles of 256 bins a side at square levels**, symmetric ones served once: yes.
5. **A grid's level: the finest whose bin is at least a device pixel** (the other way from values): yes.
6. **One value a grid cell**, coarse levels by the source's `sum` or `mean`, normalizations as `field`: yes.
7. **Held pairs get the same bin levels** (`binPairs`, one reader over held and region pairs): yes.
8. **`rows` as keyed records only** (the shape and its key index); the join (EC-224) designed elsewhere (GD's views or a
   themes amendment): yes.
9. **The `matrix` shape and a `matrix` layer** (items × states, states whole, on the rows layer's band) built in RG,
   placed in Gate 0's cut order after contact maps (MAF, contact maps, the matrix layer, the liftover): yes.
10. **`.hic` and `.mcool` runtimes**: not in the plan's formats; the Hi-C example generated in the `seq` format with two
    levels, a real format left to a later slice: yes.
11. **Images on rings as sector runs** (`seq.render.ring.image` kept for grids only): yes.
12. **Core's pairs layer draws a bin pair as one arc stroked by count** (ribbons RN's): yes.
13. **The joint small-blocks rule with GA** (runs under 2 CSS px; data read while a block is at least one device pixel;
    pairs to an unread block through their other end), in both drafts word for word: yes.
14. **A whole-assembly coarse level** (Juicebox's "All": a `whole` system in the coverage, built by the summaries for a
    one-level source) and a bound on a layer's wants together (`maxTilesPerLayer`, 1,024): yes.
15. **Bounds in bytes for every shape** (`maxTileBytes` and `summaries.maxBytes`, 16 MB each: a quarter of the cache),
    counts staying only as the parse's guards: yes.

## Review, unresolved

The independent review (2026-11-21) found two blocking issues (items unreadable past `maxTilesPerWant` at the finest
level, byte bounds counted in cells) and eight to fix; the revision took them all, and the joint small-block rule is now
in GA's draft word for word. Still open after it, for the gate:
- No `max` across a region grid's levels.
- A matrix's states are read whole: about 4,000 float or 16,000 code states at most.
- A source indexed on one end only is read whole.
- Images on rings drawn as runs are not measured yet.
- A `.hic` without a whole-assembly level still makes up to 1,024 requests.
- Reads inside runs wait for GC2's measure.
- GF's revision (genome-formats.md, "Where this differs…") asks for text here: sided pairs sources, items at tile size for every list shape, item extents and span rows with `measure`, genotype-class `codes`; taken at the gate.
