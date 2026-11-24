# Performance targets

## Intent

The author, 2026-11-17: "let's add "performance" page to the storybook with our targets? both time and memory wise".

Targets set before the code, so the design is held to them rather than measured after the fact (§performance: designed
in, guarded by tests).

## Provisional

Every number below is a first proposal, drawn from what the previous implementation measured (a 3M-nucleotide genome
loaded under 1.5 s with frames under 16 ms; first-zoom frames of 5 to 8 ms at DPR 2 with its glyph atlas; a viewer of
about 720 kB minified from an empty spec and 1,040 kB by default; grids capped at 64 MB a block pair). The author
confirms them; a target moves only by a recorded decision.

## Cost follows pixels

The author, 2026-11-17: "what is rendered depends on the pixels on the screen, not the amount of data ... so our data
structures should be chosen appropriately".

- **A frame costs what is on screen**: O(device pixels × visible layers, plus the letters and labels that fit), not
  O(data). A 3-billion-base genome and a 300-residue protein draw at the same frame time; only the levels of detail
  differ. Data size costs only at load (building summaries, once) or, for sources read by region, not at all.
- **Draw lists are bounded by the screen**: at most one primitive a device column a lane (a bar, a density cell, a
  run), or one a cell when cells are wide enough for letters; never one a data item when items are thinner than a
  pixel.
- **Every layer reads from a multi-resolution structure**, chosen per data shape so a query for a region at a
  resolution costs O(log n + what it returns at that resolution):

  | Data | Structure |
  |---|---|
  | values | a summary pyramid (min, max, mean, count a bin, each level 4× coarser), the base the typed column; bigWig's own zoom levels when read from one |
  | ranges, variants, points | sorted by start with a running maximum of ends (an interval index over typed columns), plus count pyramids by category for density |
  | letters | read only for the cells in view; at coarser levels, a class or consensus pyramid |
  | alignments | rows virtualized (only rows in view are touched); column summaries (consensus, class, gap fraction, conservation) per tile, pyramided |
  | matrix | per state, bins pyramided by majority or maximum; states virtualized past a device row each |
  | grid, dot plots | a 2D pyramid (mean, maximum), like mipmaps |
  | pairs | indexed by both ends; arcs culled by their span against the view |
  | maps | runs searched by binary search; paths precomposed at compile; projection O(log runs) |
  | text layout | lines laid out lazily, only the lines in view, with their offsets kept for scrolling |

- **A drag's work is bounded by the screen, not by the selection** (the author, 2026-11-21: selecting a rect on the
  tRNA dot plot "is really laggy", "similar lag to selecting on the genome ring with the contigs"; the close-out's
  CL11): a drag-select writes the selection at most once a frame (its moves coalesced to the viewer's frame, the full
  selection resolved on release, one history step), hovers nothing while it runs (no hit, no hover lookup, no card),
  derives its emphasis and projections in proportion to what is shown (a projection walks a part's links, never every
  system of a 3,000-block genome), costs nothing while the pointer is held (no hit, no value made), and leaves what
  follows the selection outside the frame (a story's SVG twin) to the release. Why: a pointer sends moves faster than
  frames, and each value written re-ran everything derived from the selection; measured before, a sweep across 12
  chromosomes of the 3,000 contigs drew about 300 ms a frame (the projection's walk), the 4,000 × 4,000 dot plot's SVG
  twin exported 0.5 to 1 s a move.
- **Summaries cost a fraction of the data**: a pyramid of 4× steps adds at most a third of its base; summaries are
  built once when a store entry is loaded (on the main thread in time slices until the compute port exists), or served
  by the source when its format has them; they belong to their store entry and are evicted with it.
- **What is tested**: frame time is flat across the scale fixtures (ten times the data, at most 1.2 times the frame
  time on the same view size), and proportional to the view size (twice the pixels, at most twice the time); a
  region query's cost grows with log n, measured as a ratio.

## Measured first

The numbers below are targets to confirm, not promises: each is measured on the first end-to-end slice
([version-1.md](version-1.md)) and then confirmed or moved by the author's decision. Two clarifications from the
review:

- **Draw lists are built into reused buffers**: drawing and hovering allocate nothing once warm (a layout's allocation
  is bounded by the screen, not zero: the author, 2026-11-18), and the draw lists are pure by
  contract (same inputs, same lists), not by being new objects each frame.
- **Ratios are checked where they cannot flake**: the scaling tests run in the heavy lane, each case repeated and its
  median taken, with the ratio's tolerance stated per test.

## Rules

- **The reference machine**: a mid-range laptop of 2024 (an Apple M2 or a recent x86 with integrated graphics), Chrome
  stable, DPR 2, 1440 × 900 CSS pixels; Node 22 for headless numbers. Time targets are for that machine; scaling
  targets (how time grows with size) are machine-independent ratios.
- **How each is held**: time bounds in the heavy lane (`itWithin`, `test-harness`), scaling as ratios (4× the data
  costs at most 4×, or the stated order), allocation for the hot loops only, by hand at checkpoints (the opt-in tier, [testing.md](testing.md)), heap growth by
  the heap guards at checkpoints, bundle sizes by a size check on the built entries. A target without a test is a wish.
- **A target names its path and its test when it is written**: the code path it covers (the draw loop, the layout, the
  hover, a format's decode) and the test that holds it. Why: a target without a path is read more than one way ("a
  frame allocates nothing" was read as the draw loop and as the layout too), and one without a test is not measured
  until late.
- **A quick representative when the path is built**: the stage that builds a path runs one scenario of its target
  with one bound, before later stages build on it; the full target (every scale fixture, the medians, the ratios)
  stays in the heavy lane. Why: measured only at the end, a miss is found under everything built on it.
- **Limits are configurable** ([extensibility.md](extensibility.md)): every cap below is a setting with this default.

## Levels and a frame's cost

A frame costs what is on screen only if the level a layer reads fits the screen. A source's levels are a factor `f`
apart (the store's summaries 16; a format's own, such as bigWig's zoom levels, whatever it has). Reading the finest
level no coarser than a device pixel, a layer can read up to `f` items a device column, so between levels the frame
costs up to `f` times the pixels.

- **The level choice is a config item**: a layer may read a level up to `k` times coarser than a device pixel per bin,
  so it reads at most about `f / k` items a column and its resolution drops by at most `k`. `k` defaults to 2; 1 keeps
  the finest level no coarser than a pixel (the author, 2026-11-18).
- **Placing the items of a fine level is itself a fast path**: on an identity or linear map a values layer places
  positions and aggregates per device column without mapping each position, and gives the same draw lists as the
  general path (a test compares both).
- **What is tested**: frame time at a density between two levels, not only at a level's own density (the 3M bars
  case: a synthetic source with levels 64× apart, about 38 positions a device pixel). Why: the miss came from the gap
  between levels, which a density at a level's own never shows.

## Time

| What | Target |
|---|---|
| Hover: hit-test, projection onto linked sections, overlay redraw | under 2 ms; the frame under 8 ms |
| Pan and zoom frame (letters, tracks, at the scale sizes below) | under 8 ms (120 Hz with headroom); never over 16 ms |
| First zoom into new letters (glyph atlas filling) | under 8 ms a frame |
| Selection drag, keyboard cursor step | under 4 ms a step; a drag-select's frame (its moves and the frame applying them) under 4 ms for each drag kind: a rect on the 4,000 × 4,000 dot plot, a range on its ruler, a box on the PAE grid, a sweep across 12 chromosomes of the 3,000 contigs on the ring and the row (`drag-time.test.ts`, heavy) |
| Validate and resolve a spec of 1,000 nodes (inline data aside) | under 20 ms |
| First paint of a dataset view (p53, globins family) | under 100 ms after the data is in memory |
| An edit: apply a patch, resolve again, redraw | under 16 ms; undo the same |
| Pan and zoom frame at every scale fixture | the same targets: frame time does not grow with the data |
| Open a 10,000 × 300 alignment (summaries built) | under 300 ms to the first frame |
| Open a table of 100,000 records from a file (TSV, JSON or MessagePack, the transform applied) | under 300 ms to its tables ([files.md](files.md)); 10,000 rows the quick representative |
| Open a 3M-nucleotide sequence from a file, with tracks (summaries built) | under 1 s to the first frame; then frames as for any size |
| A region arriving from a source to drawn | within the next frame after it is decoded |
| bigWig: a region at the right zoom level, local file | under 50 ms to decoded |
| Export a view to SVG (in the browser) | under 200 ms |
| Render CLI: a dataset view to SVG, cold start | under 1 s |
| Decode MessagePack | documents of columns at least 200 MB/s, typed columns without a copy where aligned; node-heavy documents (20,000 small nodes) at least 20 MB/s with validation, a regression guard (validation is not skipped; the author, 2026-11-18) |

## Memory

| What | Target |
|---|---|
| Draw and hover loops, warm | 0 bytes allocated a frame (guarded) |
| A drag-select's frame | 0 bytes while the pointer is held (no hit, no value made; the genome views' frame its own 48 B, not the drag's); a frame moving the drag's end its new value and projections, under 16 kB (`drags.quick.ts`, perf tier) |
| A hover through the structure link | at most 2 kB a lit unit (the target and its Location are new each hover; pooled if a large assembly shows it; the author, 2026-11-21) |
| Layout, warm, at a fixed screen | bounded a frame, independent of the data's size (guarded); not zero |
| Resolved model of inline data | at most 2× the data's MessagePack size, multi-resolution summaries included (a third of their base at most) |
| Draw lists | bounded by the screen (device columns × lanes), never by the data |
| Maps | never expanded past the region asked; a protein's path through a chromosome costs O(its coding sequence) |
| Lazy alignment rows (10M) | O(rows in view + rows bound to maps), never O(all rows) |
| Source cache | byte-bounded, 64 MB by default; what is in view is never evicted |
| Glyph atlas | bounded, 16 MB of pages by default; least recently used page evicted, never one in use |
| A grid or matrix block | at most 64 MB (2^24 floats) a block pair |
| A viewer retired | releases everything: no heap growth over n and 4n open and close cycles |

## Size

Bounds, not guides (the author, 2026-11-20: "reduce render first, then move the size targets", "give the bounds a
little wiggle room"): the third slice's sizes after its render reduction, with about 10 % room, checked at every
gate; a growth past one is fixed or brought to the author with the feature behind it. Headless rendering and the
default viewer hold every built-in on purpose: their size is the features; a smaller host composes the minimal viewer.
`"sideEffects": false` for `@molstar/core` and `@molstar/ui` (about 17 kB) waits for an audit of their import-time
effects (`TODO.md`).

| What | Target (minified) |
|---|---|
| Reading a spec (decode, validate, resolve), no viewer | under 265 kB (260 until TT's tooltip declaration, the orchestrator under the soft-size rule; measured 260.0) (240 until RG's region reads, the author, 2026-11-22: "Ok increase the limit"; measured 254.9 with it; 150 in the first slice, 160 in the second, 195 at the third slice's close, 230 for language 0.5.0's molecules; moved by the author, 2026-11-21, for language 0.6.0's genomes: measured 223.4 before it) |
| Headless rendering (the above, layout, draw lists, SVG) | under 980 kB (970 until TT's tooltips, the orchestrator under the soft-size rule, the tooltips behind a 26.7 kB lazy door; measured 968.1; 960 until LA4's liftover and FX17's matrix emphasis, raised by the orchestrator under the soft-size rule; measured 960.4 with them; 950 until LA3's grids over two sections and categorical axes, raised by the orchestrator under the soft-size rule, the members behind a lazy door; measured 950.8 with it; 940 until SX's section scroll, raised by the orchestrator under the soft-size rule after a trim pass recovered 0.1 kB; measured 941.2 with it; 930 until RL's rail and LA's bridges met, the author, 2026-11-22: "Raise to 940"; measured 930.5 with both; 920 until FX9's track sizing met GA3, the author, 2026-11-22: "You can raise the limits"; measured 922 with both; 910 until RG's region reads met the merged lanes, the author, 2026-11-22: "Raise again"; measured 912.1; 870 until RG's region reads, the author, 2026-11-22; measured 897.2 with it; 830 until GC2's genome axis, the author, 2026-11-21: "raise the bound, even 10% is still fine"; measured 833.6 with it; 300 until the third slice, whose layers, layouts, navigators and rows it holds; 780 until the fourth slice's query language, 820 until molecules drawn, the author, 2026-11-20; measured 820.9) |
| The viewer, default configuration | under 1,560 kB (1,550 until the close-out's CL9, the orchestrator under the author's leave (2026-11-23: "you can increase the budget if needed"): the codon snap and the version-aware validate in the viewer; measured 1,550.7; 1,545 until the close-out's CL1 and CL4, the orchestrator under the soft-size rule: the legend's eye already lazy, the row names, the genotype calls and the hidden rows eager; measured 1,545.2; 1,530 until TT's tooltips and card providers, the orchestrator under the soft-size rule; measured 1,529.9; 1,510 until HC's hover contract, FX16's rows zoom and FX17's matrix emphasis, raised by the orchestrator under the soft-size rule; measured 1,511.2 with them; 1,495 until LA3's grids over two sections and categorical axes, raised by the orchestrator under the soft-size rule; measured 1,500.5 with it; 1,480 until VL's values look C and SX's section scroll, raised by the orchestrator under the author's soft-size rule (2026-11-22), after a trim pass found no dead code (VL); measured 1,481.3 with VL; 1,470 until LA's bridge viewer, the author, 2026-11-22: "raise default to 1,480"; measured 1,470.4 with it; 1,460 until FX9's track sizing met GA3, the author, 2026-11-22; measured 1,462 with both; 1,450 until RG's region reads met the merged lanes, the author, 2026-11-22: "Raise again"; measured 1,450.1; 1,400 until RG's region reads, the author, 2026-11-22; measured 1,428 with it; 1,330 until GC2's genome axis, the author, 2026-11-21; measured 1,332.3 with it; 600 until the third slice, 1,250 until the query language, 1,290 until the structure link, the author, 2026-11-20: "Move to 1330 kB", room for the rest of the fourth slice; measured 1,290.4; the previous Seq\*: 1,040 kB) |
| The viewer, minimal (tracks, letters and ranges) | under 1,095 kB (1,090 until the close-out's CL10, the author, 2026-11-24: "you can increase the budget if needed"; CL1, CL2, CL4 and CL10 brought it to 1,090.0; 1,090 accepted as version 1's budget, the author, 2026-11-23: "accept as is"; measured 1,086.9, 368.5 gzipped; the fifth slice's +150.7 is the genome axis, region reads, bridges and the card, about 70 of it used only by genomes or reads by region; lazy doors for them, or genomes opt in for the minimal composition, wait for a host that needs a smaller viewer; 1,070 until TT's tooltips, the orchestrator under the soft-size rule, the declaration and check eager, the rest lazy; measured 1,078.8; 1,060 until LA4's liftover, raised by the orchestrator under the soft-size rule, the liftover behind the genome chunk's lazy door; measured 1,061.9 with it; 1,040 until LA's bridge met FX10's rail, the author, 2026-11-22: "raise minimal to 1,060 too"; measured 1,039.95 with LA; 1,000 until RG's region reads, the author, 2026-11-22; measured 1,026 with it; 960 until GC2's genome axis, the author, 2026-11-21; measured 964.4 with it; 880 until the query language, 920 until molecules drawn, 930 until the extension fallback, 940 until language 0.6.0's genomes, the author, 2026-11-20 and 2026-11-21: "Raise the bounds"; measured 939.4 before it) |
| A format runtime (bigWig) | under 30 kB, loaded only when composed |
| The genome ring extension (`@molstar/seqx-ext-rondo`) | under 48 kB (45 until the close-out's CL2 and CL3: the stem menu, the hole's words, the square's hint and pending state; measured 45.0; the author, 2026-11-24: "you can increase the budget if needed"; 30 in the design; 32 after RX's fix round; 38 after RX2's square; 45 after RX3's kinds of links (colour channel, legend, per-layer scale), the legend and card lazy first; the orchestrator under the soft-size rule; measured 43.9) |

## Scale

The sizes every target above holds at, each an example or a generated fixture: a 35,000-residue sequence; a 3M-nucleotide
sequence from a file; a whole human genome by region; alignments of 10,000 × 300 and 1,000 × 2,000; 10M lazy rows; 200
copies of one chain; 100 tracks in one view; a 4,000 × 4,000 grid.
