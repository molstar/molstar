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

- **Summaries cost a fraction of the data**: a pyramid of 4× steps adds at most a third of its base; summaries are
  built once at load for inline data (in a worker for large data), or served by the source when its format has them.
- **What is tested**: frame time is flat across the scale fixtures (ten times the data, at most 1.2 times the frame
  time on the same view size), and proportional to the view size (twice the pixels, at most twice the time); a
  region query's cost grows with log n, measured as a ratio.

## Measured first

The numbers below are targets to confirm, not promises: each is measured on the first end-to-end slice
([version-1.md](version-1.md)) and then confirmed or moved by the author's decision. Two clarifications from the
review:

- **Draw lists are built into reused buffers**: a frame allocates nothing once warm, and the draw lists are pure by
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
- **Limits are configurable** ([extensibility.md](extensibility.md)): every cap below is a setting with this default.

## Time

| What | Target |
|---|---|
| Hover: hit-test, projection onto linked sections, overlay redraw | under 2 ms; the frame under 8 ms |
| Pan and zoom frame (letters, tracks, at the scale sizes below) | under 8 ms (120 Hz with headroom); never over 16 ms |
| First zoom into new letters (glyph atlas filling) | under 8 ms a frame |
| Selection drag, keyboard cursor step | under 4 ms a step |
| Validate and resolve a spec of 1,000 nodes (inline data aside) | under 20 ms |
| First paint of a dataset view (p53, globins family) | under 100 ms after the data is in memory |
| An edit: apply a patch, resolve again, redraw | under 16 ms; undo the same |
| Pan and zoom frame at every scale fixture | the same targets: frame time does not grow with the data |
| Open a 10,000 × 300 alignment (summaries built) | under 300 ms to the first frame |
| Open a 3M-nucleotide sequence from a file, with tracks (summaries built) | under 1 s to the first frame; then frames as for any size |
| A region arriving from a source to drawn | within the next frame after it is decoded |
| bigWig: a region at the right zoom level, local file | under 50 ms to decoded |
| Export a view to SVG (in the browser) | under 200 ms |
| Render CLI: a dataset view to SVG, cold start | under 1 s |
| Decode MessagePack | at least 200 MB/s; typed columns without a copy where aligned |

## Memory

| What | Target |
|---|---|
| Draw and hover loops, warm | 0 bytes allocated a frame (guarded) |
| Resolved model of inline data | at most 2× the data's MessagePack size, multi-resolution summaries included (a third of their base at most) |
| Draw lists | bounded by the screen (device columns × lanes), never by the data |
| Maps | never expanded past the region asked; a protein's path through a chromosome costs O(its coding sequence) |
| Lazy alignment rows (10M) | O(rows in view + rows bound to maps), never O(all rows) |
| Source cache | byte-bounded, 64 MB by default; what is in view is never evicted |
| Glyph atlas | bounded, 16 MB of pages by default; least recently used page evicted, never one in use |
| A grid or matrix block | at most 64 MB (2^24 floats) a block pair |
| A viewer retired | releases everything: no heap growth over n and 4n open and close cycles |

## Size

| What | Target (minified) |
|---|---|
| Reading a spec (decode, validate, resolve), no viewer | under 150 kB |
| Headless rendering (the above, layout, draw lists, SVG) | under 300 kB |
| The viewer, default configuration | under 600 kB (the previous one: 1,040 kB) |
| A format runtime (bigWig) | under 30 kB, loaded only when composed |

## Scale

The sizes every target above holds at, each an example or a generated fixture: a 35,000-residue sequence; a 3M-nucleotide
sequence from a file; a whole human genome by region; alignments of 10,000 × 300 and 1,000 × 2,000; 10M lazy rows; 200
copies of one chain; 100 tracks in one view; a 4,000 × 4,000 grid.
