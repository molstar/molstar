# Drift: the zoom strip laid out again on every hull change of a drag

**Level**: code-vs-design. **Proposed outcome**: undecided.

**What**: `performance.md` (Cost follows pixels, the drag rule) bounds a drag's work by the screen, and the plan's CL11
item 3 asks the navigator's zoom strip not to be laid out again per move. The design (core `_spec` "Render",
Navigators, A navigator over blocks) makes the strip's viewport the hull of the blocks the window touches; a window
dragged across 3,000 contigs touches another hull almost every frame, and a new viewport is a new layout of the strip's
picture (the renderer's frame key holds the viewports). Measured 97 kB a frame (`navigator-blocks.quick.ts`, bounded at
256 kB), 0.7 ms headless.

**Where**: `core/src/render/navigator.ts` (`frame`: a new `reach` sets the strip's viewport and frames again),
`navigator-base.ts` (`placeBlocksX`, `setViewports`), `renderer.ts` (`frame`: any viewport change is a `layoutView`),
`layout/tracks.ts` (`layStrip`, `layOne`), `ruler.ts` (`layout`).

**Evidence**: a sampling heap profile of the case, a frame (CL11, 2026-11-24): the strip's ruler layouts 38 kB (a new
lay for each block coming into view: its tick arrays and label strings; starting the arrays at 8 instead of 64 saved
1 kB of the 97), `layStrip` 15 kB (a `LayoutContext`, a box and a lays array a strip a layer, every frame, for blocks
seen before too), `layoutView` itself 15 kB, the section cache's replayed `Set` adds 10 kB, `placeBlocks` 5 kB, a regex
`exec` 6 kB. No one record dominates: reusing the ruler's lays across blocks (the route that keeps the span rule)
recovers the arrays only, the labels are new strings per block, and the rest is the generic tracks layout's
per-strip allocation, which every layout of a section of many blocks pays.

**Options**: (a) keep the span rule and make the tracks layout reuse its per-strip records (a pool of `LayoutContext`s,
boxes and lays arrays by strip index; the ruler's lays pooled by layer type; label strings cached by value and step):
a layout-wide refactor of `layoutView`, worth it for every section of blocks (a genome row panned across contigs pays
the same), not a drag lane's; (b) change the span rule during a window's drag: the strip keeps its span while the
touched blocks lie inside it and takes the touched hull widened (a few windows' width) when they leave, so it is laid
out a few times across a sweep, the exact hull again on release; (c) leave it: 97 kB and 0.7 ms a frame, bounded by
the screen (twenty contigs), never by the genome's 3,000 blocks.

**Recommendation**: (c) for version 1, with (a) in the TODO as the layout's allocation pass (it also lowers a genome
row's pan); (b) only if a measured frame on the author's machine misses the budget, since it changes what the strip
shows mid-drag.
