# The close-out's measurements

**Goal**: the plan's measurements (the bundles against their budgets, the representatives at full size), in the shape
of the fifth slice's [measurements.md](../../fifth-slice/measurements.md) (its method, its table, its gaps), with the
close-out's new views. You measure and report. You change no code and write nothing in the repository: the orchestrator
writes `measurements.md` from your report.

**Where**: the main checkout at `packages-seq` (read-only), and Storybook at http://localhost:6006, already running.
Never start or stop a server. Use the built-in browser's tools (`mcp__Claude_Browser__*`). Call `tabs_context` first;
open a tab of your own in the foreground (`tabs_create`, `foreground: true`) and keep it in front. A hidden pane is
never evidence: if it says hidden, wait and recheck, and stop and say so if it stays hidden. The machine must be quiet:
take `uptime` before each block and wait for a 1-minute load under 3. Nothing else of the orchestrator's runs beside
you.

**Method**: the fifth slice's. An emulated 1440 × 900 at dpr 2, confirmed in the page. `requestAnimationFrame`
replaced by a queue flushed right after each synthetic input; each frame is the input plus the draws it queued, timed
with `performance.now()`. Warm 20 frames, measure 300: median, p95, max. Then one real-vsync run (real
`requestAnimationFrame`, the pane in view): the intervals and the dropped frames over 16.7 ms.

**Cases**:

1. **The fifth slice's cases again**, to compare: the ring of 12,000 links (Features › Genome Ring › The Genome At The
   Fit; Kinds Of Links): hover stems and ribbons, a wheel pan, the square's drag; a bridge (Examples › Genomes › Views ›
   Synteny, Liftover): hover, pan, zoom; GRCh38 by chromosome (Features › Genome Coordinates › Cytobands Over The
   Genome): pan; by region (Genomes › Views › Signal): pan, zoom, a jump to drawn; the 200-copy split patch (Complex ›
   Views › Copies, then Split), five samples. Its two overs (Signal's zoom, one frame at 25.7 ms; the split patch at
   16.2 ms median) named again with today's numbers.
2. **The close-out's views**: search (Examples › p53 › Search: a run of "PXXP" from typing to drawn, Next and Previous
   steps); track groups (Features › Track groups › Opening: a fold and an open, click to drawn; a pan with groups
   folded); reference lines (Genomes › Views › Manhattan: a pan and a hover); the bridge after one side moves
   (Features › Across sections › One Side Zoomed In: a pan of the zoomed side, FX28's per-view prep); the legend's
   layers menu (Features › Legends › Legend Placements: open to drawn); `height: 'content'` (Features › Tooltips › By
   Data: a pan with every lane drawn).
3. **The tiers in Node, alone and quiet**: `pnpm test:heavy packages/seq/` (after `npx tsc -b`),
   `pnpm -s -C packages/seq/core test:perf`, `pnpm -s -C packages/seq/storybook test:perf`, and
   `pnpm -s -C packages/seq/core size` (every row, and each bundle's growth over the close-out: the fifth slice's record
   in its `measurements.md`). Take the load before each.

**Targets**: hover under 2 ms; a frame under 8 ms, never over 16 ms; an edit (a patch, a fold) under 16 ms to drawn; a
search run under 50 ms to drawn for a protein and its gene.

**Report** (under 900 words): the table in the fifth slice's shape (case, story, dpr, frames, median, p95, max, target,
verdict), the Node tiers' numbers, the load at each block, the gaps (what you could not measure and why), and anything
over its target with the case named. Then set the viewport back to desktop and close your tab.

**Model**: strong. Read-only: no edits, no commits, no servers. Temporary files only in your own scratch space.
