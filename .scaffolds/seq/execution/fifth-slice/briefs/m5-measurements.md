# M5: the fifth slice's measurements

**Goal**: the plan's stage 5, "M5", in the shape of the fourth slice's
[measurements.md](../../fourth-slice/measurements.md) (its method, its table, its gaps). You measure and report. You
change no code, and you write nothing in the repository: the orchestrator writes `measurements.md` from your report.

**Where**: the main checkout at `packages-seq` (read-only), and Storybook at http://localhost:6006, already running.
Never start or stop a server. Use the built-in browser's tools (`mcp__Claude_Browser__*`). Call `tabs_context` first.
A hidden pane is never evidence: if it says hidden, wait and recheck, and stop and say so if it stays hidden. The
machine must be quiet: take `uptime` before each block and wait for a 1-minute load under 3. Nothing else of the
orchestrator's runs beside you.

**Method**: the fourth slice's. Chrome at an emulated 1440 × 900 at dpr 2, confirmed in the page. Replace
`requestAnimationFrame` by a queue flushed right after each synthetic input. Each frame is the input plus the draws it
queued, timed with `performance.now()`. Warm 20 frames, then measure 300, and give the median, p95 and max. Then, for
the real-vsync cases, a second run with the real `requestAnimationFrame` and the pane in view. Report those frames'
intervals, and the dropped frames over 16.7 ms.

**Cases**:

1. **The ring of 12,000 links at full scale, at dpr 2** (Features › Genome Ring: at the fit, and Kinds of links, two
   layers): hover over ribbons and stems, a wheel pan, and a drag of the square's window. Targets: hover under 2 ms, a
   frame under 8 ms, never over 16 ms.
2. **A bridge frame** (Examples › Genomes › Views › Synteny, and Liftover): hover, pan and zoom.
3. **A whole human genome by region, by chromosome**: the Genomes views that read by region (the Index lists them),
   panned across chromosomes. Report the frames, and the time from a move to the region drawn (pending to drawn).
4. **A VCF region of 10,000 records decoded**: say which story or test holds 10,000 or more records. Time the decode
   in the browser, or in Node through the package's own tests or perf cases (`pnpm -s -C packages/seq/core test:perf`
   lists them), and say which you used.
5. **M4's carried measurements**: a bigWig pan (the story FX1 added: find it in `/index.json`); real-vsync frames with
   the pane in view, on one of the cases above; the 200-copy split patch, from one click to drawn (Examples › Complex ›
   Views › Copies, then Split), five samples, not one.
6. **The tiers in Node, alone and quiet**: `pnpm test:heavy packages/seq/`, `pnpm -s -C packages/seq/core test:perf`,
   `pnpm -s -C packages/seq/storybook test:perf`, and `pnpm -s -C packages/seq/core size` (every row: each runtime's own
   row, the ring's own, and the minimal viewer's growth since the fourth slice's record). Take the load before each.

**Report** (under 900 words): the table in the fourth slice's shape (case, story, dpr, frames, median, p95, max,
target, verdict), the Node tiers' numbers, the load at each block, the gaps (what you could not measure and why), and
anything over its target with the case named. Then set the viewport back to desktop.

**Model**: strong. Read-only: no edits, no commits, no servers. Write temporary files only in your own scratch space.
