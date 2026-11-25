# M5 measurements (fifth slice, 2026-11-23, `packages-seq` at `6f8f5b768`)

The brief is [briefs/m5-measurements.md](briefs/m5-measurements.md); it ran on Opus, alone. The setup:
- Chrome in the Browser pane, at an emulated 1440 × 900, **dpr 2** (confirmed in the page), with the page visible.
- Load 1.4 before the browser block. No console errors.

**Method**: the fourth slice's ([M4](../fourth-slice/measurements.md)).
- `requestAnimationFrame` is replaced by a queue flushed right after each synthetic input.
- 20 warm frames, then 300 measured, timed with `performance.now()` (0.1 ms resolution).
- These are main-thread numbers. Raster and compositor time show only in the real-vsync run.

| Case | Story | Frames | Median ms | p95 ms | Max ms | Target | Verdict |
|---|---|---|---|---|---|---|---|
| ring, 12,000 links: hover stems / ribbons / wheel pan | Features/Genome Ring: The Genome At The Fit | 300 each | 0.1 / 0.0 / 0.1 | 0.2 / 0.1 / 0.2 | 0.3 / 0.2 / 0.4 | hover < 2, frame < 8, never > 16 | pass |
| the square's window drag | same | 300 | 0.1 | 0.7 | 1.8 | same | pass (release 0.3 ms) |
| two layers: hover stems / ribbons / wheel pan | Kinds Of Links | 300 each | 0.1 / 0.0 / 0.1 | 0.2 / 0.1 / 0.2 | 0.6 / 0.2 / 2.3 | same | pass |
| bridge: hover / pan / zoom | Genomes/Views: Synteny | 300 each | 0.0 / 0.5 / 0.8 | 0.2 / 0.9 / 1.9 | 0.3 / 2.5 / 5.6 | frame < 8 | pass |
| bridge: hover / pan / zoom | Liftover | 300 each | 0.1 / 1.3 / 1.6 | 0.6 / 2.0 / 3.3 | 1.4 / 2.8 / 4.4 | same | pass |
| GRCh38 by chromosome: pan across the whole genome / at 28 kb | Features/Genome Coordinates: Cytobands Over The Genome | 300 each | 0.3 / 0.3 | 0.8 / 0.4 | 3.4 / 1.1 | same | pass |
| by region (tabix bedGraph): pan / zoom | Genomes/Views: Signal | 300 each | 0.2 / 0.6 | 0.7 / 0.8 | 1.4 / **25.7** | never > 16 | **zoom: one frame over** |
| by region: a jump to drawn | same | 20 | 3.6 | 5.1 | 6.1 | | pass |
| bigWig: pan / zoom / hover | p53/Signal from a bigWig | 300 each | 0.2 / 0.1 / 0.1 | 0.4 / 0.2 / 0.2 | 0.8 / 0.2 / 0.7 | frame < 8 | pass (M4's gap closed) |
| bigWig: pending to drawn | same | 8 | 170 | 171 | 171 | | the story's own 200 ms delay; 12 of 20 jumps were drawn in the move's frame |
| 200 copies, split patch: click to drawn | Complex/Views: Copies, then Split | 5 | **16.2** | | 20.4 | edit < 16 | **over**: 12.8 / 14.1 / 16.2 / 17.2 / 20.4; the final picture in the first frame (M4's one sample, 33.6, re-measured) |
| 200 copies, merge patch | same | 6 | about 6 | | 7.3 | edit < 16 | pass |

**Real vsync** (real rAF, pane in view, The Genome At The Fit): no dropped frames (none over 33.4 ms).

| Run | Median interval ms | p95 ms | Max ms | Intervals over 17.2 ms |
|---|---|---|---|---|
| Hover stems | 16.7 | 17.8 | 18.7 | 22 |
| Hover ribbons | 16.7 | 18.1 | 18.8 | 29 |
| Wheel pan | 16.7 | 17.2 | 18.6 | 14 |
| Idle baseline | 16.7 | 17.2 | 18.6 | 13 |

**The VCF, 10,000 records**: measured in Node, by its heavy-lane bound (`vcf-time.test.ts`).
- Sites decode under 100 ms; the whole test took 74 ms.
- 100 samples' genotypes decode under 200 ms; the whole test took 1,089 ms, mostly generating the data.
- The browser can't read the fixtures (`node:zlib`), and no story holds 10,000 records.

**The Node tiers**, one at a time:

| Tier | Load before | Result |
|---|---|---|
| Heavy, `packages/seq/` | 1.21 | 174 of 176 passed |
| Core perf | 2.30 | 245 of 245 passed |
| Storybook perf | 1.92 | 22 of 22 passed |

- The heavy lane's two failures were both timing bounds: the 1M-range interval lookup (0.71 against 0.62 ms) and the pairs pyramid's slices (22–51 against 20 ms). Rerun alone, both files passed. The gate's own run at load 2.5 also passed them.
- So these bounds fail under the lane's own worker contention, not under machine load.

**Sizes** (`size`): every row passes.

| Entry | kB now | Bound | Fourth slice |
|---|---|---|---|
| Reading a spec | 260.2 | 265 | 223.4 |
| Headless rendering | 974.3 | 980 | 810.3 |
| Default viewer | 1,539.7 | 1,545 | 1,304.2 |
| Minimal viewer | 1,086.9 | 1,090 | 936.2 (+150.7) |
| The ring | 43.9 | 45 | |
| Each runtime | 9.3 to 23.1 | 30 each | |
| Bridge chunk (lazy) | 25.4 | 28 | |

The plan said "the minimal viewer not grown". It grew 150.7 kB over the slice: the genome axis, the navigator over blocks, region reads, bridges and the card. Each step raised its bound with the step's reason (in the design `performance.md`).

**Gaps**:
1. No whole human genome is read by region: the Genomes views are toy genomes, and GRCh38's cytobands are read whole.
2. Cmd or Ctrl+wheel over the ring at the fit is taken but shows no change (novice 4, FX24: the wheel pans).
3. The VCF decode is bounded in Node, not timed in a browser.
4. Raster and compositor time show only in the real-vsync intervals.

**Placed**:
- The split patch over 16 ms → TODO.
- The Signal zoom's one frame of 25.7 ms → TODO.
- The two heavy bounds that fail under the lane's contention → TODO: give them `tries` or run them isolated, not raise them.
- The minimal viewer's growth → a question for the author.
- Cmd+wheel over the ring at the fit → TODO.
