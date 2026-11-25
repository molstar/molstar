# The close-out's measurements (2026-11-24, `packages-seq` at `f290f08a5`)

Brief: [briefs/measurements.md](briefs/measurements.md); one Opus agent alone in the pane, then the Node tiers alone.
The report as given, its tables kept; the orchestrator's reading at the end.

**Setup**
- Emulated 1440 × 900 at dpr 2, confirmed in the page (`innerWidth` 1440, `devicePixelRatio` 2, `visibilityState` visible). The pane said "displayed" throughout.
- Load before the browser block was 2.50; it read 1.74 to 2.08 during the block.

**Method**: the fifth slice's: `requestAnimationFrame` replaced by a queue flushed after each synthetic input, 20 warm frames, then 300 measured. I had to add two things:
1. **A pointer-capture stub.** Synthetic `setPointerCapture` throws, so drags never engaged until I stubbed it. Any earlier drag numbers may not have measured a drag.
2. **A microtask drain.** After the flush I drain microtasks (`await null` ×30), then flush again. Much of the work runs as microtasks after the queued draw, and the fifth's method missed it. Where they differ, cells read **old method | with microtasks**.

| Case | Story | Frames | Median | p95 | Max | Target | Verdict |
|---|---|---|---|---|---|---|---|
| ring 12k: hover stems / ribbons / wheel pan | Genome Ring: At The Fit | 300 each | 0.0 / 0.0 / 0.1 \| 1.8 / 0.0 / 1.9 | 0.2 / 0.1 / 0.2 \| 2.1 / 1.9 / 2.3 | 0.4 / 0.4 / 0.3 \| 4.0 / 3.9 / 3.9 | hover < 2 | stems median 1.8, p95 2.1: **hover at its bound** |
| square drag (moves) | same | 300 | 0.8 | 1.0 | 4.0 | frame < 8 | pass |
| **square release, ring redrawn** | same | 16 | **~36** | | 40 (24 min) | frame never > 16 | **over**: ~18 ms microtask work, then a ~17 ms draw. The fifth's "release 0.3" released where it started, so nothing changed |
| two layers: stems / ribbons / pan | Kinds Of Links | 300 each | 0.1 / 0 / 0.1 \| 1.8 / 0 / 1.9 | \| 2.0 / 1.9 / 2.0 | \| 5.4 / 5.8 / 4.0 | | stems at the bound |
| bridge: hover / pan / zoom | Synteny | 300 each | 0.0 / 0.2 / 0.1 | 0.1 / 0.3 / 0.2 | 0.2 / 0.5 / 0.3 | frame < 8 | pass |
| bridge: hover / pan / zoom | Liftover | 300 each | 0.0 / 0.3 / 0.2 | 0.1 / 0.6 / 0.3 | 0.1 / 1.3 / 0.6 | | pass |
| GRCh38: pan chr6–chr12 / chr8–9 / at 44 kb | Cytobands Over The Genome | 300 each | 0.2 / 0.3 / 0.3 | 0.3 / 0.5 / 0.4 | 0.7 / 3.8 / 2.6 | | pass (I hit 44 kb, not 28 kb) |
| by region: pan / zoom / big zoom | Signal | 300 each | 0.2 / 0.1 / 0.1 | 0.4 / 0.2 / 0.5 | 0.8 / 0.3 / 0.8 | | pass; the 25.7 ms zoom frame did not recur |
| by region: jump, first frame | same | 20 | 2.3 | 2.8 | 2.8 | | pass |
| by region: settled frame (150 ms after a move) | same | 20 | ~11 | 13.9 | 13.9 | frame < 8 | **over 8** (mostly the story's SVG twin, see below) |
| **200 copies split** | Complex/Views: Copies | 5 + 5 | **19.6** \| 49 | | 21.3 \| 50 | edit < 16 | **over, worse than 16.2**; final picture in the first frame, then ~27 ms of microtasks |
| 200 copies merge | same | 6 + 5 | 8.3 \| 31 | | 8.5 \| 32 | edit < 16 | first frame passes; whole task over |
| search "PXXP": rest fires → drawn | p53/Search | 10 | ~37 | | 38.9 (19 min) | < 50 | pass; typing → drawn ~238 ms with the box's 200 ms rest |
| search Next / Previous | same | 12 each | 16.0 / 15.2 | 20.1 / 17.2 | 20.1 / 17.2 | never > 16 | **over**: draw ~3.5 ms, then 6–12 ms of microtasks |
| fold / open a group: click → drawn | Track groups: Opening | 12 | ~7.5 \| ~20 | | 9.4 \| 22.2 | edit < 16 | first frame passes; ~12 ms microtask tail, not the SVG twin (unchanged) |
| pan, groups folded | same | 300 | 0.1 | 0.3 | 1.1 | | pass |
| Manhattan: hover / pan | Genomes/Views: Manhattan | 300 each | 0.0 / 0.1 | 0.1 / 0.3 | 0.5 / 0.7 | | pass |
| Manhattan: settled frame | same | 10 | 9.3 | 10.3 | 10.3 | | over 8: canvas 2.5–4 ms, the SVG twin 6–10 ms |
| zoomed side: pan / zoom, tight loop | Across sections: One Side Zoomed In | 300 each | 0.4 / 0.4 | 0.5–0.6 / 1.1 | **33–43** / 35.8 | never > 16 | **one draw-frame spike per 300** |
| same, yielding a macrotask between frames | same | 300 each | 1.0 / 1.3 | 1.5 / 3.0 | 2.3 / 5.1 | | pass |
| bridge hover (yielding) | same | 300 | 0.4 | 0.8 | **19.6** | | **one frame over** |
| layers menu: open → in DOM | Legends: Legend Placements | 12 | work 2.6, wall ~16 | | work 3.3, wall 24 | | pass; the menu shows a macrotask later |
| `height: 'content'` pan | Tooltips: By Data | 300 | 0.5 | 0.7 | **31.8–38.3** | | spikes; yielding: 0.9 / 1.3 / 1.7 |

**Real vsync** (real rAF, pane in view, At The Fit): no dropped frames (none over 33.4 ms).

| Run | Median | p95 | Max | Over 17.2 ms |
|---|---|---|---|---|
| Idle baseline | 16.7 | 18.1 | 18.7 | 39 |
| Hover stems | 16.7 | 18.1 | 18.7 | 30 |
| Hover ribbons | 16.7 | 18.0 | 18.7 | 39 |
| Wheel pan | 16.7 | 18.1 | 18.7 | 30 |

**Node tiers**, one at a time; `npx tsc -b` first (up to date):

| Tier | Load before | Result |
|---|---|---|
| Heavy `packages/seq/` | 1.92 | **187 of 187 passed**, 66 s (fifth: 174/176; the two contention failures did not recur) |
| Core perf | 2.76 | 262 of 262 passed |
| Storybook perf | 2.64 | 33 of 33 passed. Of note: the 4,000² drag frame allocates 9,004 B per moved frame, 28 GCs accounted |
| Size | 3.81 (size is not timing-bound) | every row passes |

**Sizes** (minified kB), now / bound / fifth slice → growth:

| Entry | Now | Bound | Fifth slice | Growth |
|---|---|---|---|---|
| Reading a spec | 265.0 | 268 | 260.2 | +4.8 |
| Headless rendering | 995.2 | 1,000 | 974.3 | +20.9 |
| Default viewer | 1,590.7 | 1,592 | 1,539.7 | **+51.0** |
| Minimal viewer | 1,120.9 | 1,124 | 1,086.9 | **+34.0** |
| Ring (rondo) | 45.0 | 48 | 43.9 | +1.1 |
| Bridge chunk (lazy) | 26.0 | 28 | 25.4 | +0.6 |
| Runtimes | 9.3 to 23.4 | 30 each | 9.3 to 23.1 | |

The bounds have moved since the fifth slice (265, 980, 1,545, 1,090, 45), and every bundle now sits within 1.3 to 4.8 kB of its bound.

**Over target**
1. Split patch, 200 copies: first frame 19.6 ms median (worse than 16.2), and the main thread is busy for ~49 ms. Merge's first frame passes; the whole task is ~31 ms.
2. Square release on the ring: ~36 ms to drawn (24–40 ms).
3. Search Next/Previous: ~16 ms median, 20 ms max, mostly microtask work after a ~3.5 ms draw.
4. Ring stem hover with microtasks: median 1.8, p95 2.1 against the 2 ms bound.
5. Settled frames (150 ms after a move): 9–14 ms on Signal and Manhattan. The canvas part is 2.5–4 ms; the rest is the story's SVG twin following at rest, so it is story composition.
6. Single 30–43 ms spikes in tight synthetic loops (One Side Zoomed In, By Data), roughly every 60 frames. They vanish when frames yield a macrotask, which points to deferred GC or heap work, not proven. One 19.6 ms bridge-hover frame on One Side Zoomed In.
7. Not reproduced: Signal's 25.7 ms zoom frame.

**Gaps**
1. The 20–30 ms microtask tails (split/merge, fold/open, Next) are not attributed. No profiler was available; the SVG twin is ruled out for the fold.
2. Opening the layers menu needs a macrotask, so its wall time includes waiting, not work. I report work and wall separately.
3. The GC cause of the spikes is inferred; `performance.memory` was too coarse to confirm.
4. The GRCh38 deep pan ran at 44 kb, not 28 kb.
5. Raster and compositor time show only in the real-vsync run, done on the ring alone.
6. The fifth slice's drag and release numbers were likely taken without pointer capture working, so they are not comparable.
7. No whole genome read by region, as before.

The raw logs stayed in the orchestrator's session scratch (not kept).

## The orchestrator's reading

- **The method moved**: a pointer-capture stub (the fifth slice's drags likely measured no drag) and a microtask drain
  after each flush (work the fifth's method missed). Where the cells differ they read old | with microtasks; the fifth
  slice's numbers are comparable only in the old column.
- **Over target, real**: the 200-copy split patch (first frame 19.6 ms, the task ~49 ms); the square's release on the
  ring (~36 ms to drawn); search's Next and Previous (~16 ms median, 20 ms max, a 3.5 ms draw and a microtask tail);
  the ring's stem hover at its 2 ms bound with microtasks. The microtask tails (split, merge, fold, Next) are not
  attributed.
- **Over, the stories'**: settled frames 9 to 14 ms on Signal and Manhattan, most of it the story's SVG twin following
  at rest (story composition, not the viewer).
- **Inferred, not proven**: single 30 to 43 ms frames in tight synthetic loops (One Side Zoomed In, By Data), gone when
  frames yield; read as deferred GC.
- **Sizes**: every bundle within 1.3 to 4.8 kB of its bound; the default viewer +51.0 kB and the minimal +34.0 kB over
  the close-out.
- **Decided by the author** (2026-11-24: "Yes, one perf lane"): the four real overs to FX36 before the last gate;
  the stories' settled frames and the inferred spikes to the TODO.

## After FX36 (re-measured 2026-11-24, `packages-seq` at `b27b0a4e3`, the pane displayed, load 1.85 to 2.29)

The same method; cells read old method | with microtasks, in ms.

| Case | Frames | Median | p95 | Max | Target | Before FX36 (median, max) | Verdict |
|---|---|---|---|---|---|---|---|
| 200 copies split, click to drawn | 5 | 21.1 \| 28.7 | 22.2 \| 32.5 | 22.2 \| 32.5 | edit < 16 | 19.6 \| 49, 21.3 \| 50 | **first frame still over**; the task ~20 ms shorter, still over |
| 200 copies merge | 5 | 5.6 \| 17.5 | 8.7 \| 22.5 | 8.7 \| 22.5 | edit < 16 | 8.3 \| 31, 8.5 \| 32 | first frame passes; the task **just over** |
| ring: square release to drawn (16, each changing the chromosome) | 16 | 0.2 \| 10.7 | 0.4 \| 12.2 | 0.4 \| 13.0 | never > 16 | ~36, 40 | pass |
| ring: drag moves | 96 | 0.0 \| 0.1 | 0.4 \| 3.2 | 0.6 \| 5.1 | < 8 | 0.8, 4.0 | pass |
| p53 Search: Next / Previous | 12 each | 2.5 / 2.3 \| 3.8 / 3.9 | 3.4 / 2.8 \| 5.2 / 5.1 | 3.7 / 3.4 \| 6.5 / 5.4 | never > 16 | 16.0 / 15.2, 20.1 / 17.2 | pass |
| ring: hover stems / ribbons (the item changing each frame) | 300 each | 0.1 / 0.0 \| 0.1 / 0.0 | 0.1 / 0.1 \| 0.2 / 0.1 | 0.3 / 0.3 \| 0.5 / 0.3 | < 2 | 0.0 \| 1.8, p95 0.2 \| 2.1 | pass |

**Real vsync** (real rAF, pane in view, 16 releases): idle median 16.7, max 17.4, no interval over 33.4; the drag and
release run median 16.7, max 34.0, one dropped frame in each of two releases (the frame after it), 14 releases clean.

**Not measured**: search's typing to drawn (the box's rest timer does 0.3 ms; the run lands in a later task not
attributed); the split's and merge's remaining microtask tails (~7.6 ms, ~12 ms) not profiled.

**Still over**: the split patch's first frame (21.1 ms) and its task (28.7 ms), the merge's task (17.5 ms median, 22.5
max). FX36 found the first frame's JavaScript 6 to 8 ms headless, the rest inferred canvas drawing of 200 copies; a cut
means sharing lays between identical copies (a design change).

**Carried** (the author, 2026-11-24: "Carry it"): the split patch's first frame and task and the merge's task, to the
next plan as a design question (sharing lays between identical copies; the remaining microtask tails profiled first).
