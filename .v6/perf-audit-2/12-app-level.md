# Slice 12: app level (v2)

Cross-cutting: `plugin/core` (PluginContext, managers, behaviors, render-loop policy, playback), `plugin/headless`,
`core/src/task` and `util/async-queue`, the CPU copies kept by `graphics/src/repr`, the auto-attached props, the Prism
loading sources, `apps/` and `scripts/esbuild` (bundles), `vitest.config.ts` and `tests/src/pipeline` (the regression
guard). v1 had no slice of this name: v1 → now covers the app-level items of 07, 10, 13 and 14 and the task-code bugs
03/B7, 03/B8 and 06/B9. Measured at `0944823d4`, Node 25.6, Apple M4, headless (CPU stub graphics), load 2.4–4.7 (times
±10–20 %; MB and bytes reproduce). `graphics/src/gl` and Canvas3D are frozen (WebGPU); MVS loading/animation is in
`.v6/mvs-optimizations/research/` and is not repeated.

Since v1: #2008 and #2013 removed the StructureInfo and measurement walks (one measurement walk per change and
`setMaxIndices` on components remain); #2022 cached the foreign-resource check (partly: the first update after a build
still walks); #2005/#2012 fixed 13/P1; #2014 fixed 14/N7; #2015 added the test lanes; #2013 put the pipeline bench in
the repo (CI runs none of it); #2019 added the AssetManager byte cache, through which text CIF above 512 MB decodes to a
`ChunkedBigString` (13/B1 dropped).

## Hot paths now

| Path | When | Cost now | Where it goes |
| --- | --- | --- | --- |
| `repr:auto` over many structures | N = 500, stage mode | 13.0 → 15.8 ms per structure (N = 100 → 500, profiled) | `setMaxIndices` per created component, 0.17 → 1.24 ms per structure, 617 ms = 8 % at N = 500, O(N²) (10/P4); the rest GC and heap |
| single action (theme, remove one) | N = 500, 6,077 cells | theme:one 39 ms (64 profiled), remove:one 32 ms (53) | hierarchy walk 29 ms + measurement walk 12 ms of the 64 profiled (14/N2, 10/P5); remove also `setMaxIndices` 4 ms |
| hover / select | N = 500, 1,333 representations | hover ~5 ms (8 profiled), select 4.9 ms | `repr.mark` on every representation, half of hover is the `EveryLoci` clear broadcast (14/N1) |
| merged render objects | unshared, > 100 units | 2.0× the buffers on the CPU (flat-32 gaussian 72.8 + 73.2 MB) | merged copy plus members (07/P7) |
| per-visual JS heap | every visual | 70–100 KB, 32–50 KB of it sphere extrema | 98 Vec3 arrays per sphere, ~30 unique (P2) |
| startup | every realm | module eval 126–136 ms, 54 MB heap | ElementSymbol table 50–75 ms, 21.5 MB (P1) |
| Prism start closure | page load | 86 of 132 files, 6.31 MB raw / 1.96 MB gzip | Seq\* 1.1 MB pulled in statically by the AFDB PAE reader (P9) |
| main-thread blocks | interactions, bonds, valence | ≥2.3 s in one call at 1M atoms | `findInterUnitContacts`, bond getters, `calcValenceModel` have no checkpoint (P5) |
| idle | classic apps | rAF + Canvas3D tick at 60 Hz for the tab's life | `PauseLoopWhenIdle` false outside Prism (P4) |

Already good: many-structures F1–F8 (#2013); playback per frame (#2022); marks scale with the loci (#1995); one-chain
layer smoothing (#2005, #2012); Prism's event consumers coalesce per frame (`watched.ts`) or are O(1) per event.

## v1 → now

20 v1 items: **fixed 3, partly fixed 5, still open 10, regressed 0, dropped 2** at the audit; since, #2026 fixed 10/P4 and 10/P5 (fixed 5, partly fixed 3).

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| 07/P7 | Merge policy by count; merged objects keep a second CPU copy | still open | | `resolveRenderMerged` `n > 100` ([base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127)); flat-32 gaussian 145.9 MB for 72.8 MB of geometry |
| 10/P4 | StructureInfo work per created object, O(N²) | **fixed** | #2008, #2013, #2026 | walks gone; `setMaxIndices` still runs for every component Structure ([structure-info.ts:223](../../packages/plugin/core/src/behavior/dynamic/custom-props/structure-info.ts#L223)): 8 % of `repr:auto` at N = 500 |
| 10/P5 | Measurement sync: six walks per change | **fixed** | #2008, #2013, #2026 | one `findTagInSubtree` miss per change (cache key `tree.transforms` is new on every edit); comment at `measurement.ts:448` wrong |
| 13/P1 | Instanced-surface smoothing spans every instance | fixed | #2005, #2012 | asm-864-low one chain 12.6 s → 31 ms (slices 08, 09) |
| 13/P2 | quality `auto` by total atoms | still open | | [util.ts:148](../../packages/graphics/src/repr/util.ts#L148) scores `data.root`; blocker 13/P1 is fixed |
| 13/P3 | Unshared 5M: per-element allocation | partly fixed | #2007, #2017 | flat-864 model 2.7 → 2.1 GB, cartoon 6.0 → 5.5 GB, bonds-intra 1.15 → 0.84 GB allocated |
| 13/P4 | Recolor pays a per-visual foreign-resource walk | **partly fixed** (v1 audit said fixed) | #2022 | `createForeignResourcesCheck` remembers only a checked clean pair; the first update after a build still walks (39 % of the flat-192 first cartoon recolor, slice 07) |
| 13/P5 | Overpaint/transparency touch every visual | fixed | #1995 | known-clear early return (slice 07) |
| 13/P6 | ASA ignores symmetry | still open | | asm-864 26.8 s vs flat-864 26.2 s; on demand only (spec question) |
| 13/B1 | Text mmCIF stops at ~4.8M atoms | dropped | | app paths decode bytes with `utf8ReadLong` → `ChunkedBigString` ([utf8.ts:112](../../packages/core/src/util/utf8.ts#L112)); cost now speed, ~1.7–2× slower parse (H5) |
| 14/N1 | Hover/select visit every representation | still open | | [representation.ts:54](../../packages/plugin/core/src/behavior/dynamic/representation.ts#L54), `:132` mark without a repr |
| 14/N2 | Hierarchy rebuilt by a whole-state walk per change | still open | | [hierarchy.ts:253](../../packages/plugin/core/src/state/manager/structure/hierarchy.ts#L253); 29 of 64 ms profiled theme:one |
| 14/N3 | Scripts without a transaction are O(N²) | still open | | none mode N = 200 `repr:auto` 7.3 s / 11.4 GB; no docs line |
| 14/N4 | Cartoon meshes reserve ~1.4× vertices | still open | | estimate unchanged ([polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86)) |
| 14/N5 | Per-structure linear work (bonds, parse, trace) | dropped | | per atom, owned by slices 05 and 08 |
| 14/N6 | 13/P4 re-measured at many structures | **partly fixed** | #2022 | same code as 13/P4 |
| 14/N7 | Builder removal writes into `state.tree` (HashCollisionNode) | fixed | #2014 | `renewOwner` in `transient.ts` |
| 06/B9 | `CustomStructureProperty.attach` without in-flight dedup | still open | | [custom-structure-property.ts:70](../../packages/model/src/props/common/custom-structure-property.ts#L70) |
| 03/B7 | Late abort discards a completed result undisposed | still open | | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) |
| 03/B8 | `chunkedSubtask` sizing compounds | still open | | [chunked.ts:28](../../packages/core/src/task/util/chunked.ts#L28); no commit touched `core/src/task` |

## Hypotheses, ranked

| Rank | Id | Hypothesis | Today | What it takes | Expected win | Effort | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | H7 | Startup and bundle size per app | 126–136 ms eval; Prism start closure 6.31 MB raw / 1.96 MB gzip | P1, P9, B2; webgl variant of the Viewer | −50–75 ms and −21.5 MB per realm; −1.1 MB raw from Prism's start | S | confirmed; **P1, P9, B2 fixed** (#2026): −22 MB heap, eval ~2× faster; start 6.31 → 5.16 MB raw |
| 2 | H9 | Event storms and per-action whole-state walks | 3 `stateUpdated` per cell; classic tree O(nodes) per event; hierarchy + measurement walks per change | 10/P4, P3, 10/P5 (S); 14/N2 (M) | `repr:auto` −8 % at N = 500; single actions 39 → ~15 ms; −0.5 s per bulk action with a 5k-node tree | S–M | confirmed; **10/P4, 10/P5 fixed** (#2026): `repr:auto` −6 %, single actions −20–35 %; 14/N2, P3 open |
| 3 | H6 | CPU copies after GPU upload | merged copy 2×; 32–50 KB extrema per visual | P2 (S), merge policy 07/P7 (S, needs a GPU check) | flat-192 ball-and-stick retained 2.3 → ~1.3 GB; −80–125 MB per representation at flat-192 | S–M | **P2 refuted** (#2026: kept extrema are distinct); 07/P7 open (GPU check) |
| 4 | H8 | Render loop policy | classic apps tick at 60 Hz forever; no adaptive resolution | P4 (S); adaptive resolution is Canvas3D (frozen) | idle wakeups 60/s → 0 after 250 ms | S | confirmed (not measurable headless) |
| 5 | H10 | Perf regression guard | bench in repo, run by nothing; no baselines; errors read as speedups | B1, B4, counters in quick/medium, nightly bench vs baseline | catches R1-like (×2.1 alloc) and O(N²) regressions | M | plausible |
| 6 | H1 | Worker pool for surfaces, ASA, then bonds/interactions | no Worker anywhere | ComputePool + typed-array protocol; stage 1 surfaces/ASA | flat-192 molsurf 19.9 → ~3–4 s wall; responsive main thread | L | plausible |
| 7 | H3 | Latest-wins and abort in the state/task layer | FIFO queue, no coalescing, a queued update cannot be dropped before it starts | P5 checkpoints (M), opt-in coalescing (M) | K × full cost → 1× for API bursts; abort latency s → ≤100 ms | M | plausible |
| 8 | H2 | Progressive/LOD loading, quality per symmetry group | quality by root atoms; reprs show after all visuals | 13/P2 (S); coarse-first refine (M–L, needs H3) | better look on assemblies; first pixel ~5–10× sooner (estimate) | M | plausible |
| 9 | H4 | Persistent cache of computed props and downloads | in-memory ref-counted byte cache only | opt-in content-keyed PropertyCache (Prism first) | repeated loads only: molsurf 6.5 → ~0.2 s | L | plausible (defer) |
| 10 | H5 | Byte-level text CIF tokenizer (**done** #2036: flat-192-cif parse 644 → 310 ms median) | BCIF default; `ChunkedBigString` above 512 MB | tokenizer over Uint8Array | above 512 MB parse ~1.7–2× faster | L | confirmed (low; defer) |

### H7 Startup and bundle size

- **Exists today.** `import` of PluginContext + DefaultPluginSpec: 126–136 ms eval, 54 MB heap; `types.ts:17` 50.2 ms
  self (startup profile). Classic Viewer: one IIFE, 5.00 MB raw / 1.44 MB gzip, WebGL 367 kB + WebGPU 306 kB + mediabunny
  231 kB eager. The Viewer has single-backend builds already (`molstar.variants`, `--variants`), used only by
  `viewer-smoke`. Prism: 86 of 132 chunks static, 6.31 MB raw / 1.96 MB gzip; its static closure holds WebGPU
  (`compose-plugin` → `graphics-backends` → `backends/all`), the query-language transpilers (131 kB, via
  `default-spec` → `targets/languages/all`) and all of Seq\* (1.11 MB, P9).
- **Takes.** P1 (S); P9 (S); B2 (S) so the start closure has a budget; deploy the Viewer's `webgl` variant as the
  default page (S, no code), or an ESM Viewer build (M).
- **Win.** −50–75 ms and −21.5 MB per realm (every app, vitest worker, future worker); Prism start −1.1 MB raw (~−0.3 MB
  gzip, estimate from the closure's ratio); Viewer webgl variant −~0.3 MB raw.
- **Risks.** None for P1/P9; a lazy WebGPU backend in Prism makes backend selection async (descriptor modules in `gl`,
  frozen).
- **First step.** P1 and P9 in one PR, with B2's closure budget so the next static Seq\* import fails the build.
- **Done (#2026).** P1: lazy memo (≤ 4096 symbols). Import of PluginContext + DefaultPluginSpec (bundled, Node 25, 8
  interleaved pairs at load 6–10): heap after GC 35.3 → 13.3 MB, eval 152–230 → 79–119 ms (median 200 → 93; 128–140 ms
  before at load ~3). P9: `readPae` in `sequence/molstar/pae-matrix` (no imports), imported narrowly: start closure
  6.31 → 5.16 MB raw, 1.96 → 1.55 MB gzip. B2: the entry's budget is on its start closure (5.5 MB); the old 3.2 MB
  budget fails the build at 5.16 MB, the barrel import back (6.31 MB) fails 5.5 MB. Viewer `webgl` variant not done.

### H9 Event storms and whole-state walks

- **Exists today.** Each evaluated cell emits `pending`, `processing`, `ok` ([state.ts:852](../../packages/core/src/state/state.ts#L852)),
  plus an object event and a log entry. Classic UI: every mounted node (collapsed children stay mounted,
  `display: none`) holds 4 global subscriptions (P3). Plugin managers: hierarchy walk per non-transaction change
  (14/N2), measurement walk per change (10/P5), `setMaxIndices` per created component (10/P4). Prism: no storm (slice 10:
  2–3 % of a load at N = 300).
- **Takes.** 10/P4 root-only check (S); P3 dispatcher (S); 10/P5 group ref from cell events (S); 14/N2 incremental
  hierarchy with a differential test (M).
- **Win.** `repr:auto` −617 ms (8 %) at N = 500, ~16 % at N = 1000; theme:one 64 → ~23 ms profiled without both
  walks; bulk actions −0.5 s with a 5k-node classic tree.
- **Risks.** 14/N2 must match the walk's ref lists exactly (Prism and the classic UI consume them).
- **First step.** 10/P4 with a test counting root visits per created component (0).
- **Done (#2026): 10/P4, 10/P5.** many N = 500, stage mode, 3 interleaved pairs at load 5–6 (medians): `repr:auto` 8.89
  → 8.32 s (−6 %; profile: `setMaxIndices` 796 ms = 8.6 % → 0); theme:one 58 → 46 ms, overpaint:one 57 → 37, remove:one
  52 → 37; allocated per single action 21 → 16–17 MB. The measurement walk is gone (theme:one 5.2 → 0, remove:one 5.9
  → 0 ms profiled); the hierarchy walk (14/N2) is now 12–16 ms, 50–66 % of a profiled single action: the rest of the
  "~15 ms". 14/N2 and P3 open.

### H6 CPU copies after GPU upload

- **Exists today.** Render objects keep their ValueCell arrays (updates, context restore). Merged objects (`'auto'` above
  100 members) copy every member buffer while the members keep theirs ([merged.ts:59](../../packages/graphics/src/repr/merged.ts#L59),
  copy in `gl/renderable/merged.ts`, frozen). Per visual 70–100 KB of JS heap, 32–50 KB of it sphere extrema.
- **Takes.** P2 (S, bit-identical, measured); 07/P7 policy: `'auto'` merges only small members (S), after a browser
  draw-call measurement; releasing member arrays is the frozen alternative.
- **Win.** flat-32 probe (CPU stub, heap-snapshot diff):

  | Representation | Merged RO | Members | Array buffers | JS heap | Sphere/extrema |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | gaussian | 72.8 MB | 73.2 MB | 145.9 MB | 28 MB | 20.8 MB |
  | cartoon | 35.1 MB | 32.2 MB | 67.3 MB | 36 MB | 9.9 MB |
  | spacefill | 4.5 MB | 7.5 MB | 12.0 MB | 29 MB | 13.3 MB |

  flat-192 ball-and-stick 991 MB buffers / 2,335 MB retained → ~1.3 GB without the merged copy; flat-864 10.6 → ~6 GB.
- **Risks.** Hundreds of unmerged draw calls (GPU frame time unmeasured).
- **First step.** P2.
- **P2 refuted (#2026, measured).** The extrema spheres keep are distinct points: unit spheres include atomic radii
  (98 of 98 distinct on flat-32's 416 units), visual spheres are padded (`Sphere3D.expand`) or over mesh vertices
  (flat-32 spacefill: 834 spheres, 81,732 extrema, 81,720 distinct values). The probe's ~30 of 98 were spheres of atom
  centers without radii. Sharing equal points by reference (bit-identical, kept by copies and transforms) was built and
  measured: flat-192 retained spacefill 245 → 247, cartoon 337 → 338, ball-and-stick 2,345 → 2,348, gaussian 181 → 183
  MB (no change); not shipped. What would cut it: packed `Float64Array` extrema (the Sphere3D type, read in frozen
  `gl/renderable/util.ts`) or fewer extrema on visual spheres (14 instead of 98: other bounds, a maintainers' call).
- **07/P7 open.** The policy trades CPU memory for draw calls: it needs a browser frame-time check on a GPU, which this
  headless setup cannot run; `resolveRenderMerged(value, unitCount)` sees no member sizes, so the size input is a
  change at its callers; releasing member arrays after merging is `gl/renderable/merged.ts` (frozen).

### H8 Render loop policy

- **Exists today.** `PluginAnimationLoop` sleeps after 250 ms idle only with `PauseLoopWhenIdle` (default false,
  [config.ts:69](../../packages/plugin/core/src/config.ts#L69)); only Prism turns it on (`pauseWhenIdle` default true).
  `pixelScale` fixed; `isActivelyInteracting` gates only illumination.
- **Takes.** P4 (S). Adaptive resolution while dragging belongs inside Canvas3D (frozen): hand to the WebGPU owners.
- **Win.** Idle CPU and battery of every classic tab: 60 rAF ticks/s → 0. Not measurable headless.
- **Risks.** Embedders mutating the scene past the API without `requestDraw` would not redraw (the apps are first-party).
- **First step.** Set it in `apps/viewer` behind a Viewer option, default on.

### H10 Perf regression guard

- **Exists today.** Lanes by tag ([vitest.config.ts:23](../../vitest.config.ts#L23)); heavy lane nightly on
  ubuntu-latest, Node 24; `itWithin` only in Seq\*; Mol\* core has two heavy ratio tests. `tests/src/pipeline` (bench,
  many, playback) is run by nothing in CI; no baseline JSON; B1 crashes on coarse entries; B4 records errors as fast
  stages; the Alloc column needs Node 25 (`total_allocated_bytes`).
- **Takes.** Quick/medium, clock-free: render-buffer bytes per representation on 7qpd and 1tup-2x2x2 (exact), events per
  update, StructureInfo root visits per component (0), marks per hover, a no-threshold smoke of `bench.ts` on 7qpd and
  8zzh. Nightly: `pipeline-bench run` on 7qpd and flat-32 `--runs 5`, `many --count 100/300`, `playback`, against a
  committed baseline (time +30 % median with one rerun, alloc +10 %, buffers exact, N = 300/100 per-structure ratio ≤ 1.25,
  fail on any error stage).
- **Win.** Catches R1 (cartoon alloc ×2.1) and the 10/P4 O(N²) automatically.
- **Risks.** Shared-runner noise: time thresholds stay coarse; baselines move deliberately with each perf PR.
- **First step.** B1 + B4, then the clock-free counters.

### H1 Worker pool

- **Exists today.** No `Worker` in `packages/`, `apps/`, `extensions/` sources. Cooperative yields at the task manager's
  100 ms rate. Surface kernels already take typed arrays (`PositionData`, `calcMolecularSurface`
  [molecular-surface.ts:67](../../packages/core/src/math/geometry/molecular-surface.ts#L67), `GaussianDensityCPU`);
  Model/Structure do not transfer. No SharedArrayBuffer on GitHub Pages (no COOP/COEP): inputs copied, 16 B per atom.
- **Takes.** A `ComputePool` in core (typed-array protocol, RuntimeContext proxy for progress and abort, fallback);
  a worker entry per app (Prism ESM `new Worker(new URL(...))`, a second IIFE bundle for the classic Viewer,
  `worker_threads` headless). Stage 1: molecular/gaussian surface + marching cubes per unit, then ASA by atom chunks.
  Stage 2: bonds (CSR per unit) and interactions per unit pair.
- **Win.** Main-thread seconds now: molsurf flat-32 6.5 s, flat-192 19.9 s; ASA 5.5 s at 1M; interactions 4.6 s at
  flat-192. Per-unit parallel on 8 workers: flat-192 molsurf → ~3–4 s wall (estimate); entries with few units bounded by
  the largest unit.
- **Risks.** Two code paths, worker bootstrap in embedders, memory doubling during transfer; each worker pays module
  eval (P1 first).
- **First step.** H10's guard, P1, then a pool prototype on `calcMolecularSurface` with a ratio test.

### H3 Scheduling and cancellation

- **Exists today.** `updateTree` enqueues on a FIFO `AsyncQueue` ([state.ts:376](../../packages/core/src/state/state.ts#L376));
  every queued update runs in full. The `onAbort` dequeue (`state.ts:418`) runs only after the update's body threw at a
  checkpoint, when `handled` already removed it, and the plugin's task manager forwards an abort only through a task's
  progress callback, which a queued update never calls: a queued update cannot be dropped before it starts. App-level
  latest-wins exists case by case (snapshot `SingleTaskQueue`, animation `applyAsync`, Prism scrubber, classic forms
  disabled while busy); API, agent and MVS callers get none.
- **Takes.** P5 checkpoints (M); opt-in `updateTree(tree, { coalesce: key })` that resolves a dropped update with the
  newer result (M); P6's O(1) waiter signalling (S).
- **Win.** Script/agent bursts of same-cell updates K × → 1× (a superseded flat-32 molsurf costs 6.5 s); abort latency in
  interactions/bonds/valence seconds → ≤100 ms. Nothing for a single UI user (already guarded).
- **Risks.** Update semantics (revertIfAborted, history for dropped updates).
- **First step.** P5's checkpoints in `findInterUnitContacts` and `calcValenceModel`.

### H2 Progressive loading, quality per group

- **Exists today.** 13/P2: quality by root `elementCount` (asm-864 gets `lowest`). `auto-lod` builds all three
  representations and toggles visibility. A representation reaches the canvas after all its visuals.
- **Takes.** 13/P2 (S); presets applying `lowest` first and refining on idle (M–L, needs H3 to cancel a pending refine).
- **Win.** asm-864 at `low`: molsurf 299 vs 166 ms at `auto` (look, +0.1–0.2 s); first pixel ~5–10× sooner for surfaces
  (estimate).
- **Risks.** Every large assembly looks different; a refine doubles the work when the coarse pass is short-lived.
- **First step.** Maintainers' answer on 13/P2.

### H4 Persistent cache

- **Exists today.** `DefaultAssetManager` URL byte cache, ref-counted, in memory ([assets.ts:160](../../packages/core/src/util/assets.ts#L160));
  IndexedDB only in Prism's journal. `Structure.hashCode` is not a content key (no coordinates).
- **Takes.** Opt-in `PropertyCache` keyed by SHA-256 of the source bytes + model index + derivation + params + code
  version; typed arrays only.
- **Win.** Repeated loads only (molsurf flat-32 6.5 s → ~0.2 s, ASA at 186k 0.9 s → ~10 ms); first loads nothing.
- **Risks.** Invalidation, quotas, stale results when params miss the key. **Defer.**

### H5 Text CIF tokenizer

- **Exists today.** PDB downloads default to BCIF; `parseCifText` takes a StringLike ([parser.ts:814](../../packages/io/src/reader/cif/text/parser.ts#L814)).
- **Measured** (flat-192.cif, 124 MB, 3 interleaved rounds under load): string 688–1,355 ms, `ChunkedBigString` (2^24)
  1,319–4,073 ms; ratio 1.7–2×; UTF-8 decode 5–10 ms.
- **Win.** Only above 512 MB of text: parse ~2× faster, one fewer copy (~0.5 GB at 5M atoms). **Defer**; document BCIF
  for >1M atoms.
- **Done (#2036):** text CIF is tokenized from its bytes (`TextByteTokenizer` for the other text formats, #2040):
  flat-192-cif parse 594–726 → 285–327 ms, retained 320 → 256 MB.

## Open findings, ranked

| Rank | Id | Status vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 10/P4 | **fixed** (#2026) | `setMaxIndices` on every component Structure (O(N²)) | [structure-info.ts:223](../../packages/plugin/core/src/behavior/dynamic/custom-props/structure-info.ts#L223) | time | medium | `repr:auto` −617 ms (8 %) at N = 500, ~16 % at N = 1000 | S | confirmed | no |
| 2 | P1 | **fixed** (#2026) | ElementSymbol table of 242k strings at module load | [types.ts:17](../../packages/model/src/model/structure/model/types.ts#L17) | time+memory | medium | −50–75 ms and −21.5 MB per realm | S | confirmed | no |
| 3 | P9 | **fixed** (#2026) | AFDB PAE reader pulls all of Seq\* into Prism's start closure | [afdb.ts:9](../../packages/extensions/model-archive/src/pae/afdb.ts#L9) | bundle | medium | −1.14 MB raw (~−0.3 MB gzip) from the start closure | S | confirmed | no |
| 4 | 07/P7 | still open | Merge by count only; members keep a full CPU copy | [base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127) | memory | high | unshared 1M ball-and-stick retained −1 GB, 5M −4.5 GB | S | confirmed | no (release path is) |
| 5 | 14/N1 | still open | Hover/select mark every representation | [representation.ts:54](../../packages/plugin/core/src/behavior/dynamic/representation.ts#L54) | time | medium | hover ~5 → <1 ms, select 4.9 → <0.5 ms at N = 500 (before the frozen scene update) | M | confirmed | no |
| 6 | 14/N2 | still open | Hierarchy rebuilt by a whole-state walk per change | [hierarchy.ts:253](../../packages/plugin/core/src/state/manager/structure/hierarchy.ts#L253) | time | medium | single actions at N = 500 39 → ~15 ms | M | confirmed | no |
| 7 | P2 | **refuted** (#2026) | Sphere extrema: 98 Vec3 arrays per sphere, ~30 unique | [boundary-helper.ts:139](../../packages/core/src/math/geometry/boundary-helper.ts#L139) | memory | medium | none: kept extrema are 98 distinct points (H6) | S | refuted | no |
| 8 | P4 | new | Classic apps never pause the render loop | [plugin-spec.ts:84](../../apps/viewer/src/plugin-spec.ts#L84) | time | medium | idle 60 ticks/s → 0 | S | confirmed | no |
| 9 | P3 | new | Classic state tree: every node subscribes to every state event | [tree.tsx:79](../../packages/plugin/ui/src/state/tree.tsx#L79) | time | medium | −0.5 s per bulk action with 5k mounted nodes | S | confirmed | no |
| 10 | P5 | new | Synchronous blocks without yield/abort (inter-unit contacts, bond getters, valence) | [interactions.ts:398](../../packages/model/src/props/computed/interactions/interactions.ts#L398) | time | medium | seconds-long blocks → ≤100 ms slices; abort works | M | confirmed | no |
| 11 | 10/P5 | **fixed** (#2026) | Measurement manager walks the whole tree per change | [measurement.ts:426](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L426) | time | low | −1.2–1.7 ms per change at 4–6k cells | S | confirmed | no |
| 12 | P7 | new; **fixed** #2034 | Focus recipes return new objects: mutative deep walk + `console.warn` | [structure-focus-representation.ts:212](../../packages/plugin/core/src/behavior/dynamic/selection/structure-focus-representation.ts#L212) | time | low | one warn and one walk per focus change and playback frame | S | confirmed | no |
| 13 | P8 | new; **fixed** #2048 | `getCartoonParams` clones the param tree ~6× per update | [cartoon.ts:96](../../packages/graphics/src/repr/structure/representation/cartoon.ts#L96) | time | low | −1.6 % of `repr:auto` | S | confirmed | no |
| 14 | 14/N4 | still open | Cartoon meshes reserve ~1.4× vertices | [polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86) | memory | low | −3–5 % of `repr:auto` on many small structures | S | confirmed | no |
| 15 | P6 | new | AsyncQueue: no coalescing, O(K²) waiter signalling | [async-queue.ts:18](../../packages/core/src/util/async-queue.ts#L18) | time | low | bursts K× → 1× with opt-in coalescing | M | confirmed | no |
| 16 | 14/N3 | still open; **fixed** #2048 | Scripts without a transaction pay every walk per commit | [context.ts:522](../../packages/plugin/core/src/context.ts#L522) | time | low | ×2–3 by wrapping in `dataTransaction` | S | confirmed | no |
| 17 | 13/P2 | still open | quality `auto` by total atoms | [util.ts:148](../../packages/graphics/src/repr/util.ts#L148) | architecture | low | look on assemblies, +0.1–0.2 s | S | confirmed (question) | no |
| 18 | 13/P3 | partly fixed | Unshared 5M per-element allocation | [polymer-trace-mesh.ts:62](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L62) | memory | low | −0.5–2 s GC per stage at 5M | M | plausible | no |
| 19 | 13/P6 | still open | ASA ignores symmetry | [accessible-surface-area.ts:62](../../packages/model/src/props/computed/accessible-surface-area.ts#L62) | time | low | assemblies 26 s → ~1/copies (spec) | M | confirmed (question) | no |

### 10/P4 `setMaxIndices` on every component (medium)

- **Mechanism.** `object.created` calls `setMaxIndices()` when `isModelOrStructure(o.obj)` (`:13,223`), true for every
  StructureComponent, selection and decorator. It loops all tracked roots three times (`maxModelIndex`,
  `setModelMaxIndex`, `setStructureMaxIndex`, each through `rootObjects` and the `pending` set with `isRoot` walks). For a
  component nothing can change: the max is over root models, and only root structures get it. `object.removed` (`:229`)
  does the same.
- **Evidence.** many-100 vs many-500 `repr:auto` profiles: `setMaxIndices` 0.17 → 1.24 ms per structure, equal to
  `update2`'s per-structure growth (+1.05 ms); 617 ms of 7.9 s profiled, `rootObjects` 354 ms self. remove:one at N = 500:
  4.4 ms of 53.
- **Fix.** Call `setMaxIndices()` only when the created/removed object is a root model (`rootModels.has(o.ref)` after
  `trackIfRoot`, or the removed ref was one); a new root structure gets `setMaxIndexOf(o.obj)` (O(1)); components and
  decorators nothing. Root-model creation stays O(N) (every root's MaxIndex moves), once per loaded entry in the model
  stage. Test: root visits per created component = 0.
- **Risk.** Low; the StructureInfo tests pin the index semantics.
- **Fixed (#2026)** as above, plus a recompute when an update brings a root model back or takes it out of the roots (a model recreated after its update failed; found in review); test: index, max index reads per component and representation 0, per root structure 1. Numbers under H9.

### P1 ElementSymbol table (medium)

- **Mechanism.** `_esCache` fills an object with every 1–3 character string over [A-Za-z0-9] (62 + 62² + 62³ = 242,234
  keys, each `toUpperCase()`) when `model/types.ts` is first imported: every app, vitest worker and future worker.
  `ElementSymbol` sees ~100 distinct values in practice.
- **Evidence.** Startup profile: 50.2 ms self at `types.ts:17` of 126–136 ms eval. Isolated (`escache.mjs`, 2 runs):
  57–67 ms, 21.5 MB retained after GC; lookup 5.1–5.6 ns (table) vs 4.8–4.9 ns (lazy Map).
- **Fix.** Lazy memo: `Map<string, string>`, `toUpperCase` on a miss, insert while `size < 4096` (bounds hostile
  input). Same strings out.
- **Risk.** None.
- **Fixed (#2026)**: heap −22.0 MB, eval ~2× faster (numbers under H7).

### P9 Seq\* in Prism's start closure (medium, new)

- **Mechanism.** `packages/extensions/model-archive/src/pae/afdb.ts:9` imports `readPae` from the barrel
  `@molstar/prism-viewer/sequence/molstar` (`export *` of 14 modules, `link-views` → `@molstar/seq-viewer`). The chain
  `index.ts → spec → full-spec → extensions → model-archive/prism-viewer → pae/operations → pae/afdb → sequence/molstar`
  is all static, so Seq\* loads at start although the spec says it is its own chunk, loaded when a structure is first
  shown (`spec.ts:105`, §2.3).
- **Evidence.** Prism scratch build's metafile: removing that one edge drops 273 inputs, 1,140 kB raw (1,109 kB of it
  `sequences/viewer` + `sequences/core`) from the 6.31 MB static closure. A narrow import of `pae.ts` alone would still add
  514 kB (`seq-core/format`).
- **Fix.** Load the reader lazily where a PAE JSON is read (`await spec.sequence.runtime()` or
  `await import('@molstar/prism-viewer/sequence/molstar')` in the async PAE operation); B2's closure budget keeps it out.
- **Risk.** None (the PAE card is async already).
- **Fixed (#2026)** with a narrow import instead: `PaeMatrix` and `readPae` moved to `sequence/molstar/pae-matrix.ts`
  (no imports; `pae.ts` re-exports them), exported as `@molstar/prism-viewer/sequence/molstar/pae-matrix`. Start closure
  6.31 → 5.16 MB raw, 1.96 → 1.55 MB gzip. Test: esbuild bundle of `afdb.ts` has no Seq\* input (`scripts/esbuild/_test`).

### 07/P7 Merge policy by count (high at scale)

- **Mechanism.** `resolveRenderMerged('auto', n)` is `n > 100` whatever the member size; `mergeRenderObjectsByType`
  ([merged.ts:59](../../packages/graphics/src/repr/merged.ts#L59)) builds values that copy every member buffer; members
  keep theirs for sync (`gl/renderable/merged.ts`, frozen).
- **Evidence.** flat-32 probe: array buffers 2.0× geometry (table under H6). flat-192 ball-and-stick 991 MB buffers /
  2,335 MB retained; flat-864 4.46 / 10.6 GB (v1 13).
- **Fix.** `'auto'` also requires small members (mean member vertex bytes under a budget, e.g. 256 KB), so unit-sized
  meshes of unshared structures stay unmerged. Measure draw-call cost in a browser first. Releasing member arrays after
  merging is the frozen alternative (WebGPU owners).
- **Risk.** Draw calls for hundreds of unmerged surfaces.

### 14/N1 Hover/select over every representation (medium)

- **Mechanism.** `HighlightLoci` and `SelectLoci` providers call `canvas3d.mark` without a repr (`:54,132`); Canvas3D's
  `applyMark` loops `reprRenderObjects` and, on a change, updates every render object (`scene.update(void 0, true)`,
  frozen).
- **Evidence.** N = 500 (1,333 representations), 21 hovers: 194 ms, profiled 166 ms, all of it `repr.mark` (canvas stub);
  `clearHighlights` (the `EveryLoci` remove broadcast) is 72 of the 166. Select 4.9 ms per event (0.54 at N = 100).
- **Fix.** In plugin/core, an index root structure → Representation3D cells, kept from object events; for structure
  loci push `canvas3d.mark({ loci, repr }, action)` per representation of that root (Canvas3D already buffers several
  marks and updates the scene once). Track which representations hold a highlight, and clear only those instead of
  broadcasting `EveryLoci`. Non-structure loci keep the broadcast. The scene update over all render objects stays for the
  WebGPU owners.
- **Risk.** Loci that map across roots (child structures, superpositions): keep `Loci.remap`/`areRootsEquivalent`.

### 14/N2 Hierarchy walk per change (medium)

- **Mechanism.** Every non-transaction `changed` calls `sync(true)` → `buildStructureHierarchy`
  ([hierarchy-state.ts:22](../../packages/plugin/core/src/state/manager/structure/hierarchy-state.ts#L22)), a pre-order walk
  of the data tree with immutable.js iteration.
- **Evidence.** N = 500 profiles: theme:one 64 ms, hierarchy 29 ms (45 %); remove:one 53 ms, hierarchy 28 ms. (v1 said
  57 % for the hierarchy: that is the hierarchy and measurement walks together.)
- **Fix.** Dirty refs from cell/object events during the update; update only the affected Structure/Model/Component refs
  (`createOrUpdateRefList`); fall back to the walk when a structural root changes; a differential test against the walk.
- **Risk.** Consumers depend on the exact ref lists.

### P2 Sphere extrema (medium)

- **Mechanism.** `getSphere` stores `this.extrema.slice()`: 98 points as separate `[x, y, z]` arrays, many boxed
  (`heap number` nodes in the snapshot). Each visual keeps several spheres (`boundingSphere`, `invariantBoundingSphere`,
  the geometry's, merged members'), and `Sphere3D.clone` copies them all (52.6 ms + `Vec3.clone` 47.7 ms in many-500
  `repr:auto`).
- **Evidence.** Heap-snapshot diff on flat-32 (416 visuals), nodes on paths through bounding spheres: spacefill 13.3 MB
  (more than its 12.0 MB of typed arrays), gaussian 20.8 MB, cartoon 9.9 MB. `extrema.ts`: 98 per unit sphere, 29.6
  unique (50 units of 1tup). `v/dedupe.ts`: order-preserving dedupe of all 416 unit spheres gives a **bit-identical**
  union sphere (center, radius, extrema); the union takes 1.51 vs 8.85 ms; dedupe costs 1.6 µs per sphere.
- **Fix.** In `getSphere`, drop repeated points keeping first occurrences in order (exact compare, ≤98²/2 compares).
  The auditor's caveat does not hold: a parent helper picks the same point per direction (ties resolve to the first
  occurrence either way), and `radiusSphere` takes a max. Packed `Float64Array` extrema later (changes the Sphere3D type
  used in frozen gl code).
- **Risk.** Low; snapshot tests on bounds should not move (bit-identical in the probe).
- **Refuted (#2026).** Dropping points is not bit-identical: `Box3D.fromSphere3D` (≥ 14), `Sphere3D.expand` (moments
  over all entries) and the merged spheres' quality (`loci.ts`, `util/boundary.ts`) read the count. Sharing by
  reference keeps it, but the kept spheres have no equal points (H6): no memory to win.

### P4 Classic apps never pause (medium)

- **Mechanism.** `PauseLoopWhenIdle` defaults to false; the Viewer's spec (`plugin-spec.ts:84`) sets other General options
  but not this one, nor do docking-viewer, mesoscale-explorer and mvs-stories; `frame` asks for a frame every vsync
  ([animation-loop.ts:90](../../packages/plugin/core/src/animation-loop.ts#L90)).
- **Evidence.** Code; Prism's loop sleeps after 250 ms idle (`compose-plugin.ts:125`, prism AGENTS "idle is idle").
- **Fix.** `[PluginConfig.General.PauseLoopWhenIdle, true]` in the apps (a Viewer option, default on); keep the library
  default and document it for embedders.
- **Risk.** Code moving the camera past the API without `requestDraw` stops redrawing.

### P3 Classic tree fan-out (medium)

- **Mechanism.** Each `StateTreeNode` subscribes to `cell.stateUpdated`, `cell.created`, `cell.removed` (`:79,91,97`);
  each label to `stateUpdated` (filtered, debounced) and `currentObject` (`:199,209`). Collapsed children stay mounted
  (`display: none`, `:171`). A bulk update of K cells costs 3K × nodes callbacks.
- **Evidence.** `fanout.ts` (rxjs Subject, 2 filtered subscriptions per node), 4,800 events (theme:all at N = 500):
  500 nodes 56–58 ms, 5,000 nodes 556–570 ms (2 runs).
- **Fix.** One subscription per stream at the tree root, `Map<ref, Set<listener>>` dispatch.
- **Risk.** Low (UI only).

### P5 Synchronous blocks (medium)

- **Mechanism.** `findInterUnitContacts` is one call over every unit pair (the intra loop above it yields per symmetry
  group); `unit.bonds` and `structure.interUnitBonds` are lazy getters; `calcValenceModel` loops units without
  `ctx.update` ([valence-model.ts:366](../../packages/model/src/props/computed/chemistry/valence-model.ts#L366)).
- **Evidence.** interactions flat-192 4.8 s with the inter-unit loop 50 % (flat) / 83 % (asm) (slice 06), so ≥2.3 s in one
  block at 1M; bonds-intra 0.82 s at flat-192, 3.6 s at flat-864; valence 0.9 s at 5M.
- **Fix.** Batch unit pairs with `if (ctx.runtime.shouldUpdate) await ctx.runtime.update(...)`; same per unit in valence;
  an async `ensureBonds(ctx, structure)` that the representation and property paths call before the getters.
- **Risk.** Async signatures in model internals; deterministic order.

### Low

- **10/P5** `groupRef()` caches by `tree.transforms`, new on every edit, and `findTagInSubtree` never stops (03/B10): one
  full walk per change with no group (Prism slice: 1.15 ms at 4,201 cells; theme:one at N = 500: 12 of 64 ms profiled).
  Track the group ref from cell events (tag on created/removed cells); fix the comment at `measurement.ts:448`.
  **Fixed (#2026):** the ref from `cell.created`/`cell.stateUpdated` tags, checked O(1) per sync, one search after the
  group is removed or untagged; the comment is rewritten.
- **P7** `(old) => ({ ...old, bundle })` recipes (`:212,217,271,272`; also `volume-tools/segmentor/manager.ts:251,269,324`)
  return a non-draft object: mutative 1.3.0's `handleReturnValue` walks it recursively (typed arrays of the bundle by
  key) and `console.warn`s with no production guard. Mutate the draft or pass the params object.
- **P8** `PD.clone(CartoonParams)` + a units scan per call, ~6 calls per update: `getCartoonParams` 126 ms of 7.9 s in
  many-500 `repr:auto`. Memoize per (structure, registries).
- **14/N4** Estimate per element type (sheets 8 vertices per segment); owned with slice 08.
- **P6** `AsyncQueue` signals every waiter per `handled` (O(K²)); with K small the cost is the full re-run of superseded
  updates, which needs H3's opt-in coalescing.
- **14/N3** Document `plugin.dataTransaction` for scripts; 14/N2 and 10/P4 shrink the penalty.
- **13/P2** Score instanced representations by the largest symmetry group's atoms (question for the maintainers).
- **13/P3** Owned by slices 04, 08 (`arrayView` copies, polymer-trace temporaries); here P2.
- **13/P6** Spec decision (per-group ASA); else off-thread (H1).

## Bugs

| Id | Status vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B2 | **fixed** (#2026) | low | Prism's size budget is on the 7 KB entry stub | [package.json:67](../../apps/prism-viewer/package.json#L67) | confirmed |
| 06/B9 | still open; **fixed** #2045 | low | `CustomStructureProperty.attach` has no in-flight dedup | [custom-structure-property.ts:70](../../packages/model/src/props/common/custom-structure-property.ts#L70) | plausible |
| B4 | new | low | Bench stages record errors as fast results | [stage.ts:99](../../tests/src/pipeline/stage.ts#L99) | confirmed |
| B1 | new | low | Pipeline bench crashes on coarse-only entries | [bench.ts:324](../../tests/src/pipeline/bench.ts#L324) | confirmed |
| 03/B7 | still open | low | Late abort discards a completed result undisposed | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) | confirmed |
| B3 | new; **fixed** #2044 | low | Plugin log entries grow without bound | [context.ts:692](../../packages/plugin/core/src/context.ts#L692) | confirmed |
| 03/B8 | still open; **fixed** #2045 | low | `chunkedSubtask` sizing compounds | [chunked.ts:28](../../packages/core/src/task/util/chunked.ts#L28) | confirmed |

- **B2.** `budgets: { 'index.js': 3200000 }` is checked per output file (`scripts/esbuild/app.mjs:304`). With
  `splitting: true`, `index.js` is a 7,385 B stub of imports; the start code is in 86 static chunks (6.31 MB). The budget
  can never fail; P9 (1.1 MB of Seq\* moved into the start closure) went unnoticed. The single-file budget
  (`prism-viewer.js` 8.0 MB, now 7.11 MB) does catch total growth. Fix: sum the entry's static `import-statement`
  closure from the metafile and budget that. **Fixed (#2026):** `scripts/esbuild/start-closure.mjs` sums it; an ESM
  entry's budget is on that (Prism 5.5 MB, closure now 5.16 MB).
- **06/B9.** `attach` returns early only when a value exists, before `await builder.obtain`. Two overlapping attaches
  (a selection query's `ensureCustomProperties`, `queries/structure/query.ts:78`, outside the state queue, while an
  ASA-colored representation builds) both compute: 0.9 s at 186k, 5.5 s at 1M; the last `set` wins. Fix: a pending
  promise per (data, descriptor, params); one requester's abort must not cancel it for the others.
- **B4.** `stage()` catches, sets `info = 'ERROR: …'` and records `ms` (also `bench.ts:209`): a broken stage reads as a
  speedup in any comparison. Fix: `error: true`, non-zero exit unless `--allow-errors`.
- **B1.** `selections()` takes `s.units.find(atomic polymer)!`; on `ihm-coarse-8zzh.cif` it is undefined and `chainId`
  throws "reading 'elements'" outside `stage()`, killing the run after `theme-eval`; no `--skip` avoids it. Fix: fall back
  to `s.units[0]` or skip the component/interactive stages.
- **03/B7.** `execute` throws `Task.Aborted` after `task.f` returned when an abort was requested meanwhile; the returned
  object (texture-mesh textures, an asset) is never stored or disposed. Fix (v1): return the value, let the abort land
  at the next checkpoint.
- **B3.** `log.entries = log.entries.push(e)` keeps every entry; every evaluated transform logs unless
  `doNotLogTiming`. 339 B per entry (`logmem.ts`): ~8k entries per 500 loads, 1,600 per bulk action at N = 500; tens of
  MB over long agent sessions. Readers: the classic Log (last 10) and MVS tests (deltas). Fix: cap at e.g. 1,000.
- **03/B8.** `start` resets only after an update, so `delta` is cumulative and `currentTime`/`globalTime` sum sums; at
  the plugin's 100 ms rate the gaps reach ~200 ms. Fix: time each chunk, size by throughput.

## Checked and dropped

- Abort token shared by tasks nested in one task-manager context (`v/abort-share.ts`): aborting A leaves a concurrently
  run B intact. Not a bug.
- A queued state update cannot be dropped before it starts (`state.ts:418` dequeue runs after `handled`): a design gap,
  recorded under H3, not a bug (the update aborts at its first checkpoint).
- v1 audit "13/P4 fixed", "14/N6 fixed": corrected to partly fixed (first update after a build still walks).
- H7 "the classic Viewer carries both backends": true for the default build, but `molstar.variants` single-backend
  builds exist (not deployed).
- P2 caveat "merged spheres not bit-identical after dedupe": refuted for order-preserving dedupe (probe identical).
- P2 itself (#2026): the probe's spheres had no atomic radii; the extrema units and visuals keep are distinct points,
  and dropping duplicates would change the counts other code reads (H6).
- P1 impact "high": lowered to medium (50–75 ms of a page start, 21.5 MB heap); still rank 2 for its effort.
- B2 impact "medium": lowered to low (the single-file budget catches total growth; only moves into the start closure
  escape).

## Bench scripts

All under `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/app-level/`
(bundled with `build.mjs`: esbuild, conditions `molstar-src, import, default`):

| Script | What |
| --- | --- |
| `prof-many100-repr/`, `prof-many500-repr/`, `prof-many500/` + `prof-top.mjs`, `cmpprof.mjs` | many-structures CPU profiles (repr:auto, theme:one, remove:one, hover, select); per-structure growth N = 100 → 500 |
| `startup.ts`, `cpuprof/` | startup profile (import PluginContext + DefaultPluginSpec) |
| `escache.mjs` | ElementSymbol table build, retained heap, lookup vs lazy Map (`node --expose-gc`) |
| `h6-probe.ts`, `snapdiff.mjs`, `snaps/` | heap-snapshot diff per representation on flat-32 (merged vs members, extrema paths) |
| `extrema.ts`, `v/dedupe.ts` | extrema uniqueness; order-preserving dedupe identity and union time |
| `fanout.ts` | classic tree subscription fan-out |
| `logmem.ts` | bytes per plugin log entry |
| `bigstr.ts`, `utf8.ts` | string vs ChunkedBigString parse of flat-192.cif; UTF-8 decode |
| `bundle/app.mjs`, `bundle/meta.mjs`, `v/closure.mjs` | production scratch builds of the apps; Prism static closure; inputs reachable only through the AFDB edge |
| `v/abort-share.ts` | task-manager abort with two tasks in one context |
