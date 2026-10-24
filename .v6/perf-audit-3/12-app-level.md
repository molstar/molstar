# Slice 12: app level (v3)

Scope: `packages/plugin/core` (PluginContext, managers, behaviors, render-loop policy), `packages/plugin/headless`,
`apps/`, `scripts/esbuild` (bundles, budgets), `tests/src/pipeline` (bench, many, the regression guard) and the test
lanes; plus the plugin side of the downloads/files rework (#2039, #2035) and the task/queue items v2 kept here (03/B7,
03/B8, P6). HEAD `16277d822`, v2 base `0944823d4`. Node 26.12.1, Apple M4, headless (CPU stub graphics).

**Measurement conditions.** Load 3 at the start, then 14–34 for most of the run (≈13 other auditors), swap 12–18 GB of
18 GB used: wall times are 2–10× v2's and swing ±50 % between back-to-back runs. Everything below compares
**allocation** (`total_allocated_bytes`), retained MB, bundle bytes and counters; times are given as ranges and never
compared across runs. The many-structures A/B bundles `tests/src/pipeline/many.ts` of each tree (base `0944823d4`
extracted with `git archive`) with the bench's esbuild setup. N = 500 ran once per tree (memory pressure; a second
N = 500 pair was stopped).

PRs since v2 that touched the slice: #2026 (app-level top 3: 10/P4, 10/P5, P1, P9, B2), #2034 (P7), #2035 (ajax by
fetch), #2039 (downloads and files as Response Data), #2016 (bond providers in the registry and root structure
params), #2037 (periodic), #2043 (Prism sequences/track providers), #2044 (B3), #2045 (06/B9, 03/B8), #2046 (P6
waiters), #2048 (P8, 14/N3 doc, 13/P4).

## Hot paths now

| Path | Fixture | Cost now | Where it goes |
| --- | --- | --- | --- |
| `repr:auto` over many structures | many, stage mode | alloc N = 100 1,028 → 731 MB (−29 %), N = 500 4,525 → 2,945 MB (−35 %); 0.86–16 s at N = 100 (load 3–25) | geometry building: `ChunkedArray.add3` 12 % self, `addTube`, `_compact`, curve interpolation (slices 01/08); app-level overhead small (param clones 0.2 %) |
| single action (theme, remove one) | many N = 200 | theme:one 16–20 ms, remove:one 16 ms (load ~25) | hierarchy walk `buildStructureHierarchy` 8.2 of 9.1 ms (theme:one) and 8.5 of 10.0 ms (remove:one) of the state update, profiled (14/N2) |
| hover / select | many N = 100/200/500 | 1,579 / 3,171 / 7,933 render-object updates per hover; alloc unchanged vs v2 | `repr.mark` on every representation 30.7 of 38 ms over 21 hovers at N = 200; `clearHighlights` 10.3 ms (14/N1) |
| model stage over many structures | many N = 200 | StructureInfo 36 ms of ~145 ms profiled (~25 %) | `setMaxIndices` per created root model over all roots: O(N²), the part #2026 kept (P13) |
| text CIF opened with format `auto` | flat-192.cif 124 MB, 7qpd | sniff 2.1–2.5 s (flat-192, load ~25), 8–11 ms (7qpd) | `guessCifVariant` scans the whole file for 5 needles, on bytes since #2039 1.3–3× slower than `String.includes` (P11) |
| startup | import PluginContext + DefaultPluginSpec | heap after import 34.9 → 13.0 MB, alloc 96–107 → 52–58 MB, after init 39.8 → 17.9 MB | eval now compile-bound: `compileSourceTextModule` 481 of 816 ms profiled under load (6.8 MB bundle) |
| Prism start closure | `apps/prism-viewer` production build | 5.290 MB raw / 1.591 gzip / 1.334 brotli, 90 files; budget 5.5 MB (96 %) | WebGPU 308 kB, MVS 234, extensions 222, transpilers 87, AI threads 27 = 877 kB exclusive (P10) |
| quick test lane | `packages/plugin`, `scripts`, `tests` (86 files) | 32.5 s summed file time (4 workers, load ~25) | `typecheck-full.test.mjs` 13.0 s (5 tsc spawns), process/exports spawns 2.1 s each (P12) |

## v2 → now

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| 10/P4 | `setMaxIndices` per component (O(N²)) | **fixed** | #2026 | test pins 0 index reads per component; repr:auto alloc −29/−35 %; the per-root-model loop stays (P13) |
| 10/P5 | Measurement walk per change | **fixed** | #2026 | `groupRef` O(1) from cell events ([measurement.ts:431](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L431)); no measurement frames in the N = 200 profiles |
| P1 | ElementSymbol table at module load | **fixed** | #2026 | heap after import 34.9 → 13.0 MB, 5 interleaved pairs |
| P9 | Seq\* in Prism's start closure | **fixed** | #2026 | no `packages/sequences` input in the start closure at HEAD |
| B2 | Prism budget on the entry stub | **fixed** | #2026 | `startBytes` on the closure ([start-closure.mjs:13](../../scripts/esbuild/start-closure.mjs#L13)); 5.29 of 5.5 MB used |
| P7 | Focus recipes return new objects (mutative walk + warn) | **fixed** | #2034 | recipes mutate the draft; mutative's `finalizeAssigned` walks only with >1 draft (1 here): no walk |
| P8 | `getCartoonParams` clones params ~6× per update | **fixed** | #2048 | `createCartoonParams` 7.5 ms of 4,059 ms profiled repr:auto at N = 200 (0.2 %) |
| 14/N3 | Scripts without a transaction O(N²) | **fixed** (doc) | #2048 | `dataTransaction` doc ([context.ts:531](../../packages/plugin/core/src/context.ts#L531)) |
| 13/P4 = 14/N6 | First update after a build walks foreign resources | **fixed** | #2048 | `Visual.foreignResourcesCheck` knows the created object ([visual.ts:129](../../packages/graphics/src/repr/visual.ts#L129)); not re-measured (slice 07) |
| 06/B9 | `CustomStructureProperty.attach` no in-flight dedup | **fixed** | #2045 | pending promise per structure; a failed/aborted first requester makes the waiter compute ([custom-structure-property.ts:86](../../packages/model/src/props/common/custom-structure-property.ts#L86)) |
| B3 | Log entries unbounded | **fixed** | #2044 | `LogSize` 1000 ([context.ts:91](../../packages/plugin/core/src/context.ts#L91)) |
| 03/B8 | `chunkedSubtask` sizing compounds | **fixed** | #2045 | chunk = processed since last update / 4 ([chunked.ts:18](../../packages/core/src/task/util/chunked.ts#L18)) |
| P6 | AsyncQueue: O(K²) waiter signalling, no coalescing | partly fixed | #2046 | waiters per value, head resolved only ([async-queue.ts:12](../../packages/core/src/util/async-queue.ts#L12)); coalescing (H3) open |
| 13/P3 | Unshared 5M per-element allocation | partly fixed | #2007, #2017 | owned by 04/08; not re-measured here |
| 07/P7 | Merge by count; members keep a CPU copy | still open | | [base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127) unchanged; needs a browser GPU check |
| 14/N1 | Hover/select mark every representation | still open | | [representation.ts:54](../../packages/plugin/core/src/behavior/dynamic/representation.ts#L54), `:132`; hover alloc 18 → 19 MB at N = 500, render-object updates per hover 7,933 (= v2) |
| 14/N2 | Hierarchy walk per change | still open | | [hierarchy.ts:253](../../packages/plugin/core/src/state/manager/structure/hierarchy.ts#L253); ~85–90 % of a single action's state update at N = 200 |
| P3 | Classic tree: every node subscribes to every event | still open | | [tree.tsx:79](../../packages/plugin/ui/src/state/tree.tsx#L79) |
| P4 | Classic apps never pause the loop | still open | | `PauseLoopWhenIdle` only in Prism ([compose-plugin.ts:125](../../packages/prism/viewer/src/sources/compose-plugin.ts#L125)); [plugin-spec.ts:85](../../apps/viewer/src/plugin-spec.ts#L85) |
| P5 | Synchronous blocks (inter-unit contacts, bond getters, valence) | still open | | [interactions.ts:455](../../packages/model/src/props/computed/interactions/interactions.ts#L455), [valence-model.ts:371](../../packages/model/src/props/computed/chemistry/valence-model.ts#L371): no checkpoint |
| 14/N4 | Cartoon meshes reserve ~1.4× vertices | still open | | [polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86) |
| 13/P2 | quality `auto` by total atoms | still open | | [util.ts:195](../../packages/graphics/src/repr/util.ts#L195) scores `data.root` (question) |
| 13/P6 | ASA ignores symmetry | still open | | spec question |
| B1 | Bench crashes on coarse-only entries | still open | | reproduced: `ihm-coarse-8zzh.cif` dies after `theme-eval:residue-name`, "reading 'elements'" ([bench.ts:324](../../tests/src/pipeline/bench.ts#L324)) |
| B4 | Bench stages record errors as fast results | still open | | [stage.ts:99](../../tests/src/pipeline/stage.ts#L99), [bench.ts:209](../../tests/src/pipeline/bench.ts#L209) unchanged |
| 03/B7 | Late abort discards a completed result undisposed | still open | | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) unchanged |
| P2 | Sphere extrema dedupe | dropped | #2026 | refuted in v2 (kept extrema are distinct) |

**Counts:** fixed 13, partly fixed 2, still open 12, regressed 0, dropped 1.

**New costs or bugs from fixes since v2:**

- #2039: the CIF variant sniff now searches bytes (`asciiIndexOf`: `Uint8Array.indexOf` on an anchor letter, then a
  compare) instead of `String.includes`: 1.3–3× slower per scan, and a regular mmCIF misses all 5 needles, so every
  text CIF opened with format `auto` pays 5 full scans (P11, partly a regression).
- #2034's test for P7 and the shared `focusStructure` helper poll with real timers (100 ms fixed + 25 ms steps): slow
  by waste in the quick lane (P12).
- #2026, #2044, #2045, #2046, #2048: hunks read, no new bug. #2026's StructureInfo keeps one O(N) loop per root model
  (P13, accepted in v2, now measured). #2039's `ResponseBody`/`HeldAsset` hand-off (`textNow` → `letGoOfBytes` →
  `toSource`) checked for two Downloads of one URL and for files: no double free, no stale bytes; every in-repo
  transform that reads Download output takes `SO.Data.Response` (grep over `from:`).

## Hypotheses, re-ranked

| Rank | Id (v2 rank) | Hypothesis | Today | Next step | Expected win | Effort | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | H9 (2) | Event storms and per-action walks | 10/P4, 10/P5 done; the hierarchy walk is now most of a single action; hover marks every representation | 14/N2 incremental hierarchy with a differential test; 14/N1 per-root mark | single actions at N = 200 ~10 → ~2 ms of state update; hover ~1.8 → <0.3 ms per event | M | confirmed |
| 2 | H10 (5) | Perf regression guard | partly: #2026 added clock-free counters (index reads per component, start-closure budget); bench run by nothing in CI, no baselines; B1, B4 open; classic apps have no size budget | B1 + B4 (S); alloc baselines for `many --count 100` and 7qpd nightly; budgets for the Viewer bundle | would have flagged #2039's sniff cost and the 130 kB Prism start growth | M | plausible |
| 3 | H7 (1) | Startup and bundle per app | P1, P9, B2 done; start closure at 96 % of its budget, 17 % of it exclusive to WebGPU, MVS, extensions, transpilers | P10 lazy MVS runtime and transpilers (S–M); WebGPU backend lazy (async selection, frozen gl descriptors) | −0.3 to −0.9 MB raw at start | S–M | confirmed |
| 4 | H6 (3) | CPU copies after GPU upload | 07/P7 open | browser draw-call check, then size-aware `'auto'` | flat-192 b&s retained 2.3 → ~1.3 GB | S (after check) | confirmed |
| 5 | H8 (4) | Render loop policy | unchanged (P4) | `PauseLoopWhenIdle` on in the classic apps | idle 60 rAF/s → 0 | S | confirmed (unmeasurable headless) |
| 6 | H3 (7) | Latest-wins and abort | P6 waiters fixed (#2046), P5 open, no coalescing | P5 checkpoints | abort latency s → ≤ 100 ms | M | plausible |
| 7 | H1 (6) | Worker pool | none | after H10 | surfaces off the main thread | L | plausible |
| 8 | H2 (8) | Progressive loading | 13/P2 open | maintainers' answer | first pixel sooner on assemblies | M | plausible |
| 9 | H4 (9) | Persistent cache | none | defer | repeated loads only | L | defer |
| — | H5 (10) | Byte-level text tokenizer | **done** #2036, #2040 | — | flat-192-cif parse 2× faster (v2) | — | closed |

## Open performance findings (ranked)

| Rank | Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 14/N2 | still open | Hierarchy rebuilt by a whole-state walk per change | [hierarchy.ts:253](../../packages/plugin/core/src/state/manager/structure/hierarchy.ts#L253) | time | measured: 8.2–8.5 ms of a 9–10 ms state update at N = 200 (profile) | M | no |
| 2 | P11 | new (partly regressed, #2039) | Text CIF variant sniff scans the whole file 5×, on bytes | [provider.ts:186](../../packages/plugin/core/src/state/formats/provider.ts#L186) | time | measured: flat-192.cif 2.1–2.5 s (string path 1.0–1.8 s), 7qpd 8–11 ms (2.6–3.6), load ~25; parse of the same file ~0.3 s (v2) | S–M | no |
| 3 | 14/N1 | still open | Hover/select mark every representation | [representation.ts:54](../../packages/plugin/core/src/behavior/dynamic/representation.ts#L54) | time | measured: 7,933 render-object updates per hover at N = 500; 30.7 of 38 ms (21 hovers, N = 200) in `mark` | M | scene update yes |
| 4 | 07/P7 | still open | Merge by count; members keep a CPU copy | [base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127) | memory | v2: flat-192 b&s retained −1 GB | S | release path yes |
| 5 | P10 | new | Prism start closure: 877 kB exclusive to WebGPU, MVS, extensions, transpilers; budget at 96 % | [compose-plugin.ts](../../packages/prism/viewer/src/sources/compose-plugin.ts), [graphics-backends.ts:7](../../packages/plugin/core/src/util/graphics-backends.ts#L7), [default-registry.ts](../../packages/plugin/core/src/default-registry.ts) | bundle | measured: −308 kB (WebGPU), −234 (MVS), −222 (extensions), −87 (transpilers), −27 (AI) raw, exclusive | S–M | WebGPU part yes |
| 6 | P4 | still open | Classic apps never pause the render loop | [plugin-spec.ts:85](../../apps/viewer/src/plugin-spec.ts#L85) | time | idle 60 ticks/s → 0 (estimate) | S | no |
| 7 | P12 | new | Quick lane: process-spawning tests and real-timer waits | [typecheck-full.test.mjs:9](../../scripts/workspace/_test/typecheck-full.test.mjs#L9) | test time | measured: 13.0 + 2.1 + 2.1 s of 32.5 s summed (86 files); focus tests ≥125 ms wait each | S | no |
| 8 | P3 | still open | Classic tree subscriptions fan out | [tree.tsx:79](../../packages/plugin/ui/src/state/tree.tsx#L79) | time | v2: −0.5 s per bulk action, 5k nodes | S | no |
| 9 | P5 | still open | Synchronous blocks without yield/abort | [interactions.ts:455](../../packages/model/src/props/computed/interactions/interactions.ts#L455) | time | v2: ≥2.3 s blocks at 1M atoms | M | no |
| 10 | P13 | new | StructureInfo: one loop over all roots per created root model (O(N²) model stage) | [structure-info.ts:230](../../packages/plugin/core/src/behavior/dynamic/custom-props/structure-info.ts#L230) | time | measured: 36 ms of ~145 ms model stage at N = 200; ~0.2 s at N = 500 (estimate, ×6) | M | no |
| 11 | P6 | partly | AsyncQueue: no coalescing | [async-queue.ts:9](../../packages/core/src/util/async-queue.ts#L9) | time | bursts K× → 1× (estimate) | M | no |
| 12 | 14/N4 | still open | Cartoon vertex reserve ~1.4× | [polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86) | memory | −3–5 % repr:auto (v2 estimate) | S | no |
| 13 | 13/P2 | still open | quality `auto` by total atoms | [util.ts:195](../../packages/graphics/src/repr/util.ts#L195) | architecture | look on assemblies | S | no |
| 14 | 13/P6 | still open | ASA ignores symmetry | [accessible-surface-area.ts](../../packages/model/src/props/computed/accessible-surface-area.ts) | time | spec question | M | no |

### 14/N2 Hierarchy walk per change

- **Evidence.** many N = 200, stage mode, CPU profiles (`prof200/`): theme:one `_updateTree` 9.1 ms, of it
  `buildStructureHierarchy` 8.2 ms (pre-order over immutable.js); remove:one 10.0 ms, hierarchy 8.5 ms. With 10/P5 gone
  the walk is nearly all of a single action's manager work. Allocation per single action: theme:one 21 → 16 MB,
  remove:one 22 → 17 MB at N = 500 (base → HEAD).
- **Fix** (unchanged from v2): dirty refs from cell/object events, update the affected refs, the walk as fallback when a
  structural root changes; a differential test against the walk (Prism and the classic UI read the ref lists).

### P11 Text CIF variant sniff (new; partly #2039)

- **Mechanism.** `DataFormatRegistry.auto` asks each CIF provider; every one calls `guessCifVariant`
  ([provider.ts:151](../../packages/plugin/core/src/state/formats/provider.ts#L151); `lastSniff` makes it once per
  data). For `.cif` it runs 2 `startsWith` and up to 5 `sniffText` searches over the whole data (`:186–198`): a plain
  mmCIF matches none, so all 5 scan to the end. Since #2039 a response is held as bytes and searched with `asciiIndexOf`
  ([utf8.ts:258](../../packages/core/src/util/utf8.ts#L258)): `Uint8Array.indexOf` of the first lowercase letter, then
  a byte compare per hit.
- **Evidence** (`sniff.ts`, 3 rounds, load ~25): flat-192.cif (124 MB) bytes 2,105 / 2,510 / 2,441 ms vs string
  1,260 / 1,778 / 990 ms; per needle 220–765 ms on bytes. 7qpd: 8.4–11.1 vs 2.6–3.6 ms. Raw throughput on the same file:
  `Uint8Array.indexOf` of an absent byte 410 ms, `String.includes` 298 ms, `Buffer.indexOf` 236 ms. The parse of the
  same file is ~0.3 s (v2, #2036), so the sniff costs several parses. The bench's `parse` stage calls `CIF.parse`
  directly and does not see it; Prism's file and URL loads (format `auto`) do.
- **Fix.** (S) search the server variants in `header` only (they are told at the start); (M) tell `sfcif`/`coreCif`
  after `ParseCif`, which every CIF provider runs anyway, from the parsed category and field names (O(categories)), or
  one pass over the bytes for all needles. Either removes the 5 full scans.

### 14/N1 Hover/select

- **Evidence.** `sceneUpdate(ROs)` per 21 hovers: 33,160 (N = 100), 66,600 (N = 200), 166,600 (N = 500), linear in the
  scene, equal to v2. Profile N = 200: `lociMarkProvider` 31.6 ms, of it canvas `mark` 30.7 ms over all representations;
  `clearHighlights` (the `EveryLoci` broadcast) 10.3 ms; `areRootsEquivalent` 9.4 ms. Fix as v2 (per-root index,
  clear only the marked ones).

### P10 Prism start closure (new)

- **Evidence.** Scratch production build of `apps/prism-viewer` with metafile: start closure 90 files, 5.290 MB raw /
  1.591 gzip / 1.334 brotli (budget 5.5 MB raw: 210 kB left). The tree after #2026 (`217376451`, same build script):
  86 files, 5.160 / 1.546 / 1.295: +130 kB raw (+45 kB gzip) in 3 days (prism/viewer +39 kB, model +18, core/util
  +18, plugin/core +13, io +11: #2043, #2039, #2016). Inputs exclusive to a subtree (closure minus closure without it):
  WebGPU 308 kB (`plugin/core/util/graphics-backends.ts` imports `gl/backends/all` and `gl/backends/webgpu`), MVS
  runtime+builder 234 kB (`features/loaders/operations.ts`, `components/formats.ts`), the bundled extensions 222 kB
  (`apps/prism-viewer/src/extensions.ts` via `full-spec`), transpilers 87 kB (`plugin/core/default-registry.ts` →
  `query-language/transpilers/all`, and `prism/viewer/targets/languages/all`), `ai/src/threads` 27 kB. Together 877 kB,
  17 % of start.
- **Fix.** Lazy where the feature is first used: MVS loading (the loader is async), the transpilers (a language is
  picked by the user; `DefaultRegistry` could register them as lazy entries), extensions that only add behaviors on a
  loaded structure. WebGPU: backend selection becomes async, which touches the frozen descriptor modules: hand to the
  WebGPU owners. Raise nothing in the budget until one of these lands: at +130 kB per 3 days the 5.5 MB budget fails
  within a week or two.

### P12 Quick lane: slow by nature and slow by waste (new)

- **Evidence.** Targeted run `vitest run packages/plugin scripts tests --reporter=json --maxWorkers=4` (86 files, 540
  tests, load ~25; the whole quick lane was not run: the machine did not allow it): summed file time 32.5 s.
  `scripts/workspace/_test/typecheck-full.test.mjs` 13.0 s (2 tests, 8.3 and 4.7 s; 5 `tsc` spawns), `process.test.mjs`
  2.1 s (child processes), `exports.test.mjs` 2.2 s (one test 2.1 s, spawns Node to resolve), then
  `trajectory.test.ts` 2.4 s (18 tests), `snapshot-restore` 1.2 s, `custom-state` 1.2 s. `structure-focus-representation.test.ts`
  0.83 s for 5 tests: `focusStructure` waits a fixed 100 ms plus 25 ms polls with real timers.
- **Fix.** The three spawning files are "CLIs in a process" (AGENTS.md: heavy) or at least medium: tag them, keep a
  quick in-process representative (the config logic without `tsc`). The focus tests wait for the focus update by event
  (`plugin.state.data.events.changed` / the update promise) instead of sleeping (slow by waste: fixed, not moved).

### P13 StructureInfo per root model (new, low)

- **Evidence.** many N = 200 `model+structure` profile: `structure-info.ts` frames 36.0 ms (`setMaxIndices` 17.7)
  of ~145 ms. Each created root model loops every root model and structure (3 loops plus `asymIdOffset` in
  `handleModel`): O(N) per model, O(N²) per load of N entries. Allocation is flat per structure (0.72 → 0.59 MB at
  N = 100 → 500), time is not.
- **Fix.** A raise of the max (the usual case: the new model has the next index) sets `MaxIndex` on every root, which
  is inherent while `MaxIndex` is a per-object property; a shared max holder read by the props (or setting the max once
  per transaction, at `changed`) makes it O(1) per model. Low priority: the model stage is ~3 % of a many load.

### Smaller notes

- Startup is compile-bound now: 481 of 816 ms profiled under load in `compileSourceTextModule` for a 6.8 MB bundle;
  `normalizeSource` regexes over every GLSL program at module load are ~17 ms of samples
  ([programs.ts:80](../../packages/graphics/src/gl/webgl/programs.ts#L80), frozen).
- Classic bundles (scratch production builds): Viewer `molstar.js` 5.079 MB / 1.467 gzip (v2 5.00 / 1.44),
  docking-viewer 3.908 / 1.134, mesoscale-explorer 3.966 / 1.145, mvs-stories 4.134 / 1.205, Prism single file 7.240 /
  2.221 (budget 8.0). Only Prism has budgets (H10).

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B1 | low | still open | Pipeline bench crashes on coarse-only entries | [bench.ts:324](../../tests/src/pipeline/bench.ts#L324) | `bench.ts --fixture data/examples/ihm-coarse-8zzh.cif`: TypeError "reading 'elements'" after `theme-eval:residue-name`, exit 1, no table |
| B4 | low | still open | Bench stages record errors as fast results | [stage.ts:99](../../tests/src/pipeline/stage.ts#L99) | a throwing stage is stored with its `ms` and `info: 'ERROR: …'`, exit 0: reads as a speedup in a comparison |
| 03/B7 | low | still open | Late abort discards a completed result undisposed | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) | abort requested while `task.f` finishes: `abort()` throws after the result exists; textures/assets of it are never disposed |

- **B1** reproduced at HEAD (`out/8zzh.log`): every stage through `theme-eval:residue-name` runs, then
  `selections()` takes `s.units.find(atomic polymer)!` (undefined on 8zzh). (`asa` no longer throws on 8zzh.)
- **B4, 03/B7**: code unchanged since v2 (`git diff 0944823d4..HEAD` empty for both lines).
- No new bug confirmed in the slice. Checked: `trackGroup` + `groupRef` (a second tagged group after the first is
  removed is found by the one search), StructureInfo's updated-root-model path, the hierarchy `outer` stack for MVS
  trajectory-in-model, `ResponseBody`/`letGoOfBytes` with two Downloads of one URL, DownloadBlob's `onText` asset
  lookup (finds the held asset through `tryFindUrl`), DeflateData pass-through of already inflated downloads,
  `LogSize` (the only readers of `entries.size` deltas are two tests with small logs).

## Already good, re-checked

- No O(N²) in allocation over N structures: repr:auto per structure 10.3 → 7.3 MB (N = 100), 9.1 → 5.9 MB (N = 500);
  `remove:all` 128 → 46 MB at N = 500.
- many N = 500 base → HEAD allocation: load 2,330 → 1,471 MB, bonds 1,280 → 843, model+structure 376 → 295, theme:all
  308 → 253, select 11 → 8. Retained after repr:auto unchanged (3,464 → 3,402 MB).
- Startup: heap and allocation halved (P1); `init` allocation unchanged (8–19 MB).

## Not audited

- Wall-clock A/B of the many bench and startup at v2's precision (load 14–34, swap full); a second N = 500 pair was
  stopped to keep memory for the other auditors.
- The whole quick lane (only the slice's 86 files); `packages/plugin/headless` (no change since v2; GL capture needs a
  browser or `gl`); browser frame time (07/P7, P4).

## Bench scripts

All under `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/12/`:

| Script | What |
| --- | --- |
| `build.mjs <tree> <entry> <out>` | bundles an entry against HEAD or an extracted tree (links `@molstar/*`), as `pipeline-bench.mjs` |
| `runmany.sh`, `runmany2.sh`, `out/many*.log`, `out/*-100-*.json`, `out/*-500-a.json`, `cmp.mjs` | many A/B (base `0944823d4` vs HEAD), N = 100 ×2, N = 500 ×1 |
| `prof200/`, `prof-top.mjs`, `callers.mjs` | many N = 200 CPU profiles (repr:auto, model+structure, theme:one, hover, select, remove:one) |
| `startup/run.mjs`, `startup/{base,head}.mjs` | import heap/alloc/time and init, 5 interleaved pairs; `startup/prof/` import profile |
| `bundle/app.mjs`, `bundle/meta.mjs`, `bundle/sizes.mjs`, `bundle2026/` | scratch production builds of the apps (HEAD and `217376451`); start closure and per-package bytes |
| `sniff.ts` | CIF variant sniff, bytes vs string, flat-192.cif and 7qpd |
| `out/8zzh.log` | B1 reproduction |
| `out/vitest-plugin.json` | targeted quick-lane durations |
