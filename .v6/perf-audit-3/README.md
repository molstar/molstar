# Performance audit v3: index

## Scope

Third run, after [v2](../perf-audit-2/README.md) (2026-10-22, base `0944823d4`).

- **PRs covered**: the ~25 PRs since v2, `git log 0944823d4..16277d822` (#2024–#2048).
  - v2's top 5 and app-level top 3.
  - Slice 11 fixed.
  - Five batches of small fixes.
  - Three leftovers PRs (core, model, graphics).
  - The CIF byte tokenizer and lazy BCIF decoding.
  - Periodic trajectories.
  - Labels.
  - Bond providers.
  - Prism sequences in the structure tree.
- **Questions**:
  - Re-classify every v2 item.
  - Review each fix since v2 for new bugs or costs.
  - Find the big wins left.
  - Record every real bug met.
- **Read-only**: nothing under `packages/` was changed.
- **Conditions**:
  - HEAD `16277d822` (= master), 2026-10-25.
  - Apple M4, 10 cores, 16 GB, Node 25/26.
  - A new slice, **13-io**, audits `packages/io`, which v2 only profiled.

**Load warning.** 14 auditors ran at once on a shared machine, at load 14–41 with 12–18 GB of swap in use.

- **Wall times** are 2–10× v2's and swing up to 6× between identical runs.
- **Evidence used instead**: allocation, retained bytes, buffer bytes, counters, CPU time, and same-process interleaved A/B against
  the v2 base (`git archive 0944823d4`, the same bench bundled against it).
- **Verify pass**: it ran at load 4–9. Its reruns lowered several wall-time estimates by 3–14× (01/P13, 08/P18, 09/P19, 09/P20,
  12/P12). The estimates in the tables below are the verified ones.

Constraints:

- Canvas3D and `graphics/src/gl` are frozen while WebGPU lands. Findings there are marked `frozen`.
- MVS loading and animation are in [mvs-optimizations](../mvs-optimizations/research/README.md).
- Ids are slice-qualified (`05/P13`). A v2 item keeps its id; new ones continue the slice's numbering.

## Since this audit

Re-checked against the code at `24aa41e43` (no reruns); the tables below are marked where an item changed. Every other
item stands as written.

- **#2050** (first batch) fixed 13/P1, 08/P12 + 08/P15 (v2's R1, R2), 07/P17, 01/P14, 05/P25, 05/P23, 02/B6, 02/B7,
  12/B1 and 12/B4. Partly fixed: 01/P13 = 12/P11 (3 scans of a plain mmCIF, not 5; the header-only sniff is open).
  Still open: 01/P12 (the anchors keep the builders' maps alive; the `ChunkedArray` stores stay polymorphic across
  array kinds).
- **#2051** (graphics) fixed:
  - 07/P7 (and H6): the merged object's geometry and theme layers are the members' arrays as parts, with no copy.
    `auto` from 30 units was measured on both backends ([merged-render-objects.md](../plans/gl/merged-render-objects.md)
    steps 1, 3, 4).
  - 12/V1: the 3D volume texture is read at `Tensor.strides`
    ([direct-volume-performance.md](../plans/gl/direct-volume-performance.md) step 7). 12/V2, V4 and L1 are still open.
  - 14/N1, the scene half: a changed mark updates only the marked render objects (`Scene.updateChanged`), not the
    whole scene. The many-structures bench stub counts the changed objects too (`sceneUpdated`; N = 100: 40 per 21
    hovers, one per changed mark). `Scene.update(objects)` now clears the changed flags of the objects it updates
    (`DrawTable.dropChanged`): before, the next mark updated them again (868 per 21 hovers after a recolor of all).
  - 14/N1, the representation half: a mark without a representation visits the representations of its structure
    loci's root only (`MarkIndex`, by `Representation.getStructureRootHashes`; other loci all of them). N = 200:
    `repr.mark` 534 → 33 per mark (the equivalent copies of the fixture), hover 1.7–2.1 → 0.40 ms and select
    1.8–3.7 → 0.25 ms per event (clearing highlights marks the previous loci, so indexed too). Deselect all
    (`EveryLoci`, an action that only removes) visits the representations markers were added to (`MarkIndex.note`):
    534 → 1, 12.6 → 0.6 ms, the same objects changed.
- **Fixed after #2051**: 10/P12 (Canvas3D side): `Canvas3D.update(repr)` records that representation's renderables
  again (`PrecompileRecorder.recordAgain`), not every renderable of the scene
  ([pipeline-readiness.md](../plans/gl/pipeline-readiness.md) step 6c). The callers still update one representation at
  a time.
- **Still open in Canvas3D and `gl`**:
  - 08/P7: cylinder attributes ×6; the plan for it is [primitive-expansion.md](../plans/gl/primitive-expansion.md),
    not started.

## Method

1. **Profile** ([00-profile.md](00-profile.md)): the pipeline bench as an A/B of the v2 base against HEAD, in one session.
2. **Per-slice audit** (01–13): each auditor did four things.
   - Re-classified every v2 item against the code and a rerun of its probe.
   - Diffed `0944823d4..HEAD` over its paths and read every changed hunk for new bugs or costs.
   - Hunted new findings, each with a bench or probe.
   - Wrote the slice report.
3. **Adversarial verify**: three batched verifiers (52 items: core and io, model, graphics/Prism/app) re-read each new finding
   with its callers and re-ran the probes.
   - Verdicts: [\_verify-A](_verify-A.md), [\_verify-B](_verify-B.md), [\_verify-C](_verify-C.md).
   - Refuted: 13/B4 and 04's `areBackboneConnected` note.
   - Lowered: 01/B17, 05/B15, 06/P14, 12/P10, 12/P12.
   - Merged: 01/P13 = 12/P11, 07/P20 = 09/P17, 13/B1 = 01/B10.
   - Added: 05/P25, from a cross-slice note.

## Reports

Counts are open items after verification (perf / bugs).

| Slice                                                 | Paths audited                                                                                                     | P / bugs |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | -------: |
| [00-profile](00-profile.md)                           | pipeline stages A/B v2 base → HEAD, startup, deopts                                                               |        – |
| [01-core-data-util](01-core-data-util.md)             | `packages/core/src/data`, `util`, `algorithms`                                                                    |    6 / 7 |
| [02-core-math](02-core-math.md)                       | `packages/core/src/math` (+ molsurf visual math)                                                                  |    7 / 2 |
| [03-core-state-task](03-core-state-task.md)           | `packages/core/src/state`, `task`, `util/async-queue.ts`                                                          |    7 / 5 |
| [04-model-formats](04-model-formats.md)               | `packages/model/src/formats`, `model/src/model` minus structure, kinemage reader                                  |   12 / 5 |
| [05-model-structure](05-model-structure.md)           | `packages/model/src/model/structure` (+ `model/periodic.ts`)                                                      |   12 / 5 |
| [06-model-props-script](06-model-props-script.md)     | `packages/model/src/props`, `script`                                                                              |   11 / 5 |
| [07-graphics-repr](07-graphics-repr.md)               | `packages/graphics/src/repr`                                                                                      |   11 / 3 |
| [08-graphics-geo](08-graphics-geo.md)                 | `packages/graphics/src/geo` (not GPU density, not `gl/`)                                                          |   10 / 3 |
| [09-graphics-theme-props](09-graphics-theme-props.md) | `packages/graphics/src/theme`, `props`                                                                            |    7 / 1 |
| [10-prism-viewer](10-prism-viewer.md)                 | `packages/prism/viewer/src` (+ Canvas3D.update callers; the ai thread record)                                     |    8 / 2 |
| [11-prism-ui-assistant](11-prism-ui-assistant.md)     | `packages/prism/ui/src`, `prism/assistant/src`, a light pass over `packages/ai`                                   |    1 / 0 |
| [12-app-level](12-app-level.md)                       | `plugin/core` managers and behaviors, render loop, `plugin/headless`, `apps/`, esbuild, pipeline bench, test lanes | 14 + 9 H / 2 |
| [13-io](13-io.md) (**new**)                           | `packages/io`: CIF text and binary, msgpack, CSV, PDB/GRO/XYZ/LAMMPS and trajectory readers, gzip                 |    6 / 3 |

Not audited:

- `packages/mvs` (covered elsewhere), `query-language`, `sequences` (Seq\* has its own perf tiers), most of `packages/ai`
  (providers, narration, compaction), servers, and most extensions.
- Not measured: browser frame time, GPU memory, the WebGPU backend, and idle-machine wall times.

## Profile headlines (v2 base → now, A/B)

See [00-profile.md](00-profile.md).

**Profile A/B** (ms / MB alloc, ranges): the four small fixtures were rerun at load 5–12. Node is v26.12; v2 ran on v25.6.

- **`interactions`:**
  - flat-32: 906 → 364 ms, 579 → 194 MB.
  - asm-192: 2489 → 1065 ms, 1478 → 448 MB.
  - flat-192: 5473 → 2091 ms, 3023 → 893 MB.
- **`repr:molecular-surface`:**
  - flat-32: 8.7 → 6.2 s.
  - 7qpd: 1747 → 1337 ms.
  - 8zzh: 1829 → 632 ms.
- **Transformed structures:** `transform:bonds` on asm-192 295 → 0.7 ms; `transform:polymer-bonds` 234 → 12 ms (05/P14).
- **`model` on flat-192:** 563 → 384 ms, 520 → 152 MB alloc, 343 → 133 MB retained.
- **Text CIF `parse` on 7qpd:** 23 → 9.8 ms.
- **`lookup+boundary`:** asm-192 49 → 15 ms; flat-192 376 → 216 ms.
- **`bonds-inter`:** asm-192 288 → 168 ms; flat-192 272 → 144 ms.
- **`asa` on asm-192:** 6.0 → 6.1 s, alloc 613 → 121 MB. Still the largest stage: every copy is computed (13/P6).
- **`many`, N = 500:**
  - `represet`: 4147 → 235 ms.
  - `repr:auto`: 12.8 → 9.3 s (25.7 → 18.5 ms per structure, still super-linear).
- **Startup:** import 290–370 → 128–176 ms.
- **Unchanged:** render buffers and playback counters are identical on every fixture.

Allocation from the slice A/Bs, HEAD against the v2 base in one session:

| Stage                                | v2 base → now                                                                                               | Main cost now                                                                  |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `repr:molecular-surface`             | flat-32 alloc 1.51–1.55 GB → 0.21–0.34 GB (02/P10, #2025)                                                    | torus pass 57–76 % of surface CPU at every resolution (02/P12)                 |
| `interactions` (prop)                | flat-32 alloc 477–484 → 100–102 MB, time 2.5–3.1 → 0.83–1.0 s (loaded); asm-32 237–257 → 61–81 MB; **edges identical** | per-symmetry-pair inter loop (06/P3), line of sight 21–24 % (06/P11)            |
| `repr:interactions`                  | flat-32 alloc 148–150 → 91–100 MB; asm-192 272–278 → 128–133 MB                                             | asm: one non-instanced inter visual, 265 MB buffers on asm-192 (09/P17)        |
| `repr:ball-and-stick`                | flat-32 alloc 400–455 → 196–244 MB, buffers identical                                                        | merged copy beside members (07/P7)                                             |
| `bonds-inter` / `bonds-intra`        | asm-192 CPU 539–577 → 364–369 ms, alloc 95 → 36–39 MB; flat-32 intra 348–406 → 272–336 ms CPU                | per-symmetry-pair (05/P13), each pair found twice (05/P21)                     |
| `lookup+boundary`                    | flat-32 alloc 41 → 5–16 MB; asm-192 51 → 6 MB                                                               | 98-direction unit boundaries (05/P2)                                           |
| `model` (BCIF)                       | flat-192 alloc 487–550 → 112–179 MB, retained 320 → 124 B/atom; multi-model text 153–159 → 42–47 MB          | 124 B/atom retained (04/P16), `missingResidues` eager (04/P15)                 |
| text CIF parse                       | flat-192.cif ~2× faster (same process, #2036/#2040); decode hashes identical                                 | ByteTokens growth (13/P5); `auto` variant sniff 5 whole scans (01/P13)         |
| BCIF parse                           | flat-192.bcif retained after the model's reads 118 → 72 MB                                                  | msgpack `bin` copied byte by byte: ~80 % of `parseBinary` (13/P1)              |
| `repr:cartoon`, `repr:gaussian-surface` (7qpd) | v2's R1/R2 **still there**: alloc 71 / 77 MB, 33 / 45 MB with a module anchor                     | optimized code dropped at major GC (08/P12, 08/P15)                            |
| DSSP, valence (flat-32)              | alloc 83–97 → 14–16 MB and 79 → 52 MB, identical output                                                     | valence per unit, not per group (06/P8)                                        |
| many structures, N = 500             | `repr:auto` alloc 4.5 → 2.9 GB (−35 %), load 2.3 → 1.5 GB, bonds 1.28 → 0.84 GB                              | hover 7,933 render-object updates (14/N1); hierarchy walk ≈ all of an action (14/N2) |
| state commits                        | 2-removal commit at 10k nodes 6.3 → 0.018 ms; N = 1000 per-structure removals 51 s → 50 ms (03/P1, #2025)    | snapshot frames 24–26 ms / 28 MB at 10k (03/P15)                               |
| startup                              | heap after import 34.9 → 13.0 MB; Prism start closure 5.16 → 5.29 MB raw (96 % of the 5.5 MB budget)         | compile-bound 6.8 MB bundle; ~0.4–0.57 MB deferrable (12/P10)                  |

Biggest gains since v2:

- Molecular-surface allocation ÷5–7.
- Interactions allocation ÷4–5 and time ÷3, with identical edges.
- BCIF model ÷3–4 allocation and −60 % retained per atom.
- Text CIF parse ~2×.
- Large commit removals from seconds to ms.
- `many` −35 % allocation.
- Startup heap ÷2.7.

### Regressions

**07/P17** (#2048 merge fingerprint) and the new costs from fixes listed below.

**From v2:**

- **R1** (08/P12): 7qpd cartoon still 79 MB alloc against v1's 38. **fixed in #2050**.
- **R2** (08/P15): gaussian 69 → 71 MB. Both are unchanged against the v2 base, and `--trace-deopt` shows the same drops. **fixed in #2050**.
- **R3:** dropped.

**New profile leads** (00 only, not bisected, not verified):

- **R4:** 3k7f molecular-surface allocation ~140 → ~320 MB; the time is not worse.
  - `projectPointsRange` self allocation goes 52 → 120 MB.
  - Candidates: #2034 (`updateChunk`, compact Float64 coordinates) or #2033.
- **R5:** 1tup-2x2x2 gaussian surface 28–33 → 43–51 MB alloc (+50 %).
  - MC `interpolate` and `accumulateRange` allocate more.
  - Candidate: #2034's `updateChunk` in `gaussian-density/cpu.ts`.
- **R6 (small):** `bonds-inter` now retains 0.3 → 17 MB on flat-192, and `interactions` on flat-32 retains +11 MB.
  - Candidates: #2045's struct_conn partner cache, or the #2016 bond providers.

**L1:** allocation is not load-proof either. 7qpd molsurf allocated 115–478 MB per run at load 31, and 146–181 MB at load 6.
Any clock-free guard (12/H10) needs medians at low load.

## What changed since v2

v2 columns count the v2 items each slice re-classified. "New" counts new findings after verification.

| Slice     | Fixed | Partly | Open | Regressed | Dropped | New P | New bugs | Main PRs                                                               |
| --------- | ----: | -----: | ---: | --------: | ------: | ----: | -------: | ---------------------------------------------------------------------- |
| 01        |    16 |      1 |   10 |         0 |       0 |     3 |        2 | #2032, #2033, #2036, #2044, #2045, #2046                               |
| 02        |     9 |      2 |    6 |         0 |       1 |     1 |        2 | #2025, #2033, #2034, #2046                                             |
| 03        |    12 |      2 |    9 |         0 |       0 |     1 |        0 | #2025, #2033, #2034, #2045, #2046                                      |
| 04        |    15 |      4 |    9 |         0 |       0 |     2 |        1 | #2032, #2033, #2036, #2044, #2045, #2047                               |
| 05        |    14 |      2 |   11 |         0 |       0 |     3 |        1 | #2025, #2032, #2034, #2044, #2045, #2047                               |
| 06        |    12 |      2 |   12 |         0 |       0 |     1 |        1 | #2025, #2032–#2034, #2044, #2045, #2047                                |
| 07        |    19 |      0 |   13 |         1 |       0 |     0 |        0 | #2025, #2032–#2034, #2044, #2045, #2048                                |
| 08        |     9 |      2 |    9 |         0 |       0 |     3 |        1 | #2028, #2030, #2033, #2034, #2044, #2048                               |
| 09        |     9 |      1 |    4 |         0 |       0 |     2 |        1 | #2025, #2032–#2034, #2042, #2044, #2045, #2048                         |
| 10        |     5 |      1 |    5 |         0 |       0 |     3 |        1 | #2026, #2032, #2034, #2044                                             |
| 11        |    12 |      0 |    0 |         0 |       0 |     1 |        0 | #2029 (all)                                                            |
| 12        |    13 |      2 |   12 |         0 |       1 |     3 |        0 | #2026, #2034, #2036, #2044, #2045, #2048                               |
| 13        |     2 |      1 |    3 |         0 |       0 |     5 |        2 | #2036, #2039, #2040                                                    |
| **Total** | **147** | **20** | **103** | **1** | **2** | **28** | **12** |                                                                        |

v2 had 144 open and 22 partly fixed items, plus 61 new ones. About two thirds of them are now closed.

### New bugs or costs brought by fixes since v2

v2 found 8 such cases; this run found 12. The pattern is the same: no counter or exactness probe in CI caught them.

| From  | Item                                                                                                              |
| ----- | ----------------------------------------------------------------------------------------------------------------- |
| #2034 | **02/B6** (medium): molsurf chunk size counts voxels only (the 02/B2 fix)                                         |
| #2048 | **07/P17**: merge fingerprint watches `uLod` / `uClipObjectType`, re-versioned on every update                     |
| #2025 | 05/P24: the P14 transform fix sends the first bonds call of a small moved component through the root (~0.25 s once) |
| #2044 | 06/P14: the B6 fix widens the ASA lookup radius for every atom (opt-in `nonPolymer` only); 01/B17 doubled         |
| #2045 | 03/P18: chunk sizes collapse when several chunked tasks interleave (not reachable in one plugin)                  |
| #2046 | 09/P19: the new 16-entry distinct-colors LRU thrashes past 16 compartments (before: no cache at all)               |
| #2039 | 01/P13: the CIF variant sniff on bytes, +0.11 s on flat-192.cif; 01/P15 `HeldAsset.file` copies                   |
| #2037 | 02/B7, 05/B15, 05/P23, 04/P20, 04/P21                                                                             |
| #2040 | 13/B5, 13/P4 (`shortAscii` concatenation, also #2036)                                                             |
| #2030 | 08/P20: kerned layout allocates per label (+48 MB per 50k labels)                                                 |
| #2043 | 10/P17, 10/P19, 10/B7, 01/P14                                                                                     |
| #2029 | 11/P12: the pinch-guard fix (11/P10) still covers the turn card's scroll body                                     |

## Top open wins, ranked

Rank = impact × breadth ÷ effort; v2's still-open items and new ones together. Estimates are the verified ones.

| #   | Id                           | Status vs v2                | Title                                                                                                    | File                                                                                                                                                                         | Kind        | Estimate                                                                                                   | Effort | Frozen |
| --- | ---------------------------- | --------------------------- | -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ---------------------------------------------------------------------------------------------------------- | ------ | ------ |
| 1   | 13/P1                        | new; **fixed in #2050** | MessagePack `bin` copied byte by byte                                                                    | [decode.ts:121](../../packages/io/src/common/msgpack/decode.ts#L121)                                                                                                        | time        | measured: flat-192.bcif `decodeMsgPack` 33 → 3 ms (~80 % of `parseBinary`); every BCIF load; output identical | S      | no     |
| 2   | 08/P12 + 08/P15 (+ 01/P12)   | still open (v2 R1/R2); **fixed in #2050** but 01/P12 | Optimized builder and marching-cubes code dropped at major GC (per-run maps; ChunkedArray adders too)    | [mesh-builder.ts:120](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L120), [algorithm.ts:117](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L117), [chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59) | time+memory | measured (7qpd): cartoon alloc 71 → 33 MB, gaussian 77 → 45 MB; 3000 tubes cold 29–46 → 6.5–8.5 ms        | S      | no     |
| 3   | 07/P17                       | **regressed** (#2048); **fixed in #2050** | Merge fingerprint misses on every update                                                                 | [merged.ts:30](../../packages/graphics/src/repr/merged.ts#L30)                                                                                                              | time        | every recolor re-checks every member (416/416 on flat-32); ~35 ms of a 1.5 s flat-192 recolor               | S      | no     |
| 4   | 05/P13                       | still open                  | Inter-unit bonds recomputed per symmetry-equivalent unit pair                                            | [inter-compute.ts:406](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L406)                                                                 | time        | prototype: asm-192 ÷2.8–3.8, edges equal (30,906 pairs in 241 classes)                                     | M      | no     |
| 5   | 01/P13 = 12/P11              | new (#2039); **partly fixed in #2050** (3 scans) | CIF variant sniff scans the whole file 5× on `auto` open                                                 | [provider.ts:186](../../packages/plugin/core/src/state/formats/provider.ts#L186), [utf8.ts:271](../../packages/core/src/util/utf8.ts#L271)                                 | time        | ~0.27 s per open of flat-192.cif (124 MB), about the parse itself; 3 scans exact (−40 %), or header sniff   | S      | no     |
| 6   | 01/P14                       | new (#2043); **fixed in #2050** | `numericRange` builds a string `Set` to find min/max                                                     | [value-palette.ts:641](../../packages/core/src/util/color/value-palette.ts#L641)                                                                                            | time+memory | ~0.2 s and 91 MB per 1.1M values (Seq\* rows) → 4–17 ms                                                     | S      | no     |
| 7   | 05/P25                       | new (verify); **fixed in #2050** | Always-true `areIntersecting` per residue in `getAtomicPolymerElements`                                  | [polymer.ts:25](../../packages/model/src/model/structure/structure/util/polymer.ts#L25)                                                                                    | time        | −40–60 % of the function (flat-32 11–12 → 4.5–7 ms), identical output                                     | S      | no     |
| 8   | 05/P23                       | new (#2037); **fixed in #2050** | Periodic whole frames held by a module-level `Map`                                                       | [periodic.ts:65](../../packages/model/src/model/structure/model/periodic.ts#L65)                                                                                           | memory      | 83 MB at 100k atoms, 191 MB at 250k retained after every reference is dropped                              | S      | no     |
| 9   | 07/P2                        | still open                  | Segment visuals size markers and layers to the whole grid, per segment                                   | [segment.ts:221](../../packages/graphics/src/repr/volume/segment.ts#L221)                                                                                                   | memory+time | −152.6 MB for 20 segments in 200³; −44 ms per hover pair                                                  | S      | no     |
| 10  | 14/N2 + 14/N1                | 14/N2 still open; 14/N1 fixed in Canvas3D (#2051 scene half, `MarkIndex`) | Hierarchy rebuilt by a whole-state walk per change; hover marks every representation                     | [hierarchy.ts:253](../../packages/plugin/core/src/state/manager/structure/hierarchy.ts#L253), [representation.ts:54](../../packages/plugin/core/src/behavior/dynamic/representation.ts#L54) | time        | N2 8.2–8.5 ms of a 9–10 ms action at N = 200; N1 7,933 render-object updates per hover at N = 500          | M      | no     |
| 11  | 05/P18                       | still open                  | Fixed-point `includeConnected` is quadratic                                                              | [modifiers.ts:358](../../packages/model/src/model/structure/query/queries/modifiers.ts#L358)                                                                               | time        | 7qpd ~0.44–0.97 s → a few ms (frontier BFS)                                                                | M      | no     |
| 12  | 02/P2 + 02/P3                | still open                  | Grid query: bucket-ordered copies; 32 elements per cell                                                  | [grid.ts:360](../../packages/core/src/math/geometry/lookup3d/grid.ts#L360)                                                                                                 | time        | query −22–27 % (+12 B per element), −10–30 % more with 8 per cell (order changes)                          | S      | no     |
| 13  | 04/P20 + 13/P2               | new (#2037) + still open    | LAMMPS: eager unwrapped coordinates per frame; token store per column and frame                           | [lammps-trajectory.ts:111](../../packages/model/src/formats/structure/lammps-trajectory.ts#L111), [traj/parser.ts:47](../../packages/io/src/reader/lammps/traj/parser.ts#L47) | memory      | −12 B/atom/frame (−1.2 GB at 100k × 1000); token stores 47.5–71.5 B per atom-frame                          | S / M  | no     |
| 14  | 09/P19 + 09/P20              | new                         | Particle themes: distinct colors thrash a 16-entry LRU; compartment index rebuilt per `Theme.create`     | [particle-hierarchy.ts:239](../../packages/graphics/src/theme/color/particle-hierarchy.ts#L239), [particle-compartment.ts:75](../../packages/graphics/src/theme/color/particle-compartment.ts#L75) | time        | 17 × 100 compartments 33–45 ms per update → 0; 10M particles ~46 ms + 40 MB per update → 0                 | S      | no     |
| 15  | 03/P15                       | still open                  | Every snapshot-transition frame rebuilds (`fromJSON`) and diffs the whole tree                           | [state.ts:233](../../packages/core/src/state/state.ts#L233)                                                                                                                | time+memory | 24–26 ms and 28 MB per frame at 10k nodes → 0.1–0.3 ms with a frame-delta path                            | M      | no     |
| 16  | 09/P17 (= 07/P20)            | still open                  | Inter-unit interactions on assemblies: one non-instanced visual                                          | [interactions-inter-unit-cylinder.ts:326](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L326)                             | memory+time | asm-192 265 MB buffers → ~11–20 MB (edges between copies stay)                                             | L      | no     |
| 17  | 05/P2                        | still open                  | Unit boundaries use 98 directions up to 10k atoms                                                        | [unit.ts:435](../../packages/model/src/model/structure/structure/unit.ts#L435)                                                                                             | time        | extrema 45 % self of the flat-32 stage; `lookup+boundary` ~÷3 (v2)                                         | S      | no     |
| 18  | 02/P12 (+ 02/P14)            | still open                  | Molecular-surface torus pass ignores resolution; second `find` per atom                                  | [molecular-surface.ts:423](../../packages/core/src/math/geometry/molecular-surface.ts#L423)                                                                                | time        | torus pass 57–76 % of surface CPU; coarse skip `lowest` −75 % (v2); P14 exact early exit −6 % at res 8     | M      | no     |
| 19  | 10/P18                       | new                         | ai session record: each append copies and diffs the whole thread (O(E²) per session)                     | [state.ts:51](../../packages/ai/src/threads/state.ts#L51), [host.ts:166](../../packages/ai/src/threads/host.ts#L166)                                                       | time+memory | 285 KB / 0.1 ms per append at 100 entries → 2.1 MB / 0.83 ms at 4000                                       | M      | no     |
| 20  | 05/P21                       | still open                  | Intra `findBonds` finds every pair twice through per-atom queries                                        | [intra-compute.ts:387](../../packages/model/src/model/structure/structure/unit/bonds/intra-compute.ts#L387)                                                                | time        | pair enumeration ÷2.8; `bonds-intra` −25–40 %                                                              | M      | no     |
| 21  | 04/P16                       | partly fixed                | BCIF model retains 124 B/atom: occupancy/B Float64 twice, Cartn twice when sorted, `type_symbol` JS array | [atomic.ts:145](../../packages/model/src/formats/structure/basic/atomic.ts#L145)                                                                                          | memory      | −40 B/atom (−45 MB at 1.1M, −200 MB at 5M), estimated                                                      | S–M    | no     |
| 22  | 07/P7                        | **fixed in #2051**: no copy, `auto` measured ([merged-render-objects.md](../plans/gl/merged-render-objects.md)) | `renderMerged: 'auto'` by count; members keep their arrays beside the merged copy                         | [merged.ts:107](../../packages/graphics/src/repr/merged.ts#L107), [base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127)                                     | memory      | flat-32: b&s +168 MB beside 165 MB merged; ~2× retained; needs a GPU draw-call check                       | S      | release: yes |
| 23  | 12/P10                       | new                         | Prism start closure at 96 % of budget; MVS, extension bodies, transpilers, ai threads load at start       | [compose-plugin.ts](../../packages/prism/viewer/src/sources/compose-plugin.ts)                                                                                              | bundle      | ~0.4–0.57 MB deferrable (WebGPU stays: default backend `'auto'`); grew 130 kB in 3 days                    | S–M    | no     |
| 24  | 06/P3 + 06/P11 + 06/P7       | still open                  | Inter contacts per symmetry pair; line of sight and ASA through the structure lookup                     | [interactions.ts:464](../../packages/model/src/props/computed/interactions/interactions.ts#L464), [area.ts:98](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L98) | time        | P3 search reuse ~26 % of interactions; a shared heavy-atom grid pays off only across P11 + P7              | M–L    | no     |
| 25  | 10/P12                       | still open                  | `Canvas3D.update(repr)` re-records every renderable's signature; callers loop per representation          | [canvas3d.ts:2240](../../packages/graphics/src/canvas3d/canvas3d.ts#L2240)                                                                                                 | time        | drag at N = 100 ~12 → < 0.5 ms (v2 estimate); caller-side batching is not frozen                           | S      | **yes** |

Next tier, all S, in the slice reports:

- From the new code: 10/P17 (sequence rows of every structure per op), 10/P19, 08/P20 (labels +48 MB per 50k), 04/P21 (NAMD XSC
  quadratic), 05/P24, 13/P4, 13/P5, 13/P6, 01/P15, 11/P12 (pinch guard over the turn card).
- Carried over: 07/P12 (serial id field, 256 MB at 400³), 09/P4, 09/P3, 08/P18 (glyph SDF 5.5 ms each), 08/P19, 03/P16,
  10/P15, 10/P13, 12/P13 (StructureInfo O(N²) per load), 12/P12 (quick-lane `tsc` spawns and real-timer waits).

## Cross-cutting themes

| Theme                                         | Findings                                                                                          | Pattern                                                                                                     |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| Fixes that brought new bugs or costs          | 02/B6, 07/P17, 05/P24, 06/P14, 03/P18, 09/P19, 01/B17 (×2)                                        | a fix measured for its own case, with no counter for the neighbour case (chunk work, updates, small children) |
| New code without the v2 lessons               | #2037: 02/B7, 05/B15, 05/P23, 04/P20, 04/P21; #2043: 10/P17, 10/P19, 10/B7, 01/P14; #2030: 08/P20 | module-level caches with no bound; eager per-frame copies; whole-scene refresh per op; string keys per value |
| V8 maps and elements kinds                    | 08/P12, 08/P15, 01/P12 = 08/P17, 02/P13, 07/P19                                                   | per-run object maps embedded in optimized code; anchor the maps at module level                              |
| Byte paths slower than the string paths       | 01/P13, 13/P1, 13/P4, 13/P5                                                                       | per-byte loops where `indexOf`/`slice`/`fromCharCode` over a range would do                                  |
| Symmetry and sharing lost                     | 05/P13, 06/P3, 09/P17, 06/P8, 13/P6, 04/P2, 04/P3                                                 | per unit pair, copy or frame where a symmetry class or shared topology would do                             |
| Whole-state walks and per-op fan-out          | 14/N2, 14/N1, 03/P2, 03/P15, 10/P17, 10/P18, 10/P12, 12/P13                                       | every action rebuilds or scans everything                                                                    |
| Caches that never evict, leaks                | 01/B1, 05/P23, 03/B6, 03/B7, 08/B10                                                               | module-level Maps keyed by data; no dispose on error and abort                                               |
| Same root cause found twice                   | 01/P13 = 12/P11, 07/P20 = 09/P17, 13/B1 = 01/B10, 01/P12 = 08/P17, 01/B1 = 04/B1                  | fix once, in the owning package                                                                              |

## App-level hypotheses ([12-app-level](12-app-level.md))

Re-ranked by slice 12:

1. **H9**: event storms and whole-state walks (14/N2, 14/N1).
2. **H10**: the regression guard. Only partly done: clock-free counters and the closure budget landed in #2026, but no bench runs in
   CI and there are no baselines. The 12 fix-induced items above make the case again. Bench bugs 12/B1 and 12/B4 still block it.
3. **H7**: startup and bundle (12/P10; the budget is at 96 %).
4. **H6**: CPU copies after upload (07/P7; fixed in #2051).
5. **H8**: idle render loop (12/P4).
6. **H3**: latest-wins and abort.
7. **H1**: worker pool.
8. **H2**: LOD by symmetry group.
9. **H4**: persistent prop cache.

H5 (byte CIF tokenizer) is done (#2036, #2040).

## Order of work

| Step | What                                                                                                                                                                                                                        | Measure                                                                                            |
| ---- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| 0    | Bisect the profile leads R4, R5, R6 on a quiet machine. R4 and R5 point at #2034's `updateChunk`, as 02/B6 does | 3k7f molsurf, 1tup-2x2x2 gaussian, flat-192 `bonds-inter` retained |
| 1    | Exact S fixes:<br>- the fix-induced items: 07/P17, 02/B6, 05/P24, 06/P14, 09/P19<br>- new code: 13/P1, 01/P13, 01/P14, 05/P25, 05/P23, 09/P20, 04/P20, 04/P21, 08/P20, 10/P17, 11/P12<br>- the module anchors for 08/P12 + 08/P15 + 01/P12 | allocation A/B (7qpd cartoon/gaussian), BCIF parse, `auto` open of flat-192.cif, recolor merge-check count |
| 2    | New bugs:<br>- 02/B7 (27 shifts)<br>- 05/B15 (drop coordinate-derived data in `withAtomicConformation`)<br>- 06/B13, 04/B15, 10/B7, 08/B10, 09/B8, 13/B3<br>Then v2's medium bugs, all still open: 05/B14 (needs a definition), 08/B8, 04/B8, 04/B10, 01/B1, 10/B4, 03/B6, 03/B11, 06/B2 | probes from the slice reports become tests                                                          |
| 3    | Regression guard (12/H10):<br>- fix 12/B1 and 12/B4<br>- counters for merge checks, chunk counts, edges, retained per atom<br>- a nightly bench with baselines                                                              | the guard flags 07/P17 and 02/B6 when reverted                                                      |
| 4    | Algorithmic M with identical output: 05/P13 (then 06/P3), 05/P18, 05/P21, 03/P15, 14/N2 + 14/N1, 10/P18 (in `packages/ai`)                                                                                                  | `bonds-inter` on asm, preset sync and hover at N = 500, snapshot frames at 10k                      |
| 5    | Decisions and memory: 02/P2 + 02/P3, 05/P2, 04/P16, 07/P2, 09/P17 (L), 12/P10 (lazy MVS / extension bodies); hand 10/P12 and 08/P7 ([primitive-expansion.md](../plans/gl/primitive-expansion.md)) to the WebGPU owners (07/P7: fixed in #2051)                        | retained MB, render-buffer bytes, start closure bytes                                               |

## Open questions for the maintainers

- **Carried over from v2, unanswered:**
  - 02/P2 memory against time (+12 B per element for −22–27 % query).
  - 05/B14: the meaning of `within` / PyMOL `around`.
  - 04/B10 sigma as the standard deviation, which moves the default iso levels.
  - 05/P13 and 06/P3: per-symmetry-class reuse is exact for rigid operators only.
  - 03/B11.
  - 10/P12: caller-side batching now, or after WebGPU?
  - 12/H8: `PauseLoopWhenIdle` on by default?
  - 07/P7: who runs the GPU draw-call check? (Answered: measured in
    [merged-render-objects.md](../plans/gl/merged-render-objects.md) step 1.)
- **12/P10:** the Prism start-closure budget has 210 kB left. Should MVS and the extension bodies load lazily now? Or should
  the budget be raised, with a recorded reason?
- **12/H10:** v2 asked where the regression guard lives. This run counts 12 bugs or costs brought by fixes since v2, 8 in v2.
  Is a nightly heavy-lane bench with `.v6/baselines` the next step?
- **02/B7:** fix `minimumImage` by searching the 27 neighbouring lattice shifts, or document it as approximate for triclinic boxes?
- **05/B15:** `withAtomicConformation` copies dynamic property data on purpose, for speed. Should the data be split into coordinate-derived and
  topology-derived parts, or should periodic frames drop it?
