# Performance audit: index

Goal: find the big performance wins (algorithmic complexity, per-element allocation, recomputation that could be cached
or incremental, oversized buffers, blocking main-thread work) and the real bugs met along the way, across the v6 library
layers and Prism, ranked so they can be worked through. Read-only: nothing under `packages/` was changed; every number
below comes from scratch benches on `v6-opts-261019` (`90d717f49`), Node 25.6, Apple M4, headless.

Constraints: Canvas3D and `packages/graphics/src/gl` are frozen while WebGPU lands (findings there are recorded with
`frozen`); MVS loading and animation are in [mvs-optimizations](../mvs-optimizations/research/README.md) and not
repeated. Bugs are listed in [bugs.md](bugs.md).

## Method

1. **Profile** ([00-profile.md](00-profile.md)): the whole pipeline from file bytes (parse → model → structure → bonds →
   props → representations → themes) on six fixtures, 7qpd (11.5k atoms) to `big-bcif` (186k atoms, 416 unshared units);
   time, allocation, retention, CPU profiles mapped to TS, heap snapshots.
2. **Per-slice audit**: one auditor per slice reads the code with the profile, writes findings (P = perf, B/N = bug)
   with a bench or probe per claim.
3. **Adversarial verify**: a second pass re-reads each finding with its callers and re-runs the benches. Refuted
   findings are removed; estimates are the verifier's (several were cut, e.g. 01/P6, 03/B8, 05/P2, 07/P6).

## Reports

| Slice                                                 | Paths audited                                                              | P / bugs |
| ----------------------------------------------------- | -------------------------------------------------------------------------- | -------: |
| [00-profile](00-profile.md)                           | pipeline stages, fixtures, hotspots                                        |        – |
| [01-core-data-util](01-core-data-util.md)             | `packages/core/src/data`, `packages/core/src/util`                         |   9 / 11 |
| [02-core-math](02-core-math.md)                       | `packages/core/src/math`                                                   |    9 / 4 |
| [03-core-state-task](03-core-state-task.md)           | `packages/core/src/state`, `packages/core/src/task`                        |   9 / 10 |
| [04-model-formats](04-model-formats.md)               | `packages/model/src/formats`, `packages/model/src/model` (minus structure) |   14 / 7 |
| [05-model-structure](05-model-structure.md)           | `packages/model/src/model/structure`                                       |   11 / 8 |
| [06-model-props-script](06-model-props-script.md)     | `packages/model/src/props`, `packages/model/src/script`                    |  10 / 11 |
| [07-graphics-repr](07-graphics-repr.md)               | `packages/graphics/src/repr`                                               |  10 / 11 |
| [08-graphics-geo](08-graphics-geo.md)                 | `packages/graphics/src/geo` (not GPU gaussian density)                     |   11 / 5 |
| [09-graphics-theme-props](09-graphics-theme-props.md) | `packages/graphics/src/theme`, `packages/graphics/src/props`               |    6 / 4 |
| [10-prism-viewer](10-prism-viewer.md)                 | `packages/prism/viewer/src` (plus two plugin hot spots on its paths)       |   11 / 3 |
| [11-prism-ui-assistant](11-prism-ui-assistant.md)     | `packages/prism/ui/src`, `packages/prism/assistant/src`                    |    7 / 2 |
| [12-tensor-access](12-tensor-access.md)               | every `Tensor` `get`/`dataOffset`/`getCoords` site; strides (#1987)         |    9 / 2 |

Not audited: `packages/io` (profiled only: parse is cheap), `graphics/src/gl` and Canvas3D (frozen), `packages/plugin`
(only 10/P4, 10/P5), `packages/mvs`, `query-language`, `sequences`, `ai`, extensions, apps, servers. Not measured:
browser frame time, GPU memory, the WebGPU backend.

## Profile headlines

`big-bcif` (186k atoms, 416 units) unless noted; ms / MB.

| Stage                          | Cost                                            | Main cause                                                 |
| ------------------------------ | ----------------------------------------------- | ---------------------------------------------------------- |
| `repr:molecular-surface`       | **12.9 s**; 7qpd 1.7 s                          | `obscured` rescans the neighbour list per probe/voxel      |
| `asa`                          | 1.8–1.9 s, **7.5–8.3 GB** allocated             | per-atom arrays of arrays + comparator sort                |
| `repr:interactions`            | 1.1 s, 1.8 GB alloc, **422 MB** buffers         | dashed mesh cylinders; `Vec3.fromArray` megamorphic (41 %) |
| `interactions` (prop)          | 890 ms                                          | inter-unit scan, grid queries                              |
| `theme-apply` on surfaces      | 808 ms per recolor                              | color smoothing splat recomputed                           |
| `repr:gaussian-surface` (CPU)  | 806 ms; 7qpd 612 ms                             | marching cubes reads via `Tensor` `get` closure (45 %)     |
| `repr:ball-and-stick`          | 456 ms, 165 MB buffers (333 retained)           | cylinder impostor data ×6 per vertex (frozen); merged copy |
| `bonds-intra`                  | 242 ms                                          | `GridLookup3D.query` 48 %                                  |
| merged render objects          | retained ≈ 2.2× buffers                         | members keep arrays (frozen `gl/renderable/merged.ts`)     |
| mesh fallback (no `fragDepth`) | ball-and-stick 230 MB, spacefill 114 MB on 7qpd | headless/Node capture pays 20–200×                         |
| parse                          | 5 ms (BCIF), 90 ms (20 MB text)                 | fine                                                       |
| theme eval                     | ≤ 26 ms per 186k elements                       | fine                                                       |

## Top wins, ranked

Rank = impact × breadth ÷ effort. Several rows combine the same fix found from two slices. Estimates are verified.

| #   | Slice  | Id         | Title                                                                                 | File                                                                                                                                                                       | Kind        | Estimate                                                                             | Effort       | Frozen |
| --- | ------ | ---------- | ------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ------------------------------------------------------------------------------------ | ------------ | ------ |
| 1   | 02, 08 | P5, P2     | Marching cubes: direct strided reads instead of `Tensor` `get` closure, corner reuse  | [tensor.ts:292](../../packages/core/src/math/linear-algebra/tensor.ts#L292), [algorithm.ts:250](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L250)     | time        | MC 1.7× (4–7× polymorphic); gaussian surface CPU −25–30 %; identical output          | S            | no     |
| 2   | 02, 08 | P4, P4     | Inline array reads in `MeshBuilder.addPrimitive` (megamorphic `Vec3.fromArray`)       | [mesh-builder.ts:122](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L122)                                                                                  | time        | addPrimitive 1.5–3×; mesh interactions −30–50 % (1.1 → ~0.6 s)                       | S            | no     |
| 3   | 06     | P1         | ASA `computeRange`: flat neighbour buffers, no per-atom sort closure                  | [area.ts:64](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L64)                                                                    | time+memory | 1.6–2.4×; 8 GB → <10 MB allocated at 186k; identical areas                           | S            | no     |
| 4   | 02     | P1         | Molecular surface: compact neighbours, prune per torus, test occlusion only on change | [molecular-surface.ts:87](../../packages/core/src/math/geometry/molecular-surface.ts#L87)                                                                                  | time        | 1.4–1.7×; big-bcif 12.5 → 7.6 s; bit-identical                                       | M            | no     |
| 5   | 08     | P1         | Mesh color smoothing kernel (no sqrt on 216 cells, AABB grid, z inner)                | [color-smoothing.ts:128](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L128)                                                                            | time+memory | surface recolor 2.5–3× (big-bcif 808 → ~280 ms); grid cells −2/3                     | M            | no     |
| 6   | 07     | P1         | Skip layer smoothing on visuals the overpaint/transparency layers do not touch        | [visual.ts:303](../../packages/graphics/src/repr/visual.ts#L303)                                                                                                           | time        | one-chain overpaint 700 ms → few ms on 416-unit surface; keep merged defines uniform | S (M merged) | no     |
| 7   | 10, 11 | P1, P1     | History checkpoint: memoize size per params; skip snapshot past hard limit (11/B2)    | [checkpoint.ts:81](../../packages/prism/viewer/src/sources/checkpoint.ts#L81), [history.ts:608](../../packages/prism/ui/src/managers/history.ts#L608)                      | time        | −35–40 % of every undoable op at N=100; 20–200 ms per agent call on 5–50 MB scenes   | S            | no     |
| 8   | 09     | P1         | Interactions: cylinder-impostor path instead of dashed meshes                         | [interactions-intra-unit-cylinder.ts:164](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L164)                             | time+memory | repr 1.1 s → ~0.15 s, buffers 422 → ~80 MB; size-theme caveat                        | M            | no     |
| 9   | 02     | P2, P3     | `GridLookup3D`: bucket-ordered Float32 coords; default 8 elements per cell            | [grid.ts:346](../../packages/core/src/math/geometry/lookup3d/grid.ts#L346)                                                                                                 | time        | queries ~2× (+ up to 1.4× from cells); bonds, contacts, ASA; +2.2 MB at 186k         | S            | no     |
| 10  | 10     | P4         | StructureInfo: three `rootsOfType` walks per `object.created` (plugin)                | [structure-info.ts:133](../../packages/plugin/core/src/behavior/dynamic/custom-props/structure-info.ts#L133)                                                               | time        | "color all" at N=100 3–9× faster (100 ms → linear); loads ~20 %                      | S            | no     |
| 11  | 06     | N1, P6, P5 | Interactions inter loop: query A's feature grid; 4 per cell; fast boundary for grids  | [contacts.ts:231](../../packages/model/src/props/computed/interactions/contacts.ts#L231)                                                                                   | time        | inter loop 352 → 138 ms (−44 %); feature boundaries 60 → 3 ms                        | S            | no     |
| 12  | 05     | P1         | `structureUnion`: per-unit concat + sort instead of repeated `SortedArray.union`      | [structure-set.ts:24](../../packages/model/src/model/structure/query/utils/structure-set.ts#L24)                                                                           | time+memory | 100k singletons in one unit 8.7 s → 18 ms; residue union 36 → 2.7 ms                 | S            | no     |
| 13  | 08     | P3         | Trim underfilled single-chunk builder buffers; MeshBuilder triangle estimate          | [mesh-builder.ts:47](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L47)                                                                                    | memory      | cartoon −43 % buffers, surfaces −20 % vertex buffers (CPU and GPU)                   | S            | no     |
| 14  | 01     | P1         | `ungzip` via native `DecompressionStream('gzip')`; drops JS CRC; fixes 01/B10, B11    | [zip.ts:231](../../packages/core/src/util/zip/zip.ts#L231)                                                                                                                 | time        | ~3.5×: 275 → ~90 ms per 100 MB; every `.gz` load                                     | S            | no     |
| 15  | 04     | P1         | Sequence alignment on flat `Int32Array` (2-row + 2-bit traceback later)               | [alignment.ts:32](../../packages/model/src/model/sequence/alignment/alignment.ts#L32)                                                                                      | time+memory | 2.3–5×, 5× less memory (3000 res: 204 ms/501 MB → 41 ms/101 MB); every PDB+SEQRES    | M            | no     |
| 16  | 05     | P6, P5     | Child structures: subset the parent's intra and inter-unit bond graphs                | [unit.ts:305](../../packages/model/src/model/structure/structure/unit.ts#L305), [structure.ts:240](../../packages/model/src/model/structure/structure/structure.ts#L240)   | time        | intra 90 → 1.6 ms; inter ~50–60 → <1 ms per large child                              | M            | no     |
| 17  | 07     | P2         | Segment visuals: groupCount 1 for 'object' granularity instead of whole-grid markers  | [segment.ts:234](../../packages/graphics/src/repr/volume/segment.ts#L234)                                                                                                  | memory      | −152 MB for 20 segments in 200³; ~2.7 GB for 100 in 300³                             | S            | no     |
| 18  | 04     | P2         | Multi-frame GRO: share the hierarchy, frames as coordinates                           | [gro.ts:134](../../packages/model/src/formats/structure/gro.ts#L134)                                                                                                       | memory      | 100–200 → 15–25 B/atom/frame (8–15×); mind 04/B7                                     | S            | no     |
| 19  | 05, 02 | P2, P6     | Boundary: 14 directions for structure and loci spheres; keep 98 for unit spheres      | [boundary.ts:15](../../packages/model/src/model/structure/structure/util/boundary.ts#L15)                                                                                  | time        | structure sphere 11–14 → 3–4 ms; lookup+boundary 74 → ~55 ms                         | S            | no     |
| 20  | 05     | P4         | `findPairBonds`: grid fetch of A's atoms in B's sphere instead of full scan           | [inter-compute.ts:105](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L105)                                                                | time+memory | inter bonds 55 → ~28 ms; 90 → <10 MB allocated per call                              | M            | no     |
| 21  | 09     | P6         | Memoize the particle-hierarchy color index across `Theme.create`                      | [particle-hierarchy.ts:134](../../packages/graphics/src/theme/color/particle-hierarchy.ts#L134)                                                                            | time+memory | ~340 ms + 40–80 MB per update at 10M particles → ~0                                  | S            | no     |
| 22  | 10, 03 | P5, P2     | Measurement sync: forEach / root children instead of 6 walks; then a cached pre-order | [measurement.ts:435](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L435), [immutable.ts:134](../../packages/core/src/state/tree/immutable.ts#L134) | time        | ~2 ms per op at 1.4k cells; walks 20–50× cheaper (~25 ms per change at 10k nodes)    | S / L        | no     |

Frozen, recorded only: 08/P7 cylinder impostor attributes ×6 per vertex (big-bcif 165 → ~35 MB, L, shader/schema in both
backends); 07/P7 merged objects keep a second CPU copy (policy change is repo-side, releasing member arrays is in
`gl/renderable/merged.ts`); 07/N2 root (merged values from member 0).

Next tier (medium, in the slice reports): 04/P3 IHM sharing, 04/P4 multi-model hierarchy compare, 04/P8 grid gradients,
05/P3 bundle symmetry sharing, 05/P7, 06/P2 feature cache key, 06/P4, 06/P7 ASA world grid, 07/P3–P5 marking, 08/P5,
08/P6, 08/P8, 09/P2, 10/P2, 10/P3, 01/P2, 03/P1.

## Cross-cutting themes

| Theme                                                 | Findings                                                                                                                | Pattern                                                                                                |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| Megamorphic helpers and closures in per-element loops | 02/P4 = 08/P4, 02/P5 = 08/P2, 01/P2, 02/P2, 07/P8                                                                       | `Vec3.fromArray`, `Tensor` `get`, `dataOffset`, `OrderedSet.getAt`: 3–5× in hot loops; inline or hoist |
| Per-element JS allocation                             | 06/P1, 04/P1, 04/P4, 04/P6, 04/P11, 06/P8, 02/P9, 04/P9, 08/P10, 09/P4, 11/P5                                           | arrays of arrays, tuples, template-string keys, regex per atom; use typed scratch buffers              |
| Quadratic merges and scans                            | 05/P1, 01/P3, 01/P7, 05/P4, 06/N1, 07/P3, 07/P10, 10/P4, 03/P6, 03/P9, 04/P12, 05/P9                                    | repeated union/scan from 0 per item; batch, index, or binary search                                    |
| Recompute where cache or update would do              | 08/P1, 07/P1, 09/P6, 05/P5, 05/P6, 06/P2, 06/P5, 10/P1, 10/P3, 10/P6, 10/P7, 03/P1, 03/P2, 01/P9                        | color-independent weights, untouched visuals, parent graphs, unchanged params                          |
| Symmetry or sharing lost                              | 05/P3, 05/B1, 06/P3, 06/P8, 06/P9, 04/P3, 04/P2                                                                         | per-unit/per-model/per-component work where a group or root could be shared                            |
| Caches that never evict, leaks                        | 01/B1 = 04/B1, 06/P2, 05/P11, 08/P9, 04/P8, 03/B6, 10/B2, 03/B7                                                         | module-level Maps keyed by data objects; missing dispose                                               |
| Oversized or duplicated buffers                       | 08/P3, 07/P2, 02/P7 = 07/P6, 08/P5, 08/P6, 08/P8, 09/P2, 04/P13, 01/B2, 08/P7, 07/P7                                    | estimates in the wrong unit, single chunk never trimmed, Float64 where Float32 suffices                |
| Grid and boundary over-precision                      | 02/P3, 06/P6, 02/P6, 05/P2, 06/P5                                                                                       | 32 elements per cell and 98 EPOS directions by default; only merged spheres need 98                    |
| Costs hidden from budgets                             | 10/P1, 11/P1, 11/B2                                                                                                     | the checkpoint size pass is outside `#take`'s timing                                                   |
| Same root cause found twice                           | 01/B1 = 04/B1, 02/P4 = 08/P4, 02/P5 = 08/P2, 02/P7 = 07/P6, 01/P5 ≈ 04/P5, 02/P9 ≈ 06/P8, 03/B10 ↔ 10/P5, 02/P3 ↔ 06/P6 | fix once, in the owning package                                                                        |

## Order of work

| Step | What                                                                                                                                                                                                                                                                                                            | Measure                                                                    |
| ---- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| 0    | Move the scratch pipeline bench into the repo (`tests/src/benchmark/pipeline.ts` or a `*.bench.test.ts` behind an env flag, as for MVS), with `mkbig` for `big-bcif`, and the Prism op bench (`bench-op100`). Add exactness probes: molsurf field/idField diff, MC checksum, ASA areas, grid hits, union output | baseline row per stage                                                     |
| 1    | Exact S-effort kernels: #1 MC, #2 addPrimitive, #3 ASA, #9 grid (P2 first, P3 changes edge order), #12 union, #14 ungzip; high bugs 03/B1, 04/B2, 06/B1                                                                                                                                                         | stage table: gaussian, interactions repr, asa, bonds; bit-identical probes |
| 2    | M rewrites with identical output: #4 molsurf, #5 color smoothing, #11 interactions loop, #20 findPairBonds; then #8 impostor interactions and #6 layer smoothing (needs the merged-define rule)                                                                                                                 | molsurf, theme-apply, interactions; a browser check for #6, #8             |
| 3    | Incremental and cached: #7 checkpoint size, #10 StructureInfo, #22 measurement sync and tree index, #16 child bonds, 10/P3 sequence per frame, #21 particle index                                                                                                                                               | Prism op bench at N=25/100/300; child-structure bench                      |
| 4    | Memory: #13 builder trim, #17 segments, #18 GRO, 04/P3 IHM, 04/P8 gradients, 08/P5–P6; hand frozen items (08/P7, 07/P7, 07/N2) to the WebGPU owners                                                                                                                                                             | retained MB, render-buffer bytes                                           |

## Open questions

- **Output changes**: which are acceptable? 02/P3 (contact and edge order), 02/P8 (density within fasterExp error,
  snapshots regenerate), 08/B1 (fillEdges starts filling holes), 01/B6 (palettes shift), 09/P2 (open dash ends), 05/P2
  (looser structure spheres), 06/N2 (weak H-bonds per copy).
- **Interactions impostors** (#8): impostors take size from the size theme, not `radius()`. Is a size iterator mapping
  each group to `min(sizeA, sizeB)` needed, or is the uniform default enough?
- **06/P3**: symmetry-equivalent unit pairs are exact only without inter-unit bonds and line of sight across units: a
  spec decision. 06/P9 waits on 06/N2.
- **03/B1**: carry `cell.state` into builder commits, or give `updateCellState` a tree identity the builder detects?
  Both change tree.md §2/§4.
- **03/B5**: hash versions over tags and `dependsOn` change the MVS `VirtualStateBuilder` too (load.md §3).
- **04/B2**: `ModelSymmetry` became `asDynamic` in 6dfd2b734; per-model dynamic data, or back to static?
- **10/P8**: trajectory seeks: `record: false` or coalesce? The agent may want the frame in the record.
- **08/P8**: builder `initialCount` in items changes a public API (MVS calls `LinesBuilder.create(1024, 1024)`).
- **07/P7**: unmerging big meshes needs a GPU frame-time measurement first.
- **11/B1**: drop `key?` from `For`'s props and migrate ~40 sites to `by=`?
- **Bench home**: one repo bench for the pipeline and one for Prism ops, or per-package `*.bench.test.ts`?
