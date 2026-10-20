# Profile: where pipeline time and memory go

Measured 2026-10-19 on `v6-opts-261019` (`90d717f49`), Node 25.6, Apple M4 (10 cores), headless, no GPU. Read-only
audit: nothing in the repo was changed to measure this.

## Fixtures

| Id           | File                                                                      |    Size | Atoms / elements | Residues | Units (groups) | Kind                                           |
| ------------ | ------------------------------------------------------------------------- | ------: | ---------------: | -------: | -------------: | ---------------------------------------------- |
| `7qpd`       | `data/examples/7qpd.fw2.cif`                                              | 1.8 MB  |           11,521 |    1,554 |          7 (7) | atomic mmCIF text (no `chem_comp`)             |
| `3k7f`       | `packages/prism/viewer/src/_stories/dynamic-pdb/fixtures/3k7f-qfit.bcif`  | 571 KB  |           13,904 |    2,808 |          8 (8) | BCIF, qFit alt locs, P3₂21                     |
| `8zzh`       | `data/examples/ihm-coarse-8zzh.cif`                                       | 298 KB  |  1,764 sph + 1,764 gauss |   – |          4 (4) | IHM coarse only                                |
| `1tup-cell`  | `…/sequence/molstar/_test/fixtures/structures/1tup.bcif`, sym 0..0         | 342 KB  |           23,312 |    4,056 |        36 (9)  | symmetry copies, shared groups                 |
| `1tup-2x2x2` | same, `buildSymmetryRange` 000..111                                       | 342 KB  |          186,496 |  ~32,400 |      288 (9)   | 32 copies, 9 symmetry groups (instancing)      |
| `big-cif` / `big-bcif` | `1tup-2x2x2` written by `to_mmCIF` as one model (scratch `fixtures/`) | 20.5 MB / 3.4 MB | 186,496 | ~32,400 | 416 (416) | large entry, nothing shared between units |

`big-*` is the realistic "large deposited entry" case (a cryo-EM-sized model): every unit is its own group, so no work
is shared. `1tup-2x2x2` is the same atoms with symmetry groups (9 representatives, instanced).

## How to rerun

The bench now lives in the repo (`tests/src/pipeline/`, `node scripts/pipeline-bench.mjs`; `flat-32` is `big-bcif`),
see [13-large-scale](13-large-scale.md#how-to-rerun). The scratch scripts below are what produced this report.

Scripts in `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/`
(`bench.ts` pipeline, `build.mjs` esbuild bundle with `conditions: ['molstar-src']`, `analyze.mjs` per-stage
self/total from `.cpuprofile` mapped back to TS by the bundle's source map, `mkbig.ts` big fixtures, `snapdiff.py` /
`retainers.py` heap-snapshot diffs). `node_modules` there is a symlink to `tests/node_modules`.

```sh
cd <scratch>/perf-audit
node build.mjs bench.ts bench.mjs                       # bundles the TS sources (tsx is not installed; esbuild is)
node build.mjs mkbig.ts mkbig.mjs && node mkbig.mjs      # once: fixtures/1tup-2x2x2.{cif,bcif}
node --expose-gc bench.mjs --fixture 7qpd --runs 5 > out-7qpd.json && python3 show.py out-7qpd.json
node --expose-gc bench.mjs --fixture 7qpd --profile prof/7qpd && node analyze.mjs prof/7qpd 8
node --expose-gc bench.mjs --fixture 7qpd --runs 5 --snapshots snap/7qpd   # heap snapshots after run 2 and 5
# flags: --only a,b  --skip a,b  --mesh (no impostors)  ; big fixtures need --max-old-space-size=12000
```

Worth moving into the repo as a bench (e.g. `tests/src/benchmark/pipeline.ts`, or a `*.bench.test.ts` behind an env
flag like `MVS_BENCH`), so slice fixes can be measured before/after.

Method: each run does the whole pipeline from the file bytes; stages in order. Time = median of N runs (5 for small
fixtures, 2 for `1tup-2x2x2`/`big-bcif`; the cold first run is in the JSON). Allocated = V8 `total_allocated_bytes`
delta (JS heap only; typed-array backing stores are not counted). Retained = (`heapUsed` + `arrayBuffers`) after
`gc()` minus before, with all stage outputs kept alive until the run ends. Profiles: `node:inspector` `Profiler`
started/stopped around each stage (100 µs sampling) after one warm-up run; inspector overhead removed.

Representations: `provider.factory({ graphics: stub, …registries })` + `ensureCustomProperties` + `createOrUpdate`
with defaults and `tryUseGpu: false`. The stub reports `fragDepth` so sphere/cylinder **impostors** are used as in a
browser; its textures ignore uploads (the CPU work before the upload is measured). Without a graphics context
everything falls back to meshes (see Observations). Theme "eval" = `color(location)` for every element; theme "apply" =
`setTheme` + `createOrUpdate({})` on an existing repr with the size theme unchanged (no geometry rebuild).

## Stage table

ms (median) / MB allocated (JS heap) / MB retained (heap + array buffers)

| Stage | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 | big-bcif |
| --- | ---: | ---: | ---: | ---: | ---: |
| `parse` | 12 / 1 / 5 | 4 / 1 / 1 | 5 / 0 / 1 | 6 / 8 / 2 | 5 / 3 / 4 |
| `model` | 7 / 4 / 1 | 7 / 5 / 3 | 6 / 9 / 2 | 8 / 13 / 2 | 91 / 123 / 58 |
| `structure` (+ symmetry) | 1 / 1 / 0 | 0.8 / 0 / 0 | 0.5 / 0 / 0 | 2 / 6 / 0 | 9 / 4 / 2 |
| `lookup+boundary` | 4 / 1 / 0 | 2 / 1 / 0 | 1 / 1 / 0 | 16 / 20 / 0 | 74 / 47 / 8 |
| `bonds-intra` | 18 / 5 / 0 | 26 / 5 / 0 | 0.1 / 0 / 0 | 12 / 9 / 1 | 242 / 38 / 7 |
| `bonds-inter` | 3 / 2 / 0 | 1 / 0 / 0 | 0.4 / 1 / 0 | 60 / 99 / 2 | 51 / 64 / 0 |
| `secondary-structure` | 0.4 / 0 / 0 | 3 / 2 / 0 | 0.2 / 1 / 0 | 0.6 / 0 / 0 | 2 / 2 / 0 |
| `valence-model` | 5 / 8 / -1 | 7 / 9 / -3 | 0.1 / 1 / 0 | 48 / 92 / -12 | 48 / 93 / -30 |
| `interactions` | 36 / 28 / -4 | 98 / 80 / -1 | 0.8 / 1 / 0 | 478 / 288 / 20 | 890 / 705 / 11 |
| `asa` | 113 / **778** / 0 | 185 / **868** / 0 | **throws** | 1921 / **8283** / 2 | 1825 / **7546** / 2 |
| `repr:cartoon` | 20 / 12 / 24 | 21 / 27 / 17 | 18 / 25 / 14 | 14 / 23 / 4 | 184 / 370 / 126 |
| `repr:ball-and-stick` | 25 / 15 / 12 | 27 / 14 / 12 | 3 / 2 / 0 | 22 / 22 / 10 | 456 / 396 / 389 |
| `repr:spacefill` | 4 / 3 / 1 | 4 / 3 / 1 | 2 / 2 / 0 | 4 / 5 / 1 | 91 / 112 / 42 |
| `repr:molecular-surface` | **1691** / 85 / 27 | **2210** / 91 / 26 | 742 / 132 / 20 | 424 / 120 / 3 | **12917** / 1577 / 171 |
| `repr:gaussian-surface` (CPU) | 612 / 43 / 34 | 556 / 36 / 31 | 786 / 51 / 19 | 41 / 40 / 4 | 806 / 222 / 195 |
| `repr:interactions` | 21 / 4 / 30 | 94 / 9 / **121** | 1 / 1 / 0 | 430 / 974 / **227** | 1101 / 1847 / **643** |
| `theme-eval:element-symbol` | 0.7 / 0 / 0 | 0.7 / 0 / 0 | 0.3 / 0 / 0 | 6 / 1 / 0 | 13 / 0 / 0 |
| `theme-eval:chain-id` | 0.4 / 0 / 0 | 0.4 / 0 / 0 | 0.2 / 0 / 0 | 5 / 1 / 0 | 5 / 0 / 0 |
| `theme-eval:secondary-structure` | 1 / 0 / 0 | 1 / 0 / 0 | 0.7 / 0 / 0 | 6 / 1 / 0 | 16 / 2 / 0 |
| `theme-eval:residue-name` | 1 / 1 / 0 | 0.9 / 0 / 0 | 0.5 / 1 / 0 | 4 / 0 / 0 | 12 / 0 / 0 |
| `theme-eval:accessible-surface-area` | 2 / 1 / 0 | 2 / 1 / 0 | 0.3 / 1 / 0 | 17 / 4 / 0 | 26 / 4 / 0 |
| `theme-apply:spacefill:chain-id` | 1 / 0 / 0 | 1 / 0 / 0 | 0.7 / 0 / 0 | 1 / 0 / 0 | 21 / 5 / 0 |
| `theme-apply:cartoon:secondary-structure` | 0.7 / 0 / 0 | 0.5 / 0 / 0 | 1 / 0 / 0 | 0.8 / 0 / 0 | 16 / 6 / 0 |
| `theme-apply:cartoon:accessible-surface-area` | 0.9 / 0 / 0 | 0.6 / 0 / 0 | **throws** | 0.6 / 0 / 0 | 16 / 6 / 0 |
| `theme-apply:gaussian-surface:element-symbol` | 252 / 18 / 0 | 247 / 15 / 0 | 181 / 22 / 0 | 30 / 15 / 0 | 808 / 160 / -74 |
| `theme-apply:molecular-surface:residue-name` | 213 / 16 / 0 | 208 / 14 / 0 | 191 / 21 / 0 | 19 / 1 / 0 | 713 / 122 / -64 |

`big-cif` (text, 20.5 MB): `parse` 90 ms (52.8 MB retained token buffers), `model` 88 ms / 77 MB alloc.

Render buffers (sum of typed arrays in the render objects' values; retained above can be larger, see merged):

| Repr | 7qpd | 3k7f | big-bcif | Largest buffers |
| --- | ---: | ---: | ---: | --- |
| cartoon | 23.4 MB | 16.1 MB | 36.2 MB (retained 88) | `elements` (Uint32) 57 %, `aPosition`, `aNormal` |
| ball-and-stick (impostor) | 10.9 MB | 11.3 MB | 165 MB (retained 333) | `aMapping`, `aStart`, `aEnd` 31.5 MB each |
| ball-and-stick (`--mesh`, no impostor) | **230 MB** | – | – | `elements` 129, `aPosition`/`aNormal` 43 each |
| spacefill (impostor) | 0.5 MB | 0.5 MB | 4.5 MB (retained 12) | `tPositionGroup` |
| spacefill (`--mesh`) | 114 MB | – | – | |
| molecular-surface | 27.0 MB | 25.5 MB | 62.8 MB (retained 141) | |
| gaussian-surface | 33.6 MB | 30.2 MB | 72.8 MB (retained 165) | |
| interactions | 29.1 MB | **120.6 MB** | **421.6 MB** (retained 623) | mesh: `aPosition`/`aNormal` 143 MB each |

## Hotspots by stage

Self time unless marked (total). Percent of the stage's sampled time. Fixture in brackets.

| Stage | Function (file:line) | Share | Note |
| --- | --- | ---: | --- |
| parse (text) [big-cif 99 ms] | `eatValue` `io/src/reader/cif/text/parser.ts:65` | 33 % | + `moveNextInternal` :419 21 %, `readLoopChunk` :626 19 % |
| parse (BCIF) | msgpack `decode.ts:121/103` | – | 5 ms for 186k atoms: columns decode lazily, cost moves to `model` |
| model [big-bcif 90 ms] | `arrayView` `core/src/data/db/column.ts:485` | 15 % | + `addMissingResiduesFromEntityPolySeq` `model/src/formats/structure/basic/properties.ts:64` 10 %, `sortAtomSite` `basic/sort.ts:18` (total 27 %), `_makeBuckets` `core/src/data/util/buckets.ts:20` |
| lookup+boundary [big-bcif] | `computeSphereExtrema` → `Vec3.dot` `core/src/math/geometry/boundary-helper.ts:43` | 58 % total | EPOS boundary per unit |
| bonds-intra [big-bcif 242 ms] | `GridLookup3D.query` `core/src/math/geometry/lookup3d/grid.ts:305` | 48 % | `findBonds` `model/src/model/structure/structure/unit/bonds/intra-compute.ts:151` (total 92 %) |
| bonds-inter [1tup-2x2x2 60 ms] | `query` grid.ts:305 27 %, `findPairBonds` `unit/bonds/inter-compute.ts:50` 16 % | | linear in units (36 → 288 units: 6 → 60 ms) |
| valence-model [1tup-2x2x2 48 ms] | `calculateHydrogensCharge` `model/src/props/computed/chemistry/valence-model.ts:104` | 92 % total | `isConjugated` :41 46 % total, `getEdgeIndices` `core/src/math/graph/inter-unit-graph.ts:57` 21 % |
| interactions [big-bcif 890 ms] | `query` grid.ts:305 | 30 % | `addStructureContacts` `model/src/props/computed/interactions/contacts.ts:197` 51 % total; intra-unit contacts 34 % (`interactions.ts:424`) |
| asa [1tup-2x2x2 1.9 s, 8.3 GB alloc] | `computeRange` `model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts:33` | 54 % | per atom `neighbors = []` of 5-element arrays (:64, `// TODO reuse`) + `sort` closure :88; `lookup3d.find` 40 % total |
| repr:cartoon [big-bcif 184 ms] | `createPolymerTraceMesh` `graphics/src/repr/structure/visual/polymer-trace-mesh.ts:62` | 50 % total | `addTube` `geo/geometry/mesh/builder/tube.ts:59`; merge `createMergedValues` `gl/renderable/merged.ts:308` 16 % |
| repr:ball-and-stick [big-bcif 456 ms] | `createLinkCylinderImpostors` `graphics/src/repr/structure/visual/util/link.ts:356` | 40 % total | `CylindersBuilder.add` `geo/geometry/cylinders/cylinders-builder.ts:103` 17 %, `ChunkedArray.allocateNext` |
| repr:molecular-surface [big-bcif 12.9 s] | `obscured` `core/src/math/geometry/molecular-surface.ts:87` | 30 % | + `singleAtomObscures` :111 23 %, `projectTorus` :226 23 %, `OrderedSet.getAt` 7 % |
| repr:molecular-surface [7qpd 1.7 s] | `projectTorus` :226 28 % | | `Tensor` `get` `core/src/math/linear-algebra/tensor.ts:292` 15 % (marching cubes), `OrderedSet.indexOf` → `binarySearchRange` `core/src/data/int/impl/sorted-array.ts:130` 6 % (per voxel write, :285) |
| repr:gaussian-surface [7qpd 612 ms] | `Tensor` `get` tensor.ts:292 (called from `MarchingCubesState.get` `graphics/src/geo/util/marching-cubes/algorithm.ts:246`) | **45 %** | + `processCell` :250 14 %, `accumulateRange` `core/src/math/geometry/gaussian-density/cpu.ts:67` 21 % |
| repr:interactions [big-bcif 1.1 s, 1.8 GB alloc] | `Vec3.fromArray` `core/src/math/linear-algebra/3d/vec3.ts:90` | **41 %** | inside `MeshBuilder.addPrimitive` `geo/geometry/mesh/mesh-builder.ts:122` (77 % total) ← `addFixedCountDashedCylinder` `geo/geometry/mesh/builder/cylinder.ts:120` |
| theme-apply on surfaces [big-bcif 808 ms] | `calcMeshColorSmoothing` `graphics/src/geo/geometry/mesh/color-smoothing.ts:43` | 77 % | + `getIndex` :86 14 % |
| theme-apply group-granular | `createGroupColor` `graphics/src/geo/geometry/color-data.ts:268` | 57 % total | 21 ms for 186k atoms; fine |
| theme-eval | `color` `graphics/src/theme/color/chain-id.ts:97` etc. | | ≤ 26 ms for 186k elements; `accessible-surface-area` is the slowest (2× others) |

## Observations

1. **Molecular surface dominates everything**: 1.7–2.2 s for 11–14k atoms, 12.9 s for 186k atoms in 416 units.
   `obscured` scans the whole neighbour list (via `OrderedSet.getAt`) for every probe position of every torus;
   `projectTorus` resolves the atom index by `OrderedSet.indexOf` (binary search) per voxel write, though `a` and `b`
   are fixed per torus. Main-thread task (yields between chunks).
2. **Marching cubes walks the field against its memory layout**: density fields are `Tensor.Space(dim, [0, 1, 2])`
   (axis order slow→fast: z fastest), `MarchingCubesComputation.slice` loops `i` (x) innermost
   (`algorithm.ts:99`), so every `get` strides `dimY·dimZ` floats; each cell does 8 closure `get`s (each corner read 8×).
   45 % of the CPU gaussian surface, 15 % of the molecular surface. `geo/util/marching-cubes` is not under `gl/`.
3. **Surface color change recomputes the full color-smoothing splat**: 250 ms per theme change on a 12k-atom surface,
   800 ms on big-bcif; 216 cells/vertex with `sqrt` + `getIndex` (2 floor + modulo). The cell weights depend only on
   geometry, not on colors.
4. **ASA allocates ~40 KB per atom** (8.3 GB for 186k atoms; 778 MB for 11.5k): array-of-arrays neighbour list and
   a comparator sort per atom (`area.ts:64,88`), sphere points as `[x,y,z]` arrays destructured per point.
5. **Interactions representation builds dashed mesh cylinders** (no impostor path, unlike bonds which use
   `createLinkCylinderImpostors` with `addFixedCountDashes`): 121 MB of buffers for 3k7f, 422 MB (623 MB retained) for
   big-bcif. `graphics/src/props/computed/representations/interactions-{intra,inter}-unit-cylinder.ts`.
6. **Merged render objects keep two CPU copies**: when `renderMerged: 'auto'` merges many unit visuals (big-bcif),
   `createMergedValues` (`gl/renderable/merged.ts:308`) copies into own cells while members keep theirs; retained is
   ~2.2× the render buffers (cartoon 36 → 88 MB, ball-and-stick 165 → 333 MB, gaussian 73 → 165 MB). Frozen area.
7. **Cylinder impostors repeat per-cylinder data per vertex**: ball-and-stick `aMapping`/`aStart`/`aEnd` are 31.5 MB
   each for 180k bonds (gl schema, frozen).
8. **No-impostor fallback is 20–200× bigger**: ball-and-stick 230 MB, spacefill 114 MB for 7qpd as meshes (what
   headless/Node capture and devices without `fragDepth` get).
9. **Symmetry sharing works**: `1tup-2x2x2` (9 groups) does surfaces in 3–7 % of `big-bcif`'s time for the same
   atoms. Non-shared work that scales with units: `bonds-inter`, `interactions` (inter-unit), `asa` (computed for every
   unit), `lookup+boundary`. No super-linear scaling seen from 36 → 288 units (×8 atoms: bonds-inter ×9.4,
   interactions ×6, asa ×7.6).
10. **Parse is cheap; BCIF defers decoding to `model`** (5 ms parse, 91 ms model for 186k atoms; text 90 + 88 ms).
11. **Theme evaluation itself is cheap** (`graphics/src/theme`): ≤ 26 ms per 186k elements; group-granular apply
    ≤ 21 ms. Cost of a theme change is in the surface color smoothing (3), not in the themes.

## Bugs found

| Bug | Where | Failure scenario |
| --- | --- | --- |
| **Leak**: `getUniqueComponentNames = memoize1(...)` keyed by `BasicData` in a module-level `Map` that never evicts | `model/src/formats/structure/basic/properties.ts:177` (`core/src/util/memoize.ts:29`) | Any file without `chem_comp` (and the saccharide path, :169): the whole `BasicData` stays reachable, and through its CIF columns the full source text. Heap grows 4.1 MB per load of `7qpd.fw2.cif` (5 runs: 12.9 → 29.5 MB), ~1 MB per load of `8zzh`; heap snapshot retainer: `getUniqueComponentNames` closure → `cache` Map → `atom_site.label_asym_id` → `str10` → the 1.8 MB string. Use a `WeakMap` or compute in place |
| ASA throws on any structure with coarse units | `model/src/props/computed/accessible-surface-area.ts:62` (`isApplicable: () => true`) → `shrake-rupley/radii.ts:39` (`StructureProperties.residue.key` on a coarse unit) | `AccessibleSurfaceAreaProvider.attach` on `8zzh` (or a mixed atomic + coarse IHM entry): "Property only available for atomic models." The `accessible-surface-area` color theme (`graphics/src/props/computed/themes/accessible-surface-area.ts:90`, applicable to any structure) fails the same way |
