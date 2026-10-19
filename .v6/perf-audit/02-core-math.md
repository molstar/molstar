# Slice 02: core math

Scope: [`packages/core/src/math`](../../packages/core/src/math) (16k LOC): geometry (molecular surface, gaussian
density, `GridLookup3D`, boundary helpers, primitives), linear algebra (`Vec3`, `Tensor`, RMSD), graphs. Two fixes land
in their graphics callers (mesh builder, marching cubes, the molecular-surface visual); they are listed here because
the cost shows up in core functions. MVS loading is in [mvs-optimizations](../mvs-optimizations/research/), not here.
Fixtures and method: [00-profile.md](00-profile.md). No finding touches Canvas3D or `gl`. Measured 2026-10-19, Node
25.6, Apple M4.

## Hot paths at 186k atoms

| Path | Where it is paid | Cost (big-bcif unless noted) |
| --- | --- | --- |
| `calcMolecularSurface` | `repr:molecular-surface`; auto quality is `low` (res 1.3, 24 probes) above 100k atoms, which is the extended mode (`probeRadius < 2·res`), as is `medium` (res 0.8) | 12.3–12.9 s of 13.4 s; `obscured` + `singleAtomObscures` 53 % self, `projectTorus` 23 %, `getAt` 7 % |
| `GridLookup3D` `query` | bonds-intra/inter (`unit.lookup3d`), ASA and interactions (`StructureLookup3D`), feature grids, the molsurf neighbour search | 25–48 % self of those stages; one query per atom: 84–155 ms (r 2–4) |
| `Vec3.fromArray` (megamorphic) | `MeshBuilder.addPrimitive` → dashed cylinders | 420 ms self of 1.02 s `repr:interactions` (41 %); 172 of 410 ms on 1tup-2x2x2 |
| `Tensor.Space.get` closure | marching cubes (`MarchingCubesState.get`, 8 per cell) | 45 % of `repr:gaussian-surface` on 7qpd (303 of 677 ms) |
| `accumulateRange` | `GaussianDensityCPU` | 21 % of `repr:gaussian-surface` on 7qpd |
| `getBoundary` (98 directions) | every unit, feature set, loci, geometry under 10k points | 33–48 ms over 416 units |
| `InterUnitGraph.getEdgeIndices` | valence model, bond iteration | 12 ms of 59 ms valence-model |

Already efficient: the grid build (resolved indices, radix-sorted sparse cells, `MaxVolume` guard), symmetry-operator
`ArrayMapping`, histogram, `CentroidHelper`, `IntAdjacencyGraph` builders, instance-grid.

## Performance findings (confirmed and plausible, ranked)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | Molsurf occlusion rescans the full neighbour list through `OrderedSet.getAt` per probe and voxel | [molecular-surface.ts:87](../../packages/core/src/math/geometry/molecular-surface.ts#L87) | time | high | 1.4–1.7×: big-bcif 12.3–12.7 → 7.4–7.8 s, 7qpd 1.24 → 0.87 s (res 0.5), 0.76 → 0.44 s (res 1.3); bit-identical | M | confirmed | no |
| P2 | Grid query dereferences `OrderedSet` and scattered x/y/z per candidate | [grid.ts:346](../../packages/core/src/math/geometry/lookup3d/grid.ts#L346) | time | high | query ~2×: r 4 155 → 74 ms, r 6.4 262 → 141 ms at 186k; +12 B per element retained | S | confirmed | no |
| P4 | `Vec3.fromArray` megamorphic in `addPrimitive` and other per-vertex loops | [vec3.ts:90](../../packages/core/src/math/linear-algebra/3d/vec3.ts#L90), [mesh-builder.ts:122](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L122) | time | high | `repr:interactions` −30–40 % (big-bcif ~1.0 → ~0.65 s); 5.1× on the loop | S | confirmed | no |
| P5 | Marching cubes reads the field through a per-voxel `Tensor.Space.get` closure, x innermost | [tensor.ts:292](../../packages/core/src/math/linear-algebra/tensor.ts#L292), [algorithm.ts:246](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L246) | time | high | classification 3.1× with offsets alone, 4.1× with z innermost; `repr:gaussian-surface` −30 % (7qpd ~677 → ~470 ms) | S (offsets) / M (reorder) | confirmed | no |
| P3 | Default grid of 32 elements per cell is too coarse for 2–6 Å queries | [grid.ts:243](../../packages/core/src/math/geometry/lookup3d/grid.ts#L243) | time | medium | −15–40 % of query time on top of P2 (r 2: 39 → 28 ms; r 6.4: 141 → 134 ms) | S | confirmed | no |
| P6 | `getBoundary` uses 49 directions for every set under 10k points | [boundary.ts:22](../../packages/core/src/math/geometry/boundary.ts#L22) | time | low | 33–48 → 12–19 ms at 186k with the coarse helper, but derived spheres stop enclosing atoms; safe variant gains less | S | plausible | no |
| P8 | Gaussian density: `fasterExp` per voxel; the kernel is separable | [cpu.ts:101](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L101) | time | low | `accumulateRange` −20–25 % (7qpd res 0.5: 155–173 → 122–137 ms) | S | confirmed | no |
| P7 | Per-unit molsurf radius array sized `OrderedSet.end(indices)` | [molecular-surface.ts:48](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L48) | memory | low | 214 MB of virtual allocation vs 0.7 MB compact on big-bcif; time negligible | S | confirmed | no |
| P9 | `getEdgeIndices` allocates `[]` per miss; edges as objects in nested Maps | [inter-unit-graph.ts:58](../../packages/core/src/math/graph/inter-unit-graph.ts#L58) | time+memory | low | valence model and bond iteration −10–20 % (~−10 ms at 186k) | S | confirmed | no |

### P1: molecular-surface occlusion

- Mechanism: `obscured()` loops over the whole `lookup3d.find` result of the current atom (60–100 neighbours), with
  `OrderedSet.getAt` and scattered `px/py/pz/radius` loads per element. `projectTorus` calls it for every probe
  position of every pair (24–36 × ~50 per atom). `projectPointsRange` calls it for every voxel of the atom cube
  before checking whether the voxel would change; in extended mode that is every cube voxel. `projectTorus` also
  binary-searches `id[OrderedSet.indexOf(indices, a|b)]` per written voxel although `a`, `b` are fixed.
- Exactness: only the `-1` / not `-1` outcome of `obscured()` is used. Any atom that obscures a point on atom j's
  sphere is within `r_j + r_i` of j, so it is in j's neighbour list; the `lastClip` cache never changes the outcome.
  Probe points lie on a circle (centre `mid`, radius `rInt`), so an occluder needs `|c − mid| < rInt + r_c`: per-torus
  pruning is exact (eps 1e-4 only guards rounding).
- Evidence (`bench-molsurf`, rerun): 7qpd res 0.5 / 36: 1241–1324 → 870–966 ms; res 1.3 / 24: 759–768 → 442–470 ms;
  big-bcif (256 units > 50 atoms) res 1.3 / 24: 12.3–12.7 → 7.4–7.8 s. `fieldDiff 0`, `idDiff 0` over 16.7M / 1.0M /
  18.1M voxels.
- Fix (`molsurf-opt.ts`): after each `find`, `loadNeighbours()` into growable typed arrays (resolved index, x, y, z, r,
  r²); occlusion over those with an integer `lastClip` slot; per torus a pruned candidate list (skip a, b); hoist
  `idA`/`idB`; in `projectPointsRange` compute `dd` first and test occlusion only when `dd < data[idx]` (inside) or
  `dd > data[idx]` (extended); reuse resolved indices in `projectToriiRange`.
- Risk: none seen. Add a unit test comparing field and idField with today's output on a small structure.

### P2: grid query memory layout

- Mechanism: each candidate costs `OrderedSet.getAt(indices, bucketArray[i])` (Interval/SortedArray dispatch) plus
  random loads from `px/py/pz[/radius]`. `_build` already resolves the indices once.
- Evidence (`bench-grid`, rerun, one find per atom, same hits in the same order): 186k / 416 units: r 2 84–88 → 39–41
  ms, r 4 154–156 → 74–75 ms, r 6.4 262–275 → 141–142 ms, with radii r 3.2 149–155 → 92–97 ms; 7qpd r 4 12 → 5.5 ms.
  Build +0.1–1 ms. Variant without copies (resolve `sorted[...]` / `start + ...` in `query`, `v/grid-res.ts`): r 4
  110 ms, r 6.4 205 ms, i.e. a third of the gain at zero memory.
- Fix: in `_build`, bucket-ordered `bx/by/bz` (`Float32Array` when the source is, which conformations are, so values
  stay exact) and `br` when radii are given; `query` scans `[offset, end)` linearly, keeps `bucketArray[i]` for
  `Result.add`. `approxQueryNearest` and `queryNearest` can read the same arrays.
- Risk: +12 B per element (+4 with radii) retained per cached lookup (~2.2 MB at 186k). Copies are taken at build
  time: in-place coordinate edits after the build are not seen. `Unit.remapModel` drops `lookup3d` when the
  conformation changes; other builders (instance-grid, blob, tunnels) build and query in one go.

### P4: megamorphic `Vec3.fromArray`

- Mechanism: 134 call sites pass `Float32Array`, `Float64Array`, `number[]`, int arrays; inlined copies share the
  callee's feedback, so the three keyed loads go megamorphic in hot loops: `addPrimitive`/`addPrimitiveFlipped`
  (dashed cylinders, mesh spheres), `addTriangleFan(WithNormal)`, `addCage`, frenet frames, bounding-sphere code.
- Evidence: profile above (41 % self of `repr:interactions`); `bench-fromarray` (rerun, 20M vertices from a
  `Float32Array`): clean feedback helper 54–56 ms = inline 54 ms; after 8 array kinds helper 276–280 ms, inline 54 ms
  (5.1×).
- Fix: read the arrays inline in the mesh-builder hot loops (`x = va[i], y = va[i+1], z = va[i+2]`, transform without
  the `Vec3` round trip). Optionally a documented hot-loop helper in core (e.g. `Mat4.transformPointArray(out, oo,
  src, so, m)`). Leave `Vec3.fromArray`'s signature; skip `gl/` (frozen).
- Risk: none (same arithmetic). The gain depends on which array kinds a session has fed the helper; measure in the
  app.

### P5: marching-cubes field access

- Mechanism: `Space.get` is a closure per `Space`, called 8× per cell through `MarchingCubesState.get`, polymorphic
  across spaces. Molsurf and gaussian fields use slow-to-fast `[0,1,2]` (z fastest) while `slice()` runs x innermost,
  so each read strides `dimY·dimZ`.
- Evidence (`bench-tensor` + `v/bench-tensor2`, 160×170×180, 8-corner classification of every cell, same result):
  closure, x innermost 351–384 ms; direct offsets, x innermost 115 ms; closure, z innermost 325 ms; direct offsets, z
  innermost 85–88 ms. Most of the gain is the closure, not the order.
- Fix, step 1 (S, no output change): in `MarchingCubesState`, take `o = space.dataOffset(i, j, k)` once per cell and
  the strides `sx = dataOffset(1,0,0) − dataOffset(0,0,0)` etc. (every 3D layout is linear), read `data[o]`,
  `data[o+sx]`, …; or carry the offset along `i`. No core API change needed; an explicit `strides` on `Space` is
  optional. Step 2 (M): slice along the slowest axis with the fastest innermost and rework the two-layer
  `verticesOnEdges` cache.
- Risk: step 1 none; step 2 changes vertex order (same geometry) and the edge cache. The WebGPU work touched
  marching-cubes tests (`f74da126d`): coordinate, though `geo/util` is not frozen.

### P3: default grid cell count

- Mechanism: with no cell size, `build` targets 32 elements per cell averaged over the box; protein interiors are
  denser, so a 2–4 Å query scans cells of mostly non-hits. Users of the default: `unit.lookup3d`
  ([unit.ts:405/651](../../packages/model/src/model/structure/structure/unit.ts#L405)), feature grids, the DSSP trace
  lookup, `StructureLookup3D.unitLookup`.
- Evidence (`bench-grid`, rerun, hits identical): count 8, original query: r 2 84–88 → 50 ms, r 4 155 → 106–114 ms, r
  6.4 262–275 → 215–220 ms; on top of P2: 39 → 28, 75 → 63, 141 → 134 ms. Build unchanged (2–3 ms). Count 4 is no
  better than 8 after P2.
- Fix: default 8, or let unit lookups pass ~4 Å cells. Slice 06 P6 already sets 4 per cell for feature grids.
- Risk: result order changes with the grid, so intra-bond edge order and contact order change (equivalent graphs);
  snapshot tests may need updating.

### P6: 98-direction boundary

- Mechanism: under 10k points `getBoundary` uses the 49-direction helper (49 dots + compare/copy per point); the
  14-direction helper with `refineExtrema` gives the same box at 3× lower cost.
- Correction: the auditor's "the radius always encloses the points" holds for the unit sphere itself, not for spheres
  built from it. `includeSphere`/`radiusSphere` use only the `extrema` points, and the structure boundary, loci,
  camera and scene spheres are built that way. Measured (`v/boundary-enclose`): atoms (centre + VdW) outside the
  structure sphere by at most 0.008 Å (7qpd) / 0.017 Å (big-bcif) with 98 directions, 0.88 / 0.12 Å with the coarse
  helper (26 directions: 0.008 / 0.12 Å); per unit from its own extrema 0.045 → 0.86 Å. The coarse 98-slot extrema
  hold only 14 distinct points.
- Evidence (`bench-boundary`, rerun, 416 units): fine 33–48 ms, coarse 12–19 ms; box max diff 0; radius ratio 0.998
  (worst 1.025).
- Fix: keep the 98 directions and make the pass cheaper: flatten `dir`/`minDist`/`maxDist` into `Float64Array`s,
  track the extremum index per direction and build the extrema once at the end (prototype `v/boundary-flat.ts`: the
  pass alone 20–25 ms); dedupe the refined extrema in coarse mode (keep ≥ 14 for `Box3D.fromSphere3D` and
  `Sphere3D.expand`). Use the coarse helper (or `getFastBoundary`) only where the sphere is never merged, e.g. feature
  grids (slice 06 P5).
- Risk: lowering the threshold makes camera fit and the scene sphere miss up to ~0.9 Å; large units (> 10k) already
  do today.

### Low

- **P8** gaussian density: `exp(k(dx²+dy²+dz²)) = ex·ey·ez`; per atom 1D tables, skip rows with `dxySq > r2sq`, hoist
  `id ? id[i] : i`. Rerun (`bench-gauss`, 7qpd res 0.5): 155–173 → 122–137 ms, idField identical, max |Δ| 0.057 of
  2.5 (the `fasterExp` error; `Math.exp` is the more accurate value). Isosurfaces move by a hair: regenerate
  snapshots.
- **P7** molsurf radii: `new Float32Array(OrderedSet.end(indices))` per unit; on a one-model 416-unit entry the sum is
  214 MB (est-radius) against 0.7 MB. The arrays are calloc'd and only the unit's own tail pages are written, so
  allocation time is ~2 ms and RSS barely moves; it counts as V8 external memory (GC pressure) and can fail on very
  large entries. `getStructurePositionDataAndMaxRadius` is compact (`OrderedSet.ofRange(0, n)`), not affected. Fix:
  compact x/y/z/r over `t` with `indices = Interval(0, n)` in the graphics caller (as `GaussianDensityCPU` does for
  radii), or let P1's `loadNeighbours` accept ordinal radii.
- **P9** inter-unit graph: return the module `emptyArray` (as `getConnectedUnits` does) instead of `[]`. Longer term
  per-unit CSR (`offsets: Int32Array`, `edgeIdx: Int32Array`) and struct-of-arrays edges; that changes the public API
  (record in `breaking-v6-changes.md`). Cantor keys become heap numbers above ~46k unit-local indices.

## Bugs

| Id | Title | File | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B1 | `approxNearest` puts the nearest candidate last; callers read index 0 | [grid.ts:441](../../packages/core/src/math/geometry/lookup3d/grid.ts#L441) | low | confirmed (repro) |
| B2 | `updateChunk = 1e5 / (maxRadius⁹ · scale)`: no yield for small radii | [molecular-surface.ts:382](../../packages/core/src/math/geometry/molecular-surface.ts#L382) | low | confirmed |
| B3 | `Sphere3D.copy` / `expandBySphere` keep stale extrema | [sphere3d.ts:59](../../packages/core/src/math/geometry/primitives/sphere3d.ts#L59) | low (latent) | confirmed (code) |
| B4 | RMSD state reuse aliases and overwrites a caller's centre | [minimize-rmsd.ts:92](../../packages/core/src/math/linear-algebra/3d/minimize-rmsd.ts#L92) | low (latent) | confirmed (code) |

- **B1** `approxQueryNearest` adds every candidate in the first non-empty cell that improves `minDistSq`, so the
  result is a decreasing-distance chain with the nearest last. `StructureLookup3D._approxNearest`
  ([lookup3d.ts:305](../../packages/model/src/model/structure/structure/util/lookup3d.ts#L305)) appends the same way.
  Consumers read `[0]`: blob-surface k-means (`blob-surface.ts:182`) and the external-structure theme
  (`plugin/core/src/themes/external-structure.ts:149`). Repro (`v/b1.ts`): points at 2.9 Å and 0.3 Å in one cell,
  `approxNearest(0,0,0,5)` → `indices [0, 1]`, `d² [8.41, 0.09]`: index 0 is the farther atom. Fix: track the best
  `(index, d²)` in the cell and add it once; same in `_approxNearest` (reset, add best). Blob k-means also never
  reassigns atoms farther than `blobSize/3` from every centroid (TODO at `:180`).
- **B2** `maxRadius⁹` was probably meant as a volume. Same formula in `gaussian-density/cpu.ts:65`. `maxRadius` 0.8,
  res 0.5 → chunk 373k atoms: a structure-level surface runs in one synchronous chunk, no progress, no abort. At the
  default VdW (max ~1.8–2) chunks are 100–250 atoms; coarse beads (> 3.6 Å) give chunk 1 (cheap `shouldUpdate`
  checks). Fix: chunk by work, e.g. `max(1, ceil(2e6 / (2·maxRadius·scale + 2)³))`. Results unchanged.
- **B3** `copy(out, a)` touches `out.extrema` only when `a` has extrema; `expandBySphere` keeps `sphere`'s extrema
  when `by` has none (or keeps `out`'s own when `out` is neither). `includeSphere`, `radiusSphere`, `Box3D.fromSphere3D`
  and `Sphere3D.includes` prefer extrema, so a reused `out` reports an old shape. No exercised path found: the scene
  sphere always gets extrema from `getSphere`, the cylinders/lines callers pass spheres with extrema, `addVec3` has no
  caller. Fix: `else out.extrema = undefined` in `copy`; clear in `expandBySphere` when either input lacks extrema.
- **B4** `resetRmsdState` assigns `state.centerA = data.centerA` when given; a later call on the same state without a
  centre writes the new centroid into the caller's vector, and both paths write the centre back into `data`. Only
  tm-align reuses state and never passes centres. Fix: own `Vec3`s in the state, `Vec3.copy` into them, no
  write-back.

## Checked and dropped

- P7 "getStructurePositionDataAndMaxRadius has the same pattern": no, its indices are `OrderedSet.ofRange(0, n)`.
- P7 "−14 % of the 1549–1577 MB allocated": the stage's allocation figure is JS heap only; typed-array backing stores
  are not in it, and the radius arrays are mostly untouched virtual pages. Impact high/medium → low.
- P6 "box identical, radius always encloses, only the centroid shifts": the box is identical; spheres merged from the
  extrema miss atoms by up to 0.88 Å. Impact medium → low, verdict plausible.
- P5 "needs strides on `Space` and a loop reorder": `dataOffset` already gives the strides; the closure is 3.1× of the
  4.1×.
- P2 auditor quoted Float64 copies (24 B per element); Float32 copies from Float32 conformations are exact (12 B).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/core-math/`,
bundled with `node ../build.mjs <in.ts> <out.mjs>` (esbuild, `molstar-src` condition). Fixture arguments: a path
(`data/examples/7qpd.fw2.cif` default, `../fixtures/1tup-2x2x2.bcif` = big-bcif); big runs need
`--max-old-space-size=8000`.

| Script | What |
| --- | --- |
| `bench-molsurf.ts` + `molsurf-opt.ts` | P1: original vs optimised, field/idField diff; `node bench-molsurf.mjs <file> <res> <probes> <reps>` |
| `bench-grid.ts` + `grid-opt.ts` | P2/P3: original vs bucket-ordered copies at counts 32/8/4, r 2/4/6.4, with radii |
| `v/bench-grid2.ts` + `v/grid-res.ts` | P2 variant: resolved indices without copies |
| `bench-fromarray.ts` | P4: helper vs inline, `poison` arg feeds 8 array kinds first |
| `bench-tensor.ts`, `v/bench-tensor2.ts` | P5: closure vs offsets, x vs z innermost |
| `bench-boundary.ts` + `boundary-coarse.ts` | P6: fine vs coarse timing, box and radius diff, unique extrema |
| `v/boundary-enclose.ts`, `v/boundary-enclose26.ts`, `v/boundary-flat.ts` | P6: atom enclosure of merged spheres (14 / 26 / 98 directions), flattened 98-direction pass |
| `bench-gauss.ts` + `gauss-opt.ts` | P8: separable kernel vs `fasterExp` |
| `est-radius.ts` | P7: sum of per-unit radius arrays |
| `v/b1.ts` | B1 repro |
