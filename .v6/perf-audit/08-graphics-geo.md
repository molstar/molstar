# Slice 08: `graphics/src/geo`

CPU geometry builders that every structure and volume visual goes through: `MeshBuilder` (primitive cache,
`addPrimitive`), the `Cylinders`/`Spheres`/`Lines`/`Points`/`Text` builders on `ChunkedArray`, marching cubes
(algorithm, builders), mesh color smoothing, `Mesh` utilities (normals, `uniformTriangleGroup`, `smoothEdges`), group
mapping, color/size/marker/transform data. The GPU gaussian density (`gaussian-density/gpu.ts`) is gl-bound and was not
audited. No `_spec`. Measured on `v6-opts-261019` (`90d717f49`), Node 25.6, Apple M4; fixtures and method in
[00-profile.md](00-profile.md). Each finding was re-read with its callers and re-measured by a second pass.

## Hot paths

| Path                                                       | Callers                                                                                         | Share                                                                     |
| ---------------------------------------------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| `calcMeshColorSmoothing`                                   | `processValues` of molecular/gaussian surface mesh, on every color, overpaint, transparency, emissive or substance change | 77 % self + 14 % `getIndex` of theme-apply on surfaces (7qpd 250 ms, big-bcif 808 ms / 160 MB) |
| `MarchingCubesState.processCell` / `get` / `interpolate`   | CPU gaussian surface, molecular surface, blob surfaces, volume isosurface and segment           | Tensor `get` 45 % of CPU gaussian surface, 15 % of molecular surface      |
| `MeshBuilder.addPrimitive`                                 | `addCylinder`, `addFixedCountDashedCylinder`, `addSphere`: interactions (mesh only), every mesh link/element visual, the no-impostor fallback | 77 % of `repr:interactions` (big-bcif 1.1 s)            |
| `ChunkedArray`-backed builders                             | all of the above; they size and keep the final buffers                                          | memory: see P3, P6, P8                                                    |

Already fine: sphere primitives are indexed and cached; the tube builder caches cos/sin; `getMarkersAverage` uses a
lookup table; `createGroupMapping` is a counting sort (but see B2); the CPU gaussian-density accumulation (core) loops in
memory order.

## Performance findings

| Id  | Title                                                                    | Where                                                                                       | Kind        | Impact | Estimate                                                       | Effort | Verdict   | Frozen |
| --- | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- | ----------- | ------ | -------------------------------------------------------------- | ------ | --------- | ------ |
| P1  | Color smoothing: sqrt on 216 cells, `getIndex` per cell, z innermost, sphere-cube grid, nothing reused | [color-smoothing.ts:128](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L128) | time+memory | high   | surface recolor ~2.5–3× (7qpd ~250 → ~90 ms); grid −2/3 with an AABB box | M      | confirmed | no     |
| P2  | Marching cubes: 8 closure reads per cell, no corner reuse, polymorphic `get` | [algorithm.ts:250](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L250) | time        | high   | MC 1.7×, 4–7× when the call site is polymorphic; CPU gaussian ~−25 % | S      | confirmed | no     |
| P4  | `addPrimitive` per vertex through megamorphic `Vec3.fromArray`/`transformMat4` | [mesh-builder.ts:122](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L122) | time        | high   | `addPrimitive` 1.5× (monomorphic) to 3× (app state); interactions ~−50 % | S      | confirmed | no     |
| P3  | Single-chunk buffers never trimmed; index reserve 3 tri/vertex           | [mesh-builder.ts:47](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L47)      | memory      | high   | cartoon −43 %, surface vertex buffers −20 %, mesh ball-and-stick −36 % (CPU and GPU) | S      | confirmed | no     |
| P7  | Cylinder impostor attributes repeated 6× per cylinder                    | [cylinders-builder.ts:116](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L116) | memory      | high   | cylinder buffers ~5× smaller (big-bcif 165 → ~35 MB)           | L      | plausible | yes    |
| P5  | Cylinder primitive keeps three.js UV duplicates (seam column, 16 cap centers) | [cylinder.ts:124](../../packages/graphics/src/geo/primitive/cylinder.ts#L124)              | time+memory | medium | capped cylinder 100 → 66 vertices; interactions buffers ~−27 % (~−110 MB big-bcif) | S      | confirmed | no     |
| P6  | Constant cylinder/lines mapping + index buffers allocated per geometry    | [cylinders-builder.ts:197](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L197) | memory      | medium | ~−50 MB JS heap, big-bcif ball-and-stick impostors             | S      | confirmed | no     |
| P8  | Cylinders/Lines builders take item counts as element counts              | [cylinders-builder.ts:65](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L65) | time+memory | medium | builder −25 % time, −50 % transient allocation                | S      | confirmed | no     |
| P10 | `smoothEdges` rebuilds `number[]` neighbour maps and Map edge counts 3+ times | [mesh.ts:683](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L683)                  | time+memory | low    | ~0.9 µs/vertex today; several × cheaper with CSR; matters once B1 is fixed | M      | confirmed | no     |
| P9  | Cylinder primitive cache keyed by exact radius, 32-bit hash, never evicts | [builder/cylinder.ts:43](../../packages/graphics/src/geo/geometry/mesh/builder/cylinder.ts#L43) | memory      | low    | bounded cache; no hash-collision risk                          | S      | plausible | no     |
| P11 | Spheres `positionGroup` filled twice on create                           | [spheres.ts:339](../../packages/graphics/src/geo/geometry/spheres/spheres.ts#L339)           | time        | low    | ≤ 1–2 ms at 186k spheres                                       | S      | confirmed | no     |

### P1: mesh color smoothing

- **Mechanism**: per sampled vertex, a 6×6×6 cell box; `Math.sqrt` before the `d > p` reject although only ~33 of
  216 cells lie in the radius-2 ball; `getIndex` (two `floor`, one `%`) per accepted cell; z innermost although x is the
  texture-fast axis; `count[index / itemSize]`. Grid box = `Box3D.fromSphere3D(boundingSphere)`, a cube of side 2R.
  `data`, `count`, `grid` allocated per call. Positions, groups, stride and resolution are color-independent, yet all of
  it reruns on every theme change (and for overpaint/transparency/emissive/substance, which share the function).
- **Evidence**: `cs.ts`, 600k vertices, stride 3, fake graphics (volume path): resolution 1: 211 → 56 ms; 0.6:
  287 → 100 ms; max output difference 0 (each cell receives at most one contribution per vertex, so the summation order
  is unchanged). `box.ts`: bounding-sphere cube / AABB = 2.97 over the 7qpd units (2.1–3.4 per unit). With the kernel
  91 % of the stage, the kernel alone gives ~2.5–3× on the stage, not 3.6×.
- **Fix**: (1) loops z/y/x with `dz²` and `dy²+dz²` early rejects, `sqrt` only when `d² ≤ p²`, `zOff[z]` table so
  `index = zOff[zi] + yi*width + xi`, colors hoisted. (2) Box from the vertex AABB (cache next to the bounding sphere;
  for instance types transform the corners, or take the smaller of both boxes), padding kept. (3) Grow-only module
  scratch for `data`/`count`. (4) Optional: cache the color-independent splat per (geometry version, resolution,
  stride) as a CSR of (cell, group, weight); recolor becomes O(nnz). This retains ~250 B per sampled vertex
  (~50 MB for 600k vertices at stride 3), so make it opt-in or aggregate weights per (cell, group). Same `zOff` trick in
  `getTrilinearlyInterpolated`.
- **Risk**: output must stay bit-identical (it is, for the kernel). The AABB changes `gridDim`/`gridTransform`, which are
  uniforms already; re-check `getFittingResolution` against `maxTextureSize`. Test group and groupInstance types and the
  four other `apply*Smoothing` variants.

### P2: marching cubes reads

- **Mechanism**: `processCell` calls `this.get` 8× per cell → `scalarFieldGet` closure → index math; every grid value is
  read 8× through two calls. Fields `Tensor.Space(dim, [0,1,2])` (gaussian, molecular, blob surfaces, cube, dx, mtz) are
  z-fastest while the loop is x-innermost (stride nY·nZ). `interpolate` adds 14 more closure reads per new vertex. The
  accessor closures differ per layout (and `floodfilledGet`, Float64 data from cube/dx), so the call site goes polymorphic
  in any session that has both a surface and a CCP4 map (`[2,1,0]`).
- **Evidence**: `mc.ts`, 220³ blobby field, re-run: separate processes current 106.6 / 91.8 ms (`[0,1,2]`/`[2,1,0]`),
  optimized 64.1 / 54.7 ms; same process (polymorphic) current `[2,1,0]` 426.7 ms, optimized 54.5 ms. Vertex count,
  triangle count, vertex checksum identical. The code already says `TODO cache scalarField differences for slices`.
- **Fix**: strides from `space.dataOffset` once (`sx`, `sy`, `sz`), direct `data[i*sx + j*sy + k*sz]` reads (idField
  too); carry the 4 corners of face i+1 to the next cell. Optional: loop the field's fastest axis innermost (permute
  coordinates, flip winding for odd permutations; ~15 % more); cache central-difference gradients for the two live slices.
- **Risk**: pure refactor. Test idField (−1/−2), negative isoLevel (normal and winding flip), `bottomLeft`/`topRight`
  ranges, the lines builder (`edgeFilter`).

### P4: `addPrimitive` vertex transform

- **Mechanism**: per vertex `v3fromArray`, `v3transformMat4`, `v3fromArray`, `v3transformMat3`, 3 ChunkedArray adds.
  `Vec3.fromArray`/`transformMat4` are called app-wide with `number[]`, Float32/Float64/Int32/Uint32 arrays; their
  feedback is per function, not per inlined call site, so they are megamorphic here too. The 16 matrix entries are
  reloaded per vertex.
- **Evidence**: profile: `Vec3.fromArray` 41 % self inside `addPrimitive` (77 % total) in `repr:interactions`.
  Re-measured with a fair bench (`verify/ap2.ts`: same matrix and cached primitive for both, 300k capped 16-segment
  cylinders, 30M vertices): current 463 ms monomorphic / 952 ms after warming `Vec3` with 6 array kinds; inlined 313 /
  315 ms; same output. The auditor's `ap.ts` compared `addFixedCountDashedCylinder` against a bare loop and overstated
  the monomorphic gain slightly; the polymorphic ~3× holds.
- **Fix**: hoist `m0..m15`, `n0..n8` into locals, read `va[i..i+2]`/`na[i..i+2]` directly in `addPrimitive` and
  `addPrimitiveFlipped`; keep the w-divide only for non-affine matrices. Type `Primitive` arrays as
  `Float32Array`/`Uint32Array` (all producers already do). Same for `transformPositionArray`/`transformDirectionArray`
  in [util.ts:26](../../packages/graphics/src/geo/util.ts#L26) (`Mesh.transform` on every surface and isosurface).
- **Risk**: floating-point order only. `addMesh` passes Float32 subarrays (same map).

### P3: oversized single-chunk buffers

- **Mechanism**: `createState` reserves `initialCount*3` triangles (3 per vertex; real meshes have 1–2).
  Caller estimates overshoot too (`link.ts:139` 64 vertices per half-link where 34 are used). `getMesh` calls
  `ChunkedArray.compact(…, true)`, which returns an underfilled single chunk as is; `gl.bufferData(array)` uploads the
  whole tail. Same in the Cylinders/Spheres/Points/Lines getters and `MarchingCubesMeshBuilder.get` (first chunk up to
  262,144 vertices with no estimate). With mesh reuse the oversized buffer becomes the next initial chunk, so the waste
  persists across updates.
- **Evidence**: `cap.ts` on 7qpd (meshes, no graphics), re-run: cartoon index 5.5/13.2 MB, vertex 7.8/10.2 MB, 10 MB
  unused (43 %); molsurf 3.1 MB, gaussian 3.3 MB (20 % of vertex buffers); ball-and-stick (mesh) 83.8 MB (36 %),
  spacefill (mesh) 21.9 MB (19 %); interactions 0. `lm.ts`, 50k solid links: index capacity 6.0× triangles, vertex 1.88×.
  Correction: ball-and-stick and spacefill use impostors in a browser with `fragDepth`; their numbers apply to the mesh
  fallback (no impostor support, headless capture, export). Cartoon and surfaces are always meshes.
- **Fix**: in the `get*` of every builder, trim a single chunk when `currentIndex < 0.8 * allocatedSize`
  (`chunk.slice(0, currentIndex)`, or a `ChunkedArray.compact` trim option); default index reserve ~2 triangles per
  vertex or a caller-given triangle estimate; `link.ts:139` → `(radialSegments + 1) * 2` per half-link plus caps.
- **Risk**: one copy when trimming (~10–20 ms per 150 MB). A trimmed reused buffer that is too small grows by chunks and
  is compacted once; correct. Nothing found relying on `buffer.length` ≥ previous capacity.

### P7: per-vertex cylinder attributes (frozen)

- **Mechanism**: `add()` writes start, end, group, scale, cap, colorMode 6× (one per quad-box vertex) because the gl
  schema reads them per vertex: 240 B per cylinder for 40 B of data, plus 72 B mapping and 48 B index.
- **Evidence**: profile: ball-and-stick big-bcif `aMapping`, `aStart`, `aEnd` 31.5 MB each; 165 MB buffers.
- **Fix**: after the WebGPU work, instanced draw with a shared 6-vertex template and per-cylinder attributes (divisor 1),
  or vertex pulling by `vertexId / 6` from a texture as spheres do with `tPositionGroup`; lines likewise (4×).
- **Risk**: shader and schema change in both backends, picking, merged render objects. Record only.

### P5: cylinder primitive duplicates

- **Mechanism**: `generateTorso` emits `radialSegments + 1` vertices per row (seam at 2π duplicates θ = 0, for UVs);
  `generateCap` emits `radialSegments` center vertices (one per face, for UVs) plus `radialSegments + 1` ring vertices.
  Mol\* uses no UVs. 16 segments, both caps: 34 + 2·33 = 100 vertices where 32 + 2·17 = 66 give the same triangles.
- **Evidence**: `lm.ts` default dashed links: 200 vertices per link (2 capped dashes). Interactions are mesh-only
  (no impostor visual), 421.6 MB buffers on big-bcif, `aPosition`/`aNormal` 143 MB each. Per dash 2,800 B vertex +
  768 B index → 1,848 + 768 B.
- **Fix**: `radialSegments` vertices per torso row with indices wrapping modulo, one center vertex per cap, same triangle
  order. Primitives are cached, so this is a one-time cost.
- **Risk**: keep the seam column when `thetaLength < 2π`. Check tests or exporters that assert vertex counts.

### P6: shared mapping and index buffers

- **Mechanism**: `getCylinders` allocates `Float32Array(n*18)` mapping and `Uint32Array(n*12)` elements filled with
  the same per-cylinder pattern for every geometry; `LinesBuilder` likewise (8 and 6 per line).
- **Evidence**: big-bcif ball-and-stick impostors: `aMapping` 31.5 MB + index ~21 MB of constant data.
- **Fix**: module-level grow-only buffers filled once up to the largest count seen; hand out `subarray(0, n*18)` /
  `subarray(0, n*12)`. WebGL `bufferData(subarray)` uploads only the view.
- **Risk**: subarrays share an ArrayBuffer; no in-place writer found in graphics, plugin, mvs or extensions
  (`mappingBuffer` is only read and ValueCell-updated). GPU memory unchanged unless the renderer also shares the buffer
  (frozen area).

### P8: builder sizes in the wrong unit

- **Mechanism**: a cylinder adds 6 ChunkedArray elements, a line 4, but `create(initialCount, chunkSize)` sizes chunks in
  elements. Callers pass item counts: `link.ts:384` and `polymer-backbone-cylinder.ts:81` (×2), `nucleotide-atomic-bond`
  (×15), `link.ts:594` lines, marching-cubes lines `builder.ts:108` (`linesCount`, chunk `linesCount/10` → ~30 extra
  chunks), `lines.ts:104`. The first chunk holds 1/3–1/4 of the data, more chunks follow, then `compact` copies all.
- **Evidence**: `cb.ts`, 300k cylinders, re-run: `create(2L, L/2)` 24.3 ms, `create(6L, …)` 18.0 ms, identical output;
  transient 12L·40 B vs 6L·40 B (144 vs 72 MB at L = 300k).
- **Fix**: make `initialCount`/`chunkSize` mean items inside the builders (×6, ×4), or fix the callers. Combine with
  the P3 trim so overestimates do not stay resident.
- **Risk**: semantic change for external callers (extensions, MVS `LinesBuilder.create(1024, 1024)`): changelog note.

### Low

- **P10**: `smoothEdges` (gaussian surface with `includeParent`) builds `getNeighboursMap` 3× and `getEdgeCounts` 2×
  (per-vertex `number[]` with `arraySetAdd`; `Map` keyed by `sortedCantorPairing`, heap numbers above ~46k vertices).
  `verify/se.ts`: 60 ms for a 66k-vertex open surface (~0.9 µs/vertex) with the fill loop exiting after one pass (B1).
  Fix with B1: CSR adjacency built once, border edges by sorted 64-bit keys or a typed hash, incremental update after
  fill.
- **P9**: `getCylinderPropsKey` rounds radii to 1e-3 and hashes 9 ints into a 32-bit FNV key for a module-global
  `cylinderMap` that never evicts; link radii are per link (`link.ts`), so continuous size themes make up to ~2,000
  radii × cap combinations, ~3 KB each. Collision probability is tiny (~5e-4 at 2k keys) but silent. Fix: unit-radius
  primitives keyed by (segments, caps, theta, top/bottom ratio) with the radius in `setCylinderMat`'s scale
  (`addPrimitive` already uses the inverse transpose for normals); structured key. Not measured on a real theme.
- **P11**: `fromArrays`/`update` call `shaderData.update()` (no LOD, a plain O(n) copy into the texture) and
  `createValues` calls it again with the real props; on update only one pass runs. Make the first one lazy.

## Bugs

| Id  | Title                                                         | Where                                                                                     | Severity | Verdict   |
| --- | ------------------------------------------------------------- | ----------------------------------------------------------------------------------------- | -------- | --------- |
| B1  | `fillEdges` never adds a triangle (`length < 2` filter)        | [mesh.ts:566](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L566)                 | low      | confirmed |
| B2  | `createGroupMapping` sized by item count, not max group id    | [util.ts:142](../../packages/graphics/src/geo/util.ts#L142)                               | low      | confirmed |
| B3  | `getTrilinearlyInterpolated` 0/0 on an integer grid coordinate | [color-smoothing.ts:265](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L265) | low      | confirmed |
| B4  | Marching-cubes lines indices in a Float32Array                | [builder.ts:70](../../packages/graphics/src/geo/util/marching-cubes/builder.ts#L70)        | low      | confirmed |
| B5  | `computeGaussianDensityTexture3d` dispatches `'2d'`           | [gaussian-density.ts:60](../../packages/graphics/src/geo/gaussian-density.ts#L60)          | low      | confirmed |

- **B1**: since `7e7993f5b` (2021) the candidates are `.filter(v => borderNeighboursMap.get(v)!.length < 2)` and then
  `bnd[1]` is read; before it was `nbs.length >= 2`. On a manifold border every vertex has 2 border neighbours, so there
  are no candidates; a 1-neighbour vertex gives `vC = NaN`, `dot = NaN`, `add` stays false. `fe.ts`: border-neighbour
  histogram `{2: 219}`, `smoothEdges` changes triangles by −2 (trim only). Effect: `includeParent` surface edges are
  never filled (visual quality only). Fix: `=== 2` (as `laplacianEdgeSmoothing`), a test asserting the triangle count
  grows on an open mesh; land with P10, since up to 10 fill iterations then run; review the visual change.
- **B2**: `bucketSizes`/`bucketFill` are `Int32Array(dataCount)` and offsets are filled for `i < dataCount`; group ids
  ≥ `dataCount` are dropped silently. `verify/gm.ts`: 50 spheres with groups 50..99 → `indices.length 0`, every range
  empty; a text label of group 10 with 4 chars → empty. The only consumer is `getGroupBoundingSphere` (shape loci
  bounds: focus/zoom on a shape group goes to the origin with radius 0). MVS primitives are not affected (dense groups
  per geometry); custom shapes with sparse ids (points/spheres for a subset, few-char labels) are. Fix: size by
  `maxId + 1`, fill offsets to `maxId + 1`; unit test.
- **B3**: `v1 - v0 = 0` on any axis where the scaled coordinate is an integer → `NaN` → `0` in the Uint8Array.
  `verify/gm.ts`: a vertex at (1, 1, 1) or (1.5, 2, 1.5) gets `[0,0,0,0]` while its neighbours get 200. Only the
  no-graphics path (`kind: 'vertex'`). Fix: `vd = v1 > v0 ? (v - v0) / (v1 - v0) : 0` per axis.
- **B4**: vertex ids are stored in `ChunkedArray.create(Float32Array, 2, …)` and cast `as Uint32Array`; ids above 2²⁴
  round, so wireframes of isosurfaces with > 16.7M vertices connect wrong vertices. Fix: `Uint32Array`; write directly
  into a `LinesBuilder` sized for `linesCount` lines (P8).
- **B5**: dead export (no callers in packages, extensions, apps) that returns a 2d texture. Pass `'3d'` (guarded by
  `texture3d` support) or remove it.

## Checked and dropped

- Nothing refuted. Corrections instead: P1 stage gain ~2.5–3× (the kernel is 91 % of it), not 3.6×; P3 ball-and-stick
  and spacefill numbers are the mesh fallback, not the browser default; P4 monomorphic gain 1.5× (fair bench), not 1.35×
  against a different workload; P11 is ≤ 1–2 ms, not "a few ms"; B2 does not reach MVS primitives.
- Aside for the core slice (not counted here): `ChunkedArray._compact` copies a partial last chunk with a scalar loop
  over the whole chunk length instead of `set(lastChunk.subarray(0, currentIndex))`
  ([chunked-array.ts:148](../../packages/core/src/data/util/chunked-array.ts#L148)).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/graphics-geo/`,
bundled with `node ../build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition), run with `node X.mjs`.

| Script          | What                                                                                      |
| --------------- | ----------------------------------------------------------------------------------------- |
| `cs.ts`         | color smoothing: current vs rewritten kernel, `node cs.mjs 600000 1 3` (vertices, resolution, stride) |
| `box.ts`        | 7qpd units: bounding-sphere cube vs AABB volume                                            |
| `mc.ts`         | marching cubes current vs strided/corner reuse, `node mc.mjs 220 ["current [0,1,2]"]` (name filter = separate process) |
| `ap.ts`         | dashed cylinders vs inlined `addPrimitive` (auditor; workloads differ)                     |
| `verify/ap2.ts` | `addPrimitive` vs inlined on the same primitive and matrix, `mono`/`poly`                   |
| `cap.ts`        | buffer capacity vs use per representation on 7qpd                                          |
| `lm.ts`         | link mesh vertices/triangles vs capacity per link style                                     |
| `cb.ts`         | `CylindersBuilder` sizing, `node --expose-gc cb.mjs 300000 2` vs `… 6`                       |
| `fe.ts`         | `fillEdges` candidates and triangle delta on an open MC mesh                                |
| `verify/se.ts`  | `smoothEdges` time on a 66k-vertex open surface                                             |
| `verify/gm.ts`  | B2 and B3 reproductions                                                                     |
