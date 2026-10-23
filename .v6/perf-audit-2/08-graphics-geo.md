# Slice 08 (v2): `graphics/src/geo`

The CPU geometry every structure and volume visual uses. That is marching cubes (algorithm and builders), `MeshBuilder`
and its tube/sheet/primitive adders, the `Cylinders`/`Spheres`/`Lines`/`Points`/`Text` builders on `ChunkedArray`,
mesh and texture-mesh color/layer smoothing, `Mesh` utilities (`smoothEdges`), group mapping, the font atlas, and the
new `clip-volumes.ts` (#2021). Frozen and not re-audited: `gl/` (merged copies, `calculateTransformBoundingSphere`),
the GPU gaussian density and compute smoothing. MVS load and animation cost is in `.v6/mvs-optimizations/research/`.
Measured at `0944823d4`, Node 25.6.1, Apple M4. Fixtures and the pipeline method are in [00-profile.md](00-profile.md).
v1 is [../perf-audit/08-graphics-geo.md](../perf-audit/08-graphics-geo.md). Each finding was re-read with its callers
by a second pass. The disputed ones were re-measured: P15, P12, P16, P17, B9.

## Hot paths now

| Path                                                       | Callers                                                                                   | Cost now                                                                                   |
| ---------------------------------------------------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `MarchingCubesState.interpolate`/`stridedRow`/`processCell` | CPU gaussian, molecular and blob surfaces, isosurfaces                                     | reads are strided (#1987). Cost now comes from re-tiering after a GC (P15): 120³ steady 38–45 ms, 18 MB, against 19 ms, 4.5 MB when kept warm |
| `addTube` / `addSheet` (cartoon)                           | `createPolymerTraceMesh` (repr)                                                            | scalar loops (#1993). Cold after every major GC (P12): 7qpd `repr:cartoon` 39 ms, 78 MB    |
| `calcMeshColorSmoothing` / `splatVertex`                   | every color, overpaint, transparency, emissive or substance change on smoothed surfaces   | tight kernel (#1991), only the instances reached (#2005, #2012). Re-splats the same weights on each recolor (P1): 7qpd 66–82 ms |
| `CylindersBuilder.add` / `addFixedCountDashes`             | impostor ball-and-stick, backbone, interactions (impostors since #2004)                    | 36 polymorphic ChunkedArray calls per cylinder (P16): flat-192 `add` 324 ms self            |
| `getFontAtlas`                                             | every `TextBuilder.create`: labels, measurements, MVS labels, camera helper               | one atlas per font, the 4 most recent kept (B6, fixed)                                    |

Since v1: #1987 (MC strides, corner reuse), #1990 (inline `addPrimitive`, cylinder without UV duplicates), #1991 (smoothing
kernel, AABB grid, scratch), #1993 (scalar cartoon loops), #2005 and #2012 (smoothing only over reached instances), #2007
(trimmed builders, item-count sizing), #2004 (interactions as cylinder impostors), #2021 (`clip-volumes.ts`).

## v1 → now

Fixed 7, partly fixed 2, still open 11, regressed 1, dropped 0.

| v1 id | Title                                                                  | Status       | PR                  | Note                                                                                              |
| ----- | ---------------------------------------------------------------------- | ------------ | ------------------- | ------------------------------------------------------------------------------------------------- |
| P1    | Color smoothing kernel: sqrt on 216 cells, `getIndex`, z innermost, cube grid, nothing reused | partly-fixed | #1991, #2005, #2012 | Parts (1)–(3) landed: `splatVertex` with early rejects, `zOff`, AABB grid, grow-only scratch. Part (4), the color-independent splat cache, is open (P1) |
| P2    | MC: 8 closure reads per cell, no corner reuse, polymorphic `get`        | fixed        | #1987               | `Tensor.strides` and `stridedRow` (algorithm.ts:291, :342), with a `get` fallback (`getRow` :373)       |
| P3    | Single-chunk buffers never trimmed; index reserve 3 tri/vertex         | fixed        | #2007               | `ChunkedArray.compactTrimmed` in every builder `get*` (36 uses); reserve of 2 tri/vertex (mesh-builder.ts:126) |
| P4    | `addPrimitive` per vertex through megamorphic `Vec3`                   | fixed        | #1990               | `addTransformedPrimitive` (mesh-builder.ts:38) hoists m0..m15 and skips the w-divide for affine matrices |
| P5    | Cylinder primitive keeps three.js UV duplicates                        | fixed        | #1990               | `isFull` gives `columnCount = radialSegments` (primitive/cylinder.ts:51), one center vertex per cap |
| P6    | Constant cylinder/lines mapping and index buffers per geometry         | still-open   |                     | cylinders-builder.ts:204–208, lines-builder.ts                                                    |
| P7    | Cylinder impostor attributes repeated 6× (frozen)                      | still-open   |                     | gl schema; recorded only                                                                          |
| P8    | Cylinders/Lines builders take item counts as element counts            | fixed        | #2007               | ×6 and ×4 inside `create` (cylinders-builder.ts:66–69); Text ×2                                    |
| P9    | Cylinder primitive cache keyed by exact radius, 32-bit hash, no eviction | still-open |                     | unchanged; one more caller with continuous radii, the clip gizmos (#2021)                         |
| P10   | `smoothEdges` rebuilds neighbour maps and edge counts                  | still-open   |                     | mesh.ts:683–697 unchanged                                                                          |
| P11   | Spheres `positionGroup` filled twice on create                         | still-open   |                     | spheres.ts:113/:121 and :339                                                                       |
| P12   | Cartoon curve, tube and sheet loops through `Vec3`                     | regressed    | #1993               | big entries faster (flat-32 208 → 168 ms). Small ones regressed (R1): the code is dropped at every major GC (P12) |
| B1    | `fillEdges` never adds a triangle                                      | still-open   |                     | mesh.ts:566 `< 2`                                                                                  |
| B2    | `createGroupMapping` sized by item count                               | still-open   |                     | util.ts:201–216                                                                                    |
| B3    | `getTrilinearlyInterpolated` 0/0                                       | still-open   |                     | color-smoothing.ts:665–670                                                                         |
| B4    | MC lines indices in a Float32Array                                     | still-open   |                     | marching-cubes/builder.ts:70                                                                       |
| B5    | `computeGaussianDensityTexture3d` dispatches '2d'                      | still-open   |                     | gaussian-density.ts:60                                                                             |
| 12/M1 | MC reads at tensor strides                                             | fixed        | #1987               | the same change as P2; 220³ 151 → 73 ms                                                            |
| 13/P1 | Instanced surfaces: layer smoothing spans every instance               | fixed        | #2005, #2012        | `getSplatInstances` / `getTextureMeshSplatInstances`; asm-864-low 12.6 s → 31 ms                  |
| 13/P3 | Unshared 5M allocation, builder over-allocation part                   | partly-fixed | #2007               | the geo part is fixed (P3, P8); repr temporaries are in other slices; cold-tier boxing (P12, P15) is new |
| 14/N4 | Cartoon meshes reserve ~1.4× their vertices                            | still-open   |                     | geo trims as designed; the overshooting estimate is the caller's (repr)                            |

## Open findings, ranked

| Rank | Id  | vs v1        | Title                                                                         | Where                                                                                                      | Kind        | Impact | Estimate                                                    | Effort | Verdict   | Frozen |
| ---- | --- | ------------ | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ----------- | ------ | ----------------------------------------------------------- | ------ | --------- | ------ |
| 1    | B6  | new; **fixed** | `getFontAtlas` keys its never-evicted cache by every Text prop (43 MB each)   | [text-builder.ts:50](../../packages/graphics/src/geo/geometry/text/text-builder.ts#L50)                    | bug         | high   | 43 MB retained per distinct label prop set → 1 per font     | S      | confirmed | no     |
| 2    | P15 | new          | MC re-tiers after every major GC (the optimized code embeds the per-run state map) | [algorithm.ts:117](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L117)              | time+memory | medium | 7qpd gaussian surface alloc 76 → 40 MB, −8 % time; MC up to 2× | S   | confirmed | no     |
| 3    | P12 | regressed    | Cartoon tube/sheet code dropped at every major GC (State map)                 | [mesh-builder.ts:120](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L120)                  | time+memory | medium | 7qpd cartoon 39 → 30 ms, 78 → 47 MB (v1 base)               | S      | confirmed | no     |
| 4    | P1  | partly-fixed | Surface recolor re-splats the color-independent weights                       | [color-smoothing.ts:414](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L414)            | time+memory | medium | recolor ~3× (7qpd 66–82 → ~25 ms)                           | M      | confirmed | no     |
| 5    | B8  | new          | Smoothed `volumeInstance` layers/colors not recomputed on a transform change  | [visual.ts:744](../../packages/graphics/src/repr/visual.ts#L744)                                           | bug         | medium | correctness of smoothed overpaint/transparency on moved surfaces | S   | confirmed | no     |
| 6    | P16 | new; **fixed** #2033 | `CylindersBuilder.add`: 36 polymorphic ChunkedArray calls per cylinder        | [cylinders-builder.ts:107](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L107)   | time        | medium | builder add ~3× (300k: 22–38 → 7–12 ms); flat-192 b&s −0.3 s | S     | confirmed | no     |
| 7    | P7  | still-open   | Cylinder impostor attributes repeated 6×                                      | [cylinders-builder.ts:120](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L120)   | memory      | high   | cylinder buffers ~5× smaller (5M: 4.46 GB)                  | L      | plausible | yes    |
| 8    | P6  | still-open   | Constant mapping and index buffers per geometry                               | [cylinders-builder.ts:204](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L204)   | memory      | medium | −120 B JS heap per cylinder (5M: ~−1.4 GB)                   | S      | confirmed | no     |
| 9    | P17 | new          | ChunkedArray adders: one polymorphic call site for every builder              | [chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59)                              | architecture| low    | ~1.5× variance in builder adds with the app's element kinds  | M      | plausible | no     |
| 10   | P9  | still-open   | Cylinder primitive cache keyed by radius, 32-bit hash, never evicts           | [builder/cylinder.ts:43](../../packages/graphics/src/geo/geometry/mesh/builder/cylinder.ts#L43)            | memory      | low    | bounded cache, no collisions                                | S      | plausible | no     |
| 11   | B7  | new; **fixed** #2033 | `addFixedCountDashedCylinder` mutates the caller's `props.topCap`             | [builder/cylinder.ts:141](../../packages/graphics/src/geo/geometry/mesh/builder/cylinder.ts#L141)          | bug         | low    | mirrored aromatic dashes lose their top caps                | S      | confirmed | no     |
| 12   | P11 | still-open   | Spheres `positionGroup` filled twice on create                                | [spheres.ts:113](../../packages/graphics/src/geo/geometry/spheres/spheres.ts#L113)                         | time        | low    | one O(n) pass per spacefill create                          | S      | confirmed | no     |
| 13   | P10 | still-open   | `smoothEdges` rebuilds neighbour maps per iteration                           | [mesh.ts:683](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L683)                                  | time+memory | low    | several × with CSR; matters once B1 is fixed                | M      | confirmed | no     |
| 14   | B1  | still-open   | `fillEdges` never adds a triangle                                             | [mesh.ts:566](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L566)                                  | bug         | low    | includeParent edges filled                                  | S      | confirmed | no     |
| 15   | B2  | still-open   | `createGroupMapping` sized by `dataCount`                                     | [util.ts:201](../../packages/graphics/src/geo/util.ts#L201)                                                | bug         | low    | focus on sparse shape groups                                | S      | confirmed | no     |
| 16   | B3  | still-open; **fixed** #2034 | `getTrilinearlyInterpolated` 0/0                                              | [color-smoothing.ts:670](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L670)            | bug         | low    | no-graphics smoothing colors                                | S      | confirmed | no     |
| 17   | B4  | still-open   | MC lines vertex ids in a Float32Array                                         | [builder.ts:70](../../packages/graphics/src/geo/util/marching-cubes/builder.ts#L70)                        | bug         | low    | wireframes above 16.7M vertices                             | S      | confirmed | no     |
| 18   | B5  | still-open   | `computeGaussianDensityTexture3d` dispatches '2d'                             | [gaussian-density.ts:60](../../packages/graphics/src/geo/gaussian-density.ts#L60)                          | bug         | low    | API correctness (no callers)                                | S      | confirmed | no     |

### P15: marching cubes re-tiers after every major GC (R2)

- **Mechanism**: each `computeMarchingCubesMesh`/`Lines` makes a new `MarchingCubesState` (:117) and new builder
  closures (builder.ts:24). TurboFan code of `interpolate`, `stridedRow`, `processCell` and `get3dOffsetFromEdgeInfo`
  embeds the instance map weakly. The map is reached only through weak transitions, so it lives only while some
  instance does. The first major GC after a run clears it and deoptimizes the code. The next surface runs in Ignition,
  Sparkplug and Maglev again, and they box every double in `interpolate`. #1987 moved the boxing from the base's
  `get` closures into `interpolate`, which is R2's "interpolate 42 MB, 4× base".
- **Evidence** (re-run, `mc.ts`, 120³ blobby field, 145k vertices): run 0 28.5 ms / 9.1 MB, then 38–49 ms / 17.6–18.2 MB
  per run (`interpolate` 15.6 MB). `--trace-deopt`: "embedded weak objects cleared" for these 4 functions once per run.
  With the last state kept alive (`mc-last.mjs`): 18.6–23.6 ms, 4.4–6.9 MB. Keeping the builder or ctx alive does not
  help. Auditor's pipeline run on 7qpd with a patched bundle (`benchbuild.mjs`, 10 runs): `repr:gaussian-surface`
  354/348 → 318/328 ms median, 76/75 → 40/40 MB (v1 base 52 MB). molecular surface 155 → 120 MB, time within
  noise. Small grids (`mcs.ts`, 40 × 48³): 61 → 44 ms, 23 → 12 MB. Measured in Node. Chrome runs the same V8
  heuristic, but this was not measured in a browser.
- **Fix**: keep the map alive. Simplest is a module-level anchor: one `MarchingCubesState` built at module init on a
  1³ field. Or a small pool acquired in the `MarchingCubesComputation` constructor and released in a `finally` after
  `run()`. One shared instance is not enough, because `doSlices` awaits `ctx.update` and two surfaces interleave.
  Give every field its final representation (`isoLevel` a double, `idField` always assigned). A kept or pooled state
  must drop its `scalarField`/`idField`/builder references after a run, or it retains the last field (7 MB at 120³).
  Pin it with a perf-tier allocation check: repeated small MC builds allocate below about 2× their output buffers.
- **Risk**: V8 behaviour, so pin it with a measurement, not a comment. Output bit-identical.

### P12: cartoon tube/sheet code dropped at every major GC (R1)

- **Mechanism**: `createState` (mesh-builder.ts:120) returns a fresh literal. `addTube`/`addSheet` are optimized
  against its map, which a major GC clears once no State is alive. The next small build then runs `addTube` in the
  lower tiers, and they box the scalar locals of #1993. That is R1 ("addTube 36 MB, base 20"), and it shows only from
  the second build on.
- **Evidence** (re-run, `tube.ts`, 3000 tubes of 8×16, gc between builds): cold and warm runs alternate, 22.5–24.2 ms
  with 34.5–34.9 MB and 5.9–6.5 ms with 0 MB. Auditor: the v1 base shows 25 ms / 23 MB against 9.8 ms. With a
  module-level `Anchor = MeshBuilder.createState(1, 1)` (`tube-a.ts`), every run takes 6.4–7.2 ms / 0 MB. Pipeline
  7qpd, patched bundle (anchor plus P15): `repr:cartoon` 39.7/38.5 → 30.3/29.4 ms, 80/76 → 48/46 MB (v1 base 32 ms /
  38 MB). The rest of R1 is likely `interpolateNormals` (repr `curve-segment.ts`, slice 07), which shows the same
  deopt reason.
- **Fix**: `const StateMapAnchor = createState(1, 1)` at module scope, with a comment citing this finding. Do the same
  for the per-build objects in curve-segment (slice 07). Add an allocation check for repeated small cartoon builds.
  This also helps every other `MeshBuilder` user rebuilt between GCs: link meshes, shapes, playback frames.
- **Risk**: the anchor retains ~50 B and 4 one-element chunks. Output unchanged.

### P1: recolor re-splats the color-independent weights

- **Mechanism**: positions, groups, stride, resolution and transforms do not depend on the colors. Each theme or layer
  change still re-runs `splatVertex` for every sampled vertex: about 33 cells of the radius-2 ball, a sqrt each, the
  position read through megamorphic `Vec3.fromArray`/`transformMat4Offset` (:334–337). It then normalizes the region
  into a fresh full-grid `Uint8Array` (:444).
- **Evidence**: 00-profile theme-apply surfaces 7qpd 66–82 ms (`splatVertex` 66 % self), flat-32 273–319 ms. Auditor's
  `cs.ts`/`cs2.ts` (145k vertices, 2000 groups, the real function with a stub graphics): 35.4–36.2 ms (res 1, stride 3)
  and 77.2 ms (res 0.7, stride 2). A cached CSR recolor takes 10.2–12.6 ms and 26–29 ms. Building the CSR once costs
  33–52 ms and 13.4–19.5 MB (34 nnz per sample), and only the 432k touched cells of 2.0M are normalized. Bit-identical
  (same `c·s` products, same order).
- **Fix**: cache per (geometry version, resolution, stride, transform version for instance types) a CSR of (cell,
  weight) per sample plus the touched-cell list. Recolor then loops over nnz into a reused grid. Aggregating per
  (cell, group) gives ~5× and ~5 MB, but changes the float order. Make it opt-in or bound it by vertex count. A
  cheaper first step is to inline the position read and hoist the 12 instance-matrix entries in `splatVertex`.
- **Risk**: memory (8 B per nnz). The cache must be invalidated on geometry and transform changes, the same hook B8
  needs.

### P16: `CylindersBuilder.add` per-element adds

- **Mechanism**: `add` (:107) writes 6 vertex copies with 6 × (2 `add3` + 4 `add`) = 36 ChunkedArray calls (60
  element writes), each with its own capacity check. The adders are shared by every builder in the app (P17). This is
  the path for impostor ball-and-stick, backbone and, since #2004, interactions (`addFixedCountDashes`).
- **Evidence** (re-run, `verify/cyl2.ts`, 300k cylinders; the auditor's `local` variant skipped 4 of 6 compactions,
  this one compacts all 6). Builder `add` phase: 22.3 ms clean and 33.7–38.8 ms after warming with the app's kinds.
  The local per-cylinder writer, all compactions included, takes 7.2 ms and 12–17 ms. Output is identical. Profile
  flat-192 `repr:ball-and-stick`: `add` 324 ms self, plus `add3` 190 ms and `add` 150 ms in ChunkedArray.
- **Fix**: one capacity check per cylinder across the 6 arrays (they grow alike), 6 copies written straight into
  `currentChunk`, then bump `currentIndex`/`elementCount`. Fall back to the generic path when any chunk would cross its
  end. Keep `addFixedCountDashes` on scalar locals.
- **Risk**: chunk alignment of the six arrays (assert or fall back). Output identical.

### B8 (bug): smoothed `volumeInstance` grids go stale on a transform change

See Bugs.

### Medium and low, short

- **P7** (frozen): `add` writes start/end/group/scale/cap/colorMode per vertex, 240 B per cylinder for 40 B of data.
  13: ball-and-stick buffers 0.99 GB at 1M atoms, 4.46 GB at 5M. After the WebGPU work: instanced draw or vertex
  pulling by `vertexId / 6`. Plausible: the gain depends on a schema change.
- **P6**: `getCylinders` allocates `Float32Array(n·18)` mapping and `Uint32Array(n·12)` indices with a constant
  pattern. It reuses them only when the count fits the previous geometry. LinesBuilder does the same. 7qpd `aMapping`
  is 2.1 MB; 13 measured 850 MB at 5M. Fix: module grow-only buffers handed out as `subarray`. Nothing writes them in
  place (v1).
- **P17**: `add`/`add2`/`add3` (chunked-array.ts:59–90) have per-function feedback for all builders. The production
  element kinds are Float32, Uint32, Int32 and Uint8 (no Array or Float64 outside tests), and they cost the cylinders
  builder 22 → 34–38 ms (1.5×). The auditor's 450 ms needs Array and Float64 stores through all three adders, which
  no production path or test does (see B9, dropped). Fix: element-kind-specific writers in core, or P16's direct
  writes. Low.
- **P9**: `cylinderMap` (builder/cylinder.ts:17) rounds radii to 1e-3, hashes with 32-bit FNV and never evicts. Since
  #2021 the clip gizmos (`clip-volumes.ts`) also call `addCylinder` with continuous radii (cylinder `scale`, cone rim
  = `reach·tan`, edge `width`), 48 segments each. Dragging a gizmo leaves one cached primitive per distinct radius.
  Fix: unit-radius primitives keyed by (segments, caps, theta, top/bottom ratio), the radius in the matrix.
- **P11**: `fromArrays`/`update` run `shaderData.update()` with the default LOD, and `createValues` (:339) runs it again
  with the real props. Make the first one lazy.
- **P10**: `smoothEdges` builds `getNeighboursMap` 3× and `getEdgeCounts` 2× (number[] per vertex, a `Map` keyed by
  Cantor pairing), about 0.9 µs per vertex. Land a CSR adjacency together with B1.

## Bugs

| Id  | vs v1      | Severity | Title                                                                       | Where                                                                                              | Verdict   |
| --- | ---------- | -------- | --------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | --------- |
| B6  | new; **fixed** | high     | `getFontAtlas` caches by every Text prop: 43 MB per distinct set, never evicted | [text-builder.ts:50](../../packages/graphics/src/geo/geometry/text/text-builder.ts#L50), [font-atlas.ts:16](../../packages/graphics/src/geo/geometry/text/font-atlas.ts#L16) | confirmed |
| B8  | new        | medium   | Smoothed `volumeInstance` colors/layers not recomputed by `setTransform`    | [visual.ts:744](../../packages/graphics/src/repr/visual.ts#L744)                                   | confirmed |
| B7  | new; **fixed** #2033 | low      | `addFixedCountDashedCylinder` mutates the caller's `props.topCap`           | [builder/cylinder.ts:141](../../packages/graphics/src/geo/geometry/mesh/builder/cylinder.ts#L141)  | confirmed |
| B1  | still-open | low      | `fillEdges` never adds a triangle (`< 2` filter, then `bnd[1]`)             | [mesh.ts:566](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L566)                          | confirmed |
| B2  | still-open | low      | `createGroupMapping` sized by `dataCount`, not `maxId + 1`                  | [util.ts:201](../../packages/graphics/src/geo/util.ts#L201)                                        | confirmed |
| B3  | still-open; **fixed** #2034 | low      | `getTrilinearlyInterpolated` 0/0 on integer grid coordinates                | [color-smoothing.ts:670](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L670)    | confirmed |
| B4  | still-open | low      | MC lines builder stores vertex ids in a Float32Array                        | [builder.ts:70](../../packages/graphics/src/geo/util/marching-cubes/builder.ts#L70)                | confirmed |
| B5  | still-open | low      | `computeGaussianDensityTexture3d` dispatches '2d' (latent)                  | [gaussian-density.ts:60](../../packages/graphics/src/geo/gaussian-density.ts#L60)                  | confirmed |

- **B6**: `TextBuilder.create` merges all `Text.Params` with the caller's props and passes the result to `getFontAtlas`.
  That function keys a module object by `JSON.stringify(props)` and never clears it. Each new key builds a `FontAtlas`
  whose constructor allocates the full texture, 350 × lineHeight × maxWidth bytes: 43.1 MB at the default
  `fontQuality` 3 (6565²). Each one then redraws its glyph SDFs. Auditor's `font.ts`: 8 prop sets gave 7 distinct
  43 MB atlases. Failure scenarios:
  - Labels: drag a label's `sizeFactor`, `backgroundMargin` or `tetherLength` slider through 20 values. That is 20
    atlases, 0.86 GB retained for the session.
  - Measurements: distance labels with different `customText`.
  - MVS: the custom-label visual's props carry `items` (every inline label's text and position, load.ts:659). Every
    snapshot of a labelled story with different labels gets its own atlas, and so does every structure.
  - Fix: `getFontAtlas({ fontFamily, fontQuality, fontStyle, fontVariant, fontWeight })`, or pick the five fields
    inside `getFontAtlas`. Add a test that props differing only outside the font fields share one atlas.
  - **Fixed**: `getFontAtlas` keys by the five font fields and gives the atlas only those; the cache is an LRU of 4
    (an evicted atlas lives on in the texts built with it). Test: `text/_test/font-atlas.test.ts`.
- **B8**: for colors or layers of an instance type, the smoothing grid is built in `aTransform` space
  (`getVertexBox`, `splatVertex` :335; texture-mesh likewise), and the shader samples it at `vModelPosition`
  (= uModel·aTransform·position). Overpaint, transparency, emissive and substance smoothing is always `groupInstance`
  → `volumeInstance`, so this affects every smoothed layer on any surface, single-instance ones too. `aTransform`
  includes the representation `matrix` (`updateTransformData`). `Visual.setTransform` rewrites `aTransform` and the
  bounds, but never re-smooths. units-representation `setVisualState` also calls `setLayers` (:389) before
  `setTransform` (:395). Failure scenarios:
  - Overpaint one chain of a gaussian surface (smoothed by default) and move the structure by hand (Prism, #2021) or
    by `setState({ transform })`. The paint stays where the structure was, and outside the grid the unclamped
    `texture3dFrom2dLinear` reads other slices.
  - With `unitTransforms`, e.g. an assembly unwind animation, overpaint stays in the assembled positions.
  - The CPU (no-graphics) path is unaffected.
  - Fix: re-run the instance-type smoothings in the visual's `setTransform` (debounced during drags), and apply the
    transform before `setLayers`.
- **B7**: for an odd `segmentCount` with `stubCap` false, the last dash sets `props.topCap = false` on the caller's
  object and never restores it. Failure scenario: `link.ts` MirroredAromatic (:304, :309) passes the same
  `cylinderProps` twice. With `aromaticDashCount` 1, 3 or 5 and `linkStub` false, every dash of the mirrored row has no
  top cap. polymer-gap uses 10 segments (even), so it is not hit. Fix: a local `{ ...props, topCap: false }` for the
  last dash.
- **B1**: `.filter(v => borderNeighboursMap.get(v)!.length < 2)`, then `bnd[1]` is read. A manifold border has no
  candidates, and a 1-neighbour vertex gives NaN. So `includeParent` surface edges are never filled. Fix: `=== 2`, and
  land it with P10.
- **B2**: `bucketSizes`/`bucketFill` are `Int32Array(dataCount)`. Failure scenario: a shape of 50 spheres with groups
  50..99 gets an empty mapping, and focusing on a group goes to the origin with radius 0. Fix: size by `maxId + 1`.
- **B3**: `(v − floor v)/(ceil v − floor v)` is 0/0 at integral coordinates. A vertex at (1, 1, 1) gets [0,0,0,0] on
  the no-graphics path. Fix: `v1 > v0 ? … : 0` per axis.
- **B4**: ids stored in `ChunkedArray.create(Float32Array, 2, …)` and cast `as Uint32Array`, so ids above 2²⁴ round.
  Failure scenario: a wireframe isosurface with more than 16.7M vertices joins the wrong vertices. Fix: `Uint32Array`.
- **B5**: passes `'2d'` to `_computeGaussianDensityTexture`, so any caller gets a 2d texture. There are no callers in
  packages, extensions or apps. Fix: `'3d'` or remove the export.

## Checked and dropped

- **B9** (auditor: the Vitest `chunked-array.test.ts` makes later builders megamorphic, 20×): refuted. The test stores
  Array and Int32 through `add`/`add2` only, and Float32 through `add3`. Replayed exactly (`verify/warmtest.ts`),
  `CylindersBuilder` takes 38 ms, the same as the app's own element kinds, not 450 ms. The 450 ms needs Array and
  Float64 stores through all three adders, which the test does not do. No geometry builder has an `itWithin` or
  timing test to flake.
- Corrections: P16 is 36 calls (60 writes) per cylinder, not 54. Its gain is ~3× on the add phase in a fair bench (the
  auditor's local variant compacted 2 of 6 arrays). P17's real-app variance is ~1.5×, not 1.5–20×, so it is low. B8
  covers every smoothed layer, not just instanced surfaces. B6 is amplified by MVS custom labels (`items` in the props).
  For P15, the auditor's single shared state is not enough because surfaces interleave, and a kept state must drop its
  field references.
- `clip-volumes.ts` (new, #2021): gizmo meshes of a few hundred vertices built per drag; the cost is fine apart from
  the P9 cache growth.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/graphics-geo/`.
Bundle with `node build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition) and run with `node --expose-gc X.mjs`.
`alloc.ts` samples the JS heap through the inspector, and `warm.ts` warms ChunkedArray (`WARM=real|real4|1`).

| Script                         | What                                                                                     |
| ------------------------------ | ---------------------------------------------------------------------------------------- |
| `mc.ts`, `mc-last.mjs`, `mc-keep.mjs`, `mcopt/` | P15: MC 120³ repeated runs, current vs the last or every state kept alive vs a pooled state; add `--trace-deopt` |
| `mcs-cur.ts`, `mcs-opt.ts`     | P15: 40 small 48³ grids per rep, gc between reps                                          |
| `tube.ts`, `tube-a.ts`, `tube-k.ts`, `tube-l.ts`, `tube-v*.ts` | P12: 3000 tubes per build, current vs a module anchor or the last state kept vs the v1 base |
| `benchbuild.mjs`, `bb/`        | P12 and P15: pipeline bench bundle, base and patched (esbuild onLoad patch)               |
| `cs.ts`, `cs2.ts`              | P1: `calcMeshColorSmoothing` vs a cached-CSR recolor on a 145k-vertex MC mesh             |
| `cyl.ts`, `verify/cyl2.ts`, `verify/cyl3.ts`, `verify/warmtest.ts` | P16, P17, B9: CylindersBuilder vs a per-cylinder writer; `WARM=test` replays chunked-array.test.ts |
| `font.ts`                      | B6: distinct atlases per Text prop set (stub canvas)                                      |
