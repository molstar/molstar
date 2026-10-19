# Slice 12: tensor field access (`Space.get`, `set`, `dataOffset`, `getCoords`)

Every call site that reads or indexes a `Tensor` through its `Space` accessors in a loop, and where reading at strides
(`Tensor.strides`, #1987) would pay. Follows 08/P2 and 02/P5 (marching cubes), which found most of the cost in the
accessor closure and the rest in the memory order (02/P5: offsets alone 3.1×, offsets with z innermost 4.1×). A
monomorphic site inlines and is cheap, but a site that sees several spaces (two layouts, a floodfilled or wrapped space,
Float64 data from cube/dx) goes polymorphic and pays 4–8× until the page reloads. Unlike slices 00–11 this one was
written with the change: #1987 implements M1, G1 and G2. Measured on `90d717f49` + #1987, Node 25.6, Apple M4.

## The API

`Tensor.strides(space)` returns `[s0, s1, …]` when `get(data, c…) === data[Σ cᵢ·sᵢ]` and `dataOffset(c…) === Σ cᵢ·sᵢ`,
for any space made by `Tensor.Space` (all layouts of all ranks are linear); `undefined` otherwise. It is keyed by the
Space object, so the spread copies that override the accessors never have strides and stay on the accessors:
`Tensor.createFloodfilled` (`get`), `createWrappedTensor` and `createPeriodicVolume`
([util.ts:417](../../packages/graphics/src/repr/volume/util.ts#L417), all four accessors, modulo the dimensions). A
`get` or `dataOffset` replaced in place also drops them.

Pattern for a hot loop: take the strides once, keep the accessor loop as a separate function (G1 with both branches in
one loop ran at 55 ms, slower than the plain closure loop's 39; split, 27), carry the offset along the innermost axis.
Looping in memory order gains more on z-fastest data but changes the order of any accumulation; G1 keeps i innermost so
its magnitude stats stay bit-identical.

## Done in #1987

| Id  | Site                                                                                                        | Reads                                          | Measured                                                                                          |
| --- | ----------------------------------------------------------------------------------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| M1  | marching cubes, [algorithm.ts](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts)            | 8 `get` per cell + 14 per vertex               | 220³: z fastest 151 → 73 ms, x fastest 130 → 61, polymorphic 519 → 61, floodfilled 190–244 → 126 |
| G1  | `Grid.getGradients` → `computeGradients`, [grid.ts:320](../../packages/model/src/model/volume/grid.ts#L320) | 6 `get` + 1 `dataOffset` per voxel, whole grid | 160³: z fastest 39 → 27 ms, second layout in the session 143 → 18                                 |
| G2  | `Grid.makeGetInterpolatedGradient`, [grid.ts:445](../../packages/model/src/model/volume/grid.ts#L445)       | 8 `dataOffset` per sample                      | streamline trace 160³ (729 lines, 326k points): 38 → 38 ms monomorphic, 48 → 40 ms polymorphic    |

G1 is cached per grid; streamlines are its only caller today. G2 is one call per tracing step. The trace's own cost
dominates: a `Map` per line, a `Vec3` per point, `trilinearlyInterpolateVec3` and `lerp`, and one `dataOffset` per step
for the visit key ([basic.ts:76](../../packages/model/src/props/volume/streamlines/basic.ts#L76)), plus one per sphere
offset per point when marking occupancy (:212). Inlining the strided lookup into a strided `trilinearlyInterpolateVec3`
and keying by the strided offset is worth at most 10–20 % of the trace (S1).

## Candidates

Ranked by work per call; "whole grid" means one pass over every voxel each time it runs. Estimates are by analogy with
G1 (the same 7-read central-difference stencil), not measured.

| Id  | Site | Reads | When | Estimate | Effort |
| --- | ---- | ----- | ---- | -------- | ------ |
| F1  | `Tensor.floodfill3` / `floodfill2`, [tensor.ts:494](../../packages/core/src/math/linear-algebra/tensor.ts#L494), :435; `floodfilledGet`, :588 | 1 `dataOffset` per voxel, plus 1 `getCoords` and 6 `dataOffset` per visited voxel; `floodfilledGet` 1 `dataOffset` per read | every floodfilled molecular/Gaussian surface or isosurface (on the unfilled tensor, which has strides) | the best remaining one: whole grid, ~500 ms at 220³ (measured with #1987); neighbours at `o ± sᵢ` and coordinates carried instead of `getCoords` | S |
| V1  | `createVolumeTexture3d`, [util.ts:302](../../packages/graphics/src/repr/volume/util.ts#L302) | 7 `dataOffset` per voxel (value + normal), whole grid | direct volume (WebGL2/WebGPU), every rebuild | ~1.4× monomorphic, 5–8× when polymorphic | S |
| V2  | `createVolumeTexture2d`, [util.ts:199](../../packages/graphics/src/repr/volume/util.ts#L199) | 1 (`data`) or 7 (`normals`) `dataOffset` per voxel, whole grid | GPU isosurface (`data`), direct volume on WebGL1 (`normals`) | `normals` as V1; `data` small (the 2D-atlas index math dominates) | S |
| V3  | `createSegmentTexture2d`, [util.ts:367](../../packages/graphics/src/repr/volume/util.ts#L367), and `getSegmentCells`, [segment.ts:133](../../packages/graphics/src/repr/volume/segment.ts#L133) | 7 `dataOffset` per voxel of the segment's box, each through `set.includes`; `getSegmentCells` also a `segSet` (:149) per voxel | each segment of a segmentation | the `set.includes` linear scan (7 per voxel) costs more than the closure: a `Uint8Array` lookup by label first, strides second | S |
| T1  | `Grid.makeGetTrilinearlyInterpolated` → `trilinearlyInterpolate`, [grid.ts:157](../../packages/model/src/model/volume/grid.ts#L157), :191 | 8 `get` per sample through a `getValue` closure | `volume-value` color/size themes per vertex or element, volume slice per pixel, external-volume theme | per-vertex coloring of large meshes; a strided interpolator (base offset + 3 deltas), which must wrap periodic indices (12/B3) | S |
| V4  | `createVolumeSphereImpostor` / `createVolumeSphereMesh` / `createVolumePoint`, [dot.ts:178](../../packages/graphics/src/repr/volume/dot.ts#L178), :276, :372 | `get` + `dataOffset` per voxel, whole grid | volume dots and points | small: per-voxel `Vec3.transformMat4` and builder work dominate | S |
| B1  | `Grid.getIsosurfaceBoundingSphere`, [grid.ts:85](../../packages/model/src/model/volume/grid.ts#L85) (= 04/P7); segmentation bounds, [segmentation.ts:104](../../packages/model/src/formats/volume/segmentation.ts#L104) | `getCoords` (div/mod closure) per matching voxel, whole grid in memory order | CPU isosurface bounding sphere; loading a segmentation | strides do not help: loop the coordinates in memory order; 04/P7 measured 3–10× | S |
| E1  | volume-tools: `dataOffset` per voxel of a box or the grid in [soft-mask.ts:89](../../extensions/volume-tools/src/soft-mask.ts#L89), [mask-compute.ts:54](../../extensions/volume-tools/src/mask/internal/mask-compute.ts#L54), [segmentor mask-compute.ts:46](../../extensions/volume-tools/src/segmentor/internal/mask-compute.ts#L46), [voxel-labels.ts:42](../../extensions/volume-tools/src/voxel-labels.ts#L42), the flood fill in [volume-ops.ts:56](../../extensions/volume-tools/src/volume-ops.ts#L56); `getCoords` per selected voxel in [view-selection.ts:64](../../extensions/volume-tools/src/view-selection.ts#L64), [soft-mask.ts:32](../../extensions/volume-tools/src/soft-mask.ts#L32), [label-ops.ts:66](../../extensions/volume-tools/src/segmentor/internal/label-ops.ts#L66) | one accessor call per voxel | interactive masking and segmentation | moderate; same pattern | S |
| L1  | [slice.ts:511](../../packages/graphics/src/repr/volume/slice.ts#L511), :328, :534 | `get` / `dataOffset` per slice pixel | volume slice | low: 2D, and `color(location)` per pixel dominates | – |
| L2  | CIF tensor fields, reader [data-model.ts:375](../../packages/io/src/reader/cif/data-model.ts#L375) (`set` per value) and writer [encoder.ts:328](../../packages/io/src/writer/cif/encoder.ts#L328); CCP4 writer [ccp4.ts:99](../../packages/model/src/formats/writer/ccp4/ccp4.ts#L99); selection bounds [volume.ts:359](../../packages/model/src/model/volume/volume.ts#L359); periodic instance mapping [volume.ts:562](../../packages/model/src/model/volume/volume.ts#L562) | one accessor call per value or element | parse and export of CIF tensors, CCP4 export; per selected element or picked cell | low | – |

Already strided by hand (no accessor calls): the WebGPU isosurface and segment kernels get
`[dataOffset(1,0,0), dataOffset(0,1,0), dataOffset(0,0,1)]`
([isosurface.ts:335](../../packages/graphics/src/repr/volume/isosurface.ts#L335),
[segment.ts:324](../../packages/graphics/src/repr/volume/segment.ts#L324)), taken before wrapping, so correct today;
`Tensor.strides` says the same without relying on unit vectors through `dataOffset`.

## Layouts hand-coded outside the Tensor

These build a `Tensor.Space` over data they lay out by hand. Not a cost (they loop in memory order), but the layout is
the Tensor's to own: they should take their strides from `Tensor.strides(space)`.

- Surface fields, `[0, 1, 2]` with `iu = dimZ, iv = dimY`:
  [gaussian-density/cpu.ts:45](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L45),
  [molecular-surface.ts:371](../../packages/core/src/math/geometry/molecular-surface.ts#L371),
  [blob-surface.ts:510](../../packages/core/src/math/geometry/blob-surface.ts#L510),
  [sb-ncbr tunnels](../../extensions/sb-ncbr/src/tunnels/algorithm.ts#L588).
- FFT densities: [mtz.ts:205](../../packages/model/src/formats/volume/mtz.ts#L205),
  [structure-factors.ts:93](../../packages/model/src/formats/volume/structure-factors.ts#L93).
- Orbital grids: [alpha-orbitals data-model.ts:124](../../extensions/alpha-orbitals/src/data-model.ts#L124).

## Bugs

| Id  | Title | File | Failure scenario | Verdict |
| --- | ----- | ---- | ---------------- | ------- |
| B3  | Periodic trilinear interpolation reads unwrapped indices | [grid.ts:203](../../packages/model/src/model/volume/grid.ts#L203) | `trilinearlyInterpolate` skips the bounds check when `isPeriodic` but passes `i`, `ii`… unwrapped to `get`; the `volume-value` color/size and external-volume themes use the raw grid, so on MTZ or structure-factor maps (`periodicity: 'xyz'`) atoms outside the unit cell read the wrong voxel, or `undefined` → NaN → default color for negative indices (probe: `(5.5, 1, 1)` and `(-1.5, 1, 1)` on 4³ give 8 out-of-range reads). The slice is unaffected (it reads a `createPeriodicVolume` copy). Fix: wrap the indices modulo the dimensions | confirmed |
| B1  | The isosurface id field is read through the floodfilled `get` | `isosurface.ts` mesh and wireframe, before #1987 (now [getIsosurfaceFields](../../packages/graphics/src/repr/volume/isosurface.ts#L88)) | Floodfill `inside`: the id field was built with the floodfilled `scalarField.space`, so `floodfilledGet` also ran on the ids: a cell inside the surface (never visited) whose id is below `isoLevel` read `2 * isoLevel`, a wrong group for picking and coloring; large absolute iso values reach many cells. `outside` is not observable (the surfaces left only touch unvisited cells). Fixed in #1987: the id field keeps the unfilled space and is wrapped alike (`getIsosurfaceFields`) | fixed |
| B2  | `getCoords` of rank ≥ 4 writes a stray key | [tensor.ts](../../packages/core/src/math/linear-algebra/tensor.ts) `function getCoords` | Also writes `coords[ao[d + 1]]`, an index past the axis order (key `undefined`); the coordinates themselves are right | confirmed |

Not a bug, but misleading: `Grid.getGradients().magnitude` has `min`/`max` of |g| but `mean` of |g|² and
`sigma = sqrt(mean of |g|⁴)`; nothing reads it today.
