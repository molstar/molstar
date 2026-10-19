# Slice 12: tensor field access (`Space.get`, `dataOffset`, `getCoords`)

Every call site that reads or indexes a `Tensor` through its `Space` accessors, and where reading at strides
(`Tensor.strides`, #1987) would pay. Follows 08/P2 and 02/P5 (marching cubes), which showed the cost is the accessor
closure, not the memory order: a monomorphic site inlines and is cheap, but one site that sees several spaces (two
layouts, a floodfilled or wrapped space, Float64 data from cube/dx) goes polymorphic and pays 4–8× until the page
reloads. Measured on `90d717f49` + #1987, Node 25.6, Apple M4.

## The API

`Tensor.strides(space)` returns `[s0, s1, …]` when `get(data, c…) === data[Σ cᵢ·sᵢ]` and `dataOffset(c…) === Σ cᵢ·sᵢ`,
for any space made by `Tensor.Space` (all layouts of all ranks are linear); `undefined` otherwise. It is keyed by the
Space object, so the spread copies that override the accessors never have strides and stay on the accessors:
`Tensor.createFloodfilled` (`get`), `createWrappedTensor` and `createPeriodicVolume`
([util.ts:417](../../packages/graphics/src/repr/volume/util.ts#L417), all four accessors, modulo the dimensions). A
`get` or `dataOffset` replaced in place also drops them. Pattern for a hot loop: take the strides once, keep the
accessor loop as a separate function (G1 with both branches in one loop ran at 55 ms, slower than the plain closure
loop's 39; split, 27), carry the offset along the innermost axis.

## Done in #1987

| Id  | Site                                                                                                                   | Reads                                                         | Measured                                                                                      |
| --- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| M1  | marching cubes, [algorithm.ts](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts)                       | 8 `get` per cell + 14 per vertex                              | 220³: z fastest 151 → 73 ms, x fastest 130 → 61, polymorphic 519 → 61, floodfilled 190–244 → 126 |
| G1  | `Grid.getGradients` → `computeGradients`, [grid.ts:320](../../packages/model/src/model/volume/grid.ts#L320)            | 6 `get` + 1 `dataOffset` per voxel, whole grid                | 160³: z fastest 39 → 27 ms, second layout in the session 143 → 18                             |
| G2  | `Grid.makeGetInterpolatedGradient`, [grid.ts:441](../../packages/model/src/model/volume/grid.ts#L441)                  | 8 `dataOffset` per sample                                     | streamline trace 160³ (729 lines, 326k points): 38 → 38 ms monomorphic, 48 → 40 ms polymorphic |

G1 is cached per grid; streamlines are its only caller today. G2 is one call per tracing step; the trace's own cost (a `Map` per line, a `Vec3` per point,
`trilinearlyInterpolateVec3` and `lerp`) dominates, so inlining the strided lookup into a strided
`trilinearlyInterpolateVec3` is worth at most 10–20 % of the trace (T1 below).

## Candidates

Ranked by work per call; "whole grid" means one pass over every voxel each time the visual is built. Estimates are by
analogy with G1 (the same 7-read central-difference stencil), not measured.

| Id  | Site                                                                                                                                                                                            | Reads                                                                         | When                                                          | Estimate                                                                    | Effort |
| --- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------------------- | ------ |
| V1  | `createVolumeTexture3d`, [util.ts:302](../../packages/graphics/src/repr/volume/util.ts#L302)                                                                                                    | 7 `dataOffset` per voxel (value + normal), whole grid                         | direct volume (WebGL2/WebGPU), every rebuild                  | ~1.4× monomorphic, 5–8× when polymorphic                                     | S      |
| V2  | `createVolumeTexture2d`, [util.ts:199](../../packages/graphics/src/repr/volume/util.ts#L199)                                                                                                    | 1 (`data`) or 7 (`normals`) `dataOffset` per voxel, whole grid                | GPU isosurface (`data`), direct volume on WebGL1 (`normals`)  | `normals` as V1; `data` small (the 2D-atlas index math dominates)            | S      |
| V3  | `createSegmentTexture2d`, [util.ts:367](../../packages/graphics/src/repr/volume/util.ts#L367), and its copy in [segment.ts:133](../../packages/graphics/src/repr/volume/segment.ts#L133)        | 7 `dataOffset` per voxel of the segment's box, each through `set.includes`    | each segment of a segmentation, GPU path                      | the `set.includes` linear scan (7 per voxel) costs more than the closure: a `Uint8Array` lookup by label first, strides second | S      |
| V4  | `createVolumeSphereImpostor` / `createVolumeSphereMesh` / points, [dot.ts:178](../../packages/graphics/src/repr/volume/dot.ts#L178), :276, :372 (spheres impostor, spheres mesh, points)                                                 | `get` + `dataOffset` per voxel, whole grid                                    | volume dots and points                                        | small: per-voxel `Vec3.transformMat4` and builder work dominate               | S      |
| T1  | `Grid.makeGetTrilinearlyInterpolated`, [grid.ts:157](../../packages/model/src/model/volume/grid.ts#L157), and `trilinearlyInterpolateVec3`                                                      | 8 `get` per sample through a `getValue` closure                               | `volume-value` color/size themes per vertex or element, volume slice per pixel, external-volume theme | per-vertex coloring of large meshes; a strided interpolator (base offset + 3 deltas) | S      |
| B1  | `Grid.getIsosurfaceBoundingSphere`, [grid.ts:85](../../packages/model/src/model/volume/grid.ts#L85); segmentation bounds, [segmentation.ts:104](../../packages/model/src/formats/volume/segmentation.ts#L104) | `getCoords` (div/mod closure) per matching voxel, whole grid in memory order  | isosurface bounding sphere; loading a segmentation            | strides do not help here: loop the coordinates in memory order (or derive them from the offset with the strides) | S      |
| E1  | volume-tools: [soft-mask.ts:89](../../extensions/volume-tools/src/soft-mask.ts#L89), [mask-compute.ts:54](../../extensions/volume-tools/src/mask/internal/mask-compute.ts#L54), [segmentor mask-compute.ts:46](../../extensions/volume-tools/src/segmentor/internal/mask-compute.ts#L46), [voxel-labels.ts:42](../../extensions/volume-tools/src/voxel-labels.ts#L42), [volume-ops.ts:19](../../extensions/volume-tools/src/volume-ops.ts#L19) | `dataOffset` per voxel of a box or the grid; flood fill with `getCoords`      | interactive masking and segmentation                           | moderate; same pattern                                                       | S      |
| L1  | [slice.ts:511](../../packages/graphics/src/repr/volume/slice.ts#L511), :328, :534                                                                                                              | `get` / `dataOffset` per slice pixel                                          | volume slice                                                  | low: 2D, and `color(location)` per pixel dominates                           | –      |
| L2  | CCP4 writer [ccp4.ts:99](../../packages/model/src/formats/writer/ccp4/ccp4.ts#L99); CIF tensor fields [encoder.ts:328](../../packages/io/src/writer/cif/encoder.ts#L328); selection bounds [volume.ts:359](../../packages/model/src/model/volume/volume.ts#L359) | one accessor call per value                                                   | export; per selected element                                  | low                                                                          | –      |

Already strided by hand (no accessor calls): the WebGPU isosurface and segment kernels get
`[dataOffset(1,0,0), dataOffset(0,1,0), dataOffset(0,0,1)]`
([isosurface.ts:333](../../packages/graphics/src/repr/volume/isosurface.ts#L333),
[segment.ts:324](../../packages/graphics/src/repr/volume/segment.ts#L324)), taken before wrapping, so correct today;
`Tensor.strides` says the same without relying on unit vectors through `dataOffset`.

## Layouts hand-coded outside the Tensor

The field producers build `Tensor.Space(dim, [0, 1, 2])` and then index the data by hand (`iu = dimZ, iv = dimY`):
[gaussian-density/cpu.ts:53](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L53),
[molecular-surface.ts](../../packages/core/src/math/geometry/molecular-surface.ts),
[blob-surface.ts](../../packages/core/src/math/geometry/blob-surface.ts),
[sb-ncbr tunnels](../../extensions/sb-ncbr/src/tunnels/algorithm.ts). Not a cost (they loop in memory order), but the
layout is the Tensor's to own: they should take their strides from `Tensor.strides(space)`.

## Not caused by strides, found on the way

- [isosurface.ts:188](../../packages/graphics/src/repr/volume/isosurface.ts#L188) builds the id field with the
  floodfilled `scalarField.space`, so the serial ids are read through `floodfilledGet`: ids below `isoLevel` in filled
  cells become `2 * isoLevel`. Build the id field from the unfilled (possibly wrapped) space.
- `getCoords` of rank ≥ 4 spaces also writes `coords[ao[d + 1]]`, an index past the axis order (a stray `undefined`
  key; the coordinates themselves are right), [tensor.ts](../../packages/core/src/math/linear-algebra/tensor.ts)
  `function getCoords`.
