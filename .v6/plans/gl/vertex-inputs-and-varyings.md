# Plan: Vertex inputs and varyings

## Status

- **Current**: steps 1 to 6 done (the baseline, the transform as 3 vec4, segment data by member in part,
  per-primitive varyings flat and packed, layer scalars packed, mesh normals oct-encoded), see Measurements.
- **Proposed**: an analysis of what the vertex stages read (vertex buffers, per-instance attributes, data arrays) and
  pass to the fragment stages (varyings, inter-stage variables), with the reductions worth measuring. Follows
  [primitive-expansion.md](primitive-expansion.md) (per-item data, done for every impostor kind) and the instance data
  of [merged-render-objects.md](merged-render-objects.md) (`dInstanceData`, 7 texels an instance); the transform as 3
  texels was held back from there for this plan.

## Question

Per vertex the shaders read per-instance attributes (a transform, an instance id, merged objects' segment data) and,
for meshes, per-vertex attributes; they pass positions, colors, theme layers and per-primitive constants to the
fragment stage. Which of these are larger or more numerous than needed, and which reductions pay in memory, bandwidth
and interpolation?

## Current state

Per instance (vertex buffers, instance step; WebGL merged objects without base instances read the same from the
instance data texture, `dInstanceData`):

| attribute        | type      | bytes | where                                                       |
| ---------------- | --------- | ----- | ----------------------------------------------------------- |
| `aTransform`     | 4 × vec4  | 64    | every kind                                                  |
| `aInstance`      | f32       | 4     | every kind: the instance id (color, marker, layer indices)  |
| `aSegment`       | 3 × f32   | 12    | merged (`dSegmented`): group count, group offsets per member |
| `aSegmentSphere` | vec4      | 16    | merged: the member's invariant bounding sphere              |
| `aSegmentLod`    | vec4      | 16    | merged: the member's level scales                           |

68 bytes an instance, 112 merged; the GPU culler (`InstanceCopySchema`, `cell-cull.ts`) writes all 112 per visible
instance per frame into the compacted buffers (`transforms`, `instances`, `segments`, `spheres`, `lods`). The segment
data is the same for every instance of a member.

Per vertex: meshes `aPosition` (vec3), `aNormal` (vec3), `aGroup` (f32), 28 bytes; image and direct-volume their own;
spheres, points, cylinders, lines and text no vertex buffers (their data arrays, 16 to 64 bytes an item).

Varyings (`vertex_output_wgsl`, the GLSL counterparts): `vModelPosition`, `vViewPosition` (vec3 each) always declared;
`vColor` (vec4, the alpha never read: materials take `uAlpha`); one location each for `vMarker`, `vTransparency`,
`vEmissive`, `vPaletteV`, `vClipping` (scalars; marker and clipping `flat` with `dVaryingGroup`), `vOverpaint` and
`vSubstance` (vec4), `vPickIds` (flat vec2, pick); per kind:

| kind      | own varyings                              | read by the fragment stage                                      |
| --------- | ----------------------------------------- | --------------------------------------------------------------- |
| mesh      | `vNormal`                                 | model and view position, color, normal                          |
| spheres   | `vRadius`, `vPoint`, `vPointViewPosition` | radius, point, center, color (model position: `dApproximate`)   |
| cylinders | `vStart`, `vEnd`, `vSize`, `vCap`         | model position, start, end, size, cap, color                    |
| points    | `vPointCoord`                             | point coordinate, color                                         |
| text      | `vTexCoord`                               | texture coordinate, color                                       |

- Outputs the fragment stage doesn't read cost nothing: Dawn compiles the vertex stage with the fragment's inputs only
  (the spheres' HLSL signatures have neither model nor view position), GL linkers drop unread varyings.
- Per-primitive constants are interpolated: the spheres' radius and center, the cylinders' start, end, size and cap are
  the same at every corner of their quad or box, yet smooth and one location each (spheres 3, cylinders 4 locations).
- Every theme layer scalar takes a location of its own: with overpaint, transparency, emissive, substance, marker and
  clipping a mesh passes 9 locations besides positions and normal; WebGPU guarantees 16
  (`maxInterStageShaderVariables`).

## Measurements

Baseline (step 1): `pnpm benchmark --gpu-timing`, Chrome, 1280x720 unless noted; GPU bytes by `getByteCounts` after
the setup (`byteNotes`, now also on the instance, mesoscale and Gaussian surface cases, storage included).

| case                                                   | GPU bytes (WebGL2 / WebGPU)                                   |
| ------------------------------------------------------ | ------------------------------------------------------------- |
| `merging/instances-1000-spacefill-merged-rotating`     | attributes 10.8 MB (100k instances, 112 bytes each) / same    |
| `merging/instances-1000-spacefill-merged-gpu-rotating` | — / storage 33.8 MB (the instance data and its compacted copies) |
| `merging/postsyn-lod-pixel3-zoom-nohiz`                | attributes 55.4, textures 146 MB / storage 338 MB              |
| `merging/4v6x-gaussian-merged`                         | attributes 43.0 MB (28 bytes a vertex) / same                  |

Opaque pass (median, ms), NVIDIA RTX 4070 Ti / AMD integrated Radeon (RDNA2, 2 CUs), WebGL2 and WebGPU:

| case                                                   | NVIDIA WebGL2 | NVIDIA WebGPU | AMD WebGL2 | AMD WebGPU |
| ------------------------------------------------------ | ------------- | ------------- | ---------- | ---------- |
| `passes/spheres`                                       | 0.43          | 0.46          | 27.5       | 33.7       |
| `passes/cylinders`                                     | 0.15          | 0.11          | 7.7        | 8.5        |
| `passes/mesh`                                          | 0.11          | 0.08          | 3.7        | 4.1        |
| `cylinders/4v6x-bonds-merged`                          | 0.32          | 0.31          | 20.0       | 15.8       |
| `merging/instances-1000-spacefill-merged-rotating`     | 0.58          | 0.39          | 20.0       | 20.3       |
| `merging/instances-1000-spacefill-merged-gpu-rotating` | 0.68          | 0.56          | 20.0       | 20.1       |
| `merging/postsyn-lod-pixel3-zoom-nohiz`                | 5.07          | 4.84          | 362        | 342        |
| `merging/4v6x-gaussian-merged`                         | 0.29          | 0.39          | 15.9       | 16.3       |

- Per-instance data dominates the mesoscale scene: about 500k instances, 55 MB of instance data (112 bytes each); on
  WebGPU 338 MB of storage, about six times that (the instance data and the culler's compacted copies, sized per lod
  level; the 1000-instance case, one level, three times). Steps 2 and 3 halve all of it. On WebGL2 it is the
  attributes plus the instance data texture.
- The NVIDIA passes are well under a millisecond: varying changes (steps 4, 5) can only show on the iGPU, which is bound
  by the fragment stage (spheres, the mesoscale scene).
- Mesh vertices: 43 MB for 4V6X's Gaussian surfaces, 30.7 MB with oct-encoded normals (step 6).
- Also seen: 116 MB of elements in the mesoscale scene, the shared quad indices
  ([primitive-expansion.md](primitive-expansion.md), step 6: 16-bit blocks).

## Proposal

Measure first (step 1), then each change gated on its own measurement; images within the usual tolerance.

1. **Measurements** (done, see Measurements): GPU timing (`pnpm benchmark --gpu-timing`, NVIDIA and the AMD iGPU, both backends) of a
   fragment-bound impostor pass (`passes/spheres`, `cylinders/4v6x-bonds-*`), many instances
   (`merging/instances-1000-spacefill-merged-rotating`, with and without GPU culling), a mesoscale scene zoomed
   (`merging/postsyn-lod-pixel3-zoom-nohiz`), and a large surface (`merging/4v6x-gaussian-*`); memory by
   `getByteCounts`. The baseline for every step.
2. **Transform as 3 vec4** (done): attributes of item size 12 are affine transforms, on the GPU their first three rows
   (`affine-rows.ts`: `packAffineRows`, `aTransform0` to `aTransform2`; WGSL `affineTransform`, GLSL `affineFromRows`
   and an `aTransform` macro), on the CPU still a `Mat4` per instance: `TransformData` and its readers (bounding
   spheres, instance grids, position iterators, culling, geo-export) keep their layout, and the render items take the
   rows from a derived cell (`affineRowsCell`, in `splitValues`) that repacks when the transforms change (rarely: tumble
   is in the shaders). The instance data 6 texels (was 7, `InstanceDataVec4s`), the culler copies 3 vec4, direct
   volumes pass the rows as 3 flat varyings (was 4). GPU bytes: 100k merged instances 10.8 → 9.3 MB, GPU-culled 33.8 →
   29.2 MB of storage; the mesoscale scene WebGL2 attributes 55.4 → 47.5 MB and textures 146 → 138 MB, WebGPU storage
   338 → 291 MB. Opaque pass on the AMD iGPU, WebGL2: the mesoscale scene 362 → 348 ms, 100k instances 20.0 → 19.0,
   Gaussian surfaces 15.9 → 14.8; WebGPU within 1 % (Gaussian unmerged 16.3 → 15.0); NVIDIA unchanged (runs of the
   same build vary by up to 60 % under a millisecond). Images unchanged on every profile.
3. **Segment data by member** (done, in part): per instance `aSegment` (vec4: the segment's group count and offsets,
   read by every color, marker and layer index, and the segment's index in `w`); per segment its invariant bounding
   sphere and level scales, read only by tumble, instance clipping and the level scales, in a data array (`tSegment`,
   two elements a segment: storage on WebGPU, a data texture on WebGL, the shaders' `aSegmentSphere`, `aSegmentLod`). 68
   instead of 112 bytes an instance (with step 2), the instance data 5 texels, the culler's copy 68 bytes and 8 storage
   buffers (was 10): culling on the GPU now fits the core limits (`webgpu-core-limits.test.ts`). Tried first, the groups
   in the data array too (56 bytes an instance): a dependent read per vertex where every index needs it, on the AMD iGPU
   the merged cases' opaque pass 8–12 % slower (4V6X's merged bonds 19.4 → 21.0 ms on WebGL2, 15.6 → 17.2 ms on WebGPU),
   NVIDIA unchanged; with the groups per instance as they are now, the common paths are not slower. The paths that read
   the segment's row (on-demand cases `merging/instances-1000-spacefill-merged-{tumble,clip}-rotating`, culled on the
   GPU too, and the mesoscale scene's levels of detail), AMD iGPU WebGL2: tumble 46.4 → 47.5 ms (+2 %), instance
   clipping 24.3 → 25.9 ms (+7 %), the mesoscale scene 344 → 355 ms (+3 %); WebGPU within 1 %, NVIDIA within noise. GPU
   bytes: 100k merged instances 9.3 → 6.5 MB, GPU-culled 29.2 → 22.3 MB of storage; the mesoscale scene WebGPU storage
   291 → 214 MB (338 before step 2), WebGL2 attributes 47.5 → 33.6 MB and textures 138 → 130 MB; CPU merged copies 55.2
   → 41.4 MB. Opaque pass on the AMD iGPU against step 2 within 1 % (WebGL2 slightly faster). Found on the way: since
   step 2 the culler wrote the identity region (views not culled on the GPU) with the transforms' 16 floats where it
   holds 12: now packed as rows (`webgpu-gpu-culling.test.ts` draws that region against the CPU culling's image). Images
   unchanged on every profile.
4. **Per-primitive varyings flat and packed** (done): spheres pass their view-space center (negated) and radius as one
   flat vec4 (`vSphere`, was `vPointViewPosition` and `vRadius`), cylinders their start with size and end with cap as
   two (`vCylinderStart`, `vCylinderEnd`, were four); the point stays interpolated. Flat on WebGL2 and WebGPU,
   interpolated on WebGL1 (no flat there), packed on all. Cylinders' color variants no longer exceed WebGPU's 16
   inter-stage variables (the known exceedances of `minimum-limits.test.ts` are none). Images: 1/255 on a few pixels
   (the constants exact, interpolated they rounded), a few silhouette pixels of spheres on WebGL2 (up to 23 of
   `geometry/spheres/grid`), WebGL1 identical; the suite passes on every profile. Opaque pass on the AMD iGPU, 1–6 %
   faster: spheres 27.5 → 27.0 ms (WebGL2), 33.7 → 33.3 (WebGPU), cylinders 8.47 → 7.95 (WebGPU), 4V6X's bonds 14.66 →
   14.44 (WebGPU, three runs each), the mesoscale scene 360 → 355 (WebGL2), 341 → 337 (WebGPU); NVIDIA within noise.
5. **Layer scalars packed** (done, WGSL): the palette coordinate in the color's alpha (unread: materials take `uAlpha`),
   the marker with the clipping flags (`vMarkerClipping`, both flat by group), the transparency with the emissive
   (`vTransparencyEmissive`): with every layer, 2–3 locations fewer (spheres and cylinders 12 → 9 in color variants, 10
   in pick; meshes, points and text 11 → 8; lines 10 → 7), 6 or more under the 16 inter-stage variables. GLSL unchanged:
   its compilers pack varyings into vec4 rows (GLSL ES 3.00 packing, `getPackedVaryingRows`); WebGPU counts a location
   per varying whatever its size. `minimum-limits.test.ts` now counts the WGSL programs' output locations (it took the
   GLSL varyings for them). WebGPU images identical to before (the GPU smooth-color surface differs between runs of the
   same build). Timing not measured: no benchmark case draws theme layers, and the marker alone keeps its location.
6. **Mesh normals oct-encoded** (done): `aNormal` an `oct16` attribute (`oct-normals.ts`), `snorm16x2` on WebGPU and
   normalized `SHORT` × 2 on WebGL1 and 2, 28 → 20 bytes a vertex; the shaders decode (`octDecode`, branchless, the
   callers normalize). On the CPU the normals stay three floats (`Mesh.normalBuffer`: the builders, picking and
   geo-export unchanged); the render items take them from a derived cell (`octNormalsCell`, in `splitValues`) that
   encodes when they change, merged meshes part by part (a segmented array of `Int16Array` parts, only changed parts
   encoded and uploaded again). GPU-built meshes (`texture-mesh`) keep theirs. GPU bytes: 4V6X's Gaussian surfaces
   attributes 43.0 → 30.7 MB, 3J3Q's 7.0 → 5.0 MB. Opaque pass: NVIDIA within noise or faster (3J3Q WebGPU 0.24 →
   0.19 ms); the AMD iGPU at 1280x720, 4V6X's Gaussian surfaces (three runs each), about 2 % slower (vertex bound
   there: the decode): WebGL2 merged 14.75 → 15.10 ms, unmerged 15.0 → 15.4, WebGPU merged 16.1 → 16.4, unmerged 14.9 →
   15.2 (a decode with a branch and normalize 4 % on WebGPU); at 2560x1440 merged about 9 % faster (one run);
   `passes/mesh` unchanged. Images: about half differ slightly, most within 6/255, at most 21 on 76 pixels
   (`renderer.illumination.multisample`); the suite passes on every profile.

Not proposed: f16 varyings ([shader-f16.md](shader-f16.md): slower on NVIDIA), the view position derived from the model
position in the fragment stage (a matrix product per fragment for one location), removing unread outputs (stripped
already).

## Open questions

- Whether the number of interpolated locations shows in the timings of desktop GPUs at all, or only on the iGPU
  (step 4 and 5 hinge on it).
- Shape transforms that are not affine: none known; the assertion of step 2 finds them.
- The dependent read of step 3 against the attribute fetch it replaces, per GPU.
