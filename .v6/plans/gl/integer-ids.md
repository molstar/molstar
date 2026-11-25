# Plan: integer ids and indices

## Status

- **In progress**: steps 1 (integer instance ids and data indices), 2 (32-bit pick ids) and 3 (the yeast cell in the
  Viewer: marking and picking instances past 2^24 checked by hand) done; step 4 (group and object ids in 32 bits) done:
  the vertex stages, picking, quantized sites and every geometry's group storage as integers (spheres, points, cylinders,
  lines, text, images; volume dots past 2^24). Step 5 (isosurface groups past 2^24) done: the GPU paths (WebGL2, WebGPU) and meshes' `aGroup`. Found on the yeast cell ([particle-pack.md](../particle-pack.md) step 6): marking instances past 2^24 hits
  their neighbours.

## Why

- **Floats carry ids**, exact up to 2^24 (16,777,216): the instance id `aInstance` (a `float32` attribute; `CulledData`'s
  `aInstance`, `culled-data.ts`; WebGL's instance data texture `tInstanceData`), the segment data `aSegment` (group
  count, groupInstance offset, group offset, segment index: `float32` x 4), the index arithmetic of the data layers
  (GLSL in `float`: `granularityIndex`, `groupInstanceIndex`, `groupIndex`; WGSL in `f32`, then `i32`) and picking
  (`vPickIds`, `vec2<f32>`; the instance id packed into 24 bits, the RGB of an `rgba8` target).
- **The yeast cell's merged objects hold 31.5 M instances** (its sites 20.3 M groups): past 2^24 an instance's marker
  index rounds (marking highlights or selects a neighbour), a member's group index (`aSegment.z + group`) rounds (its
  per-group data, colors and sizes by group, read the neighbours'), and picking cannot name the instance at all.
- WebGL1 has no integer attributes: it keeps floats, and the limit (the yeast cell's scale is out of its reach anyway).

## Design

- **Integer attributes**: `AttributeKind` gains `uint32` and `int32`; WebGL2 binds them with `vertexAttribIPointer`
  (`in uint`, `in ivec4`), WebGPU as `uint32` and `sint32x4`; WebGL1 as floats (the same arrays converted where bound).
  Why: an attribute is the only per-instance input of the vertex stage besides storage; floats cannot hold the ids.
- **`aInstance` `uint32`** (`Uint32Array`: the transform data, merged objects' own, the instance storage),
  **`aSegment` `int32` x 4** (a groupInstance offset is negative: the member's base minus its instance base times its
  group count).
- **Index arithmetic in integers**: an `int` (`i32`) instance id and segment in the vertex stage, the data accessors'
  integer overloads (`read_<name>(int)`); GLSL macros of one meaning on both versions (`instanceId`, `segmentData`:
  `int(...)` of the float attributes on WebGL1). A groupInstance index (`offset + instance * groupCount + group`) stays
  within `int32` as the layers it indexes do (a layer of more elements cannot be bound).
- **Ids in float buffers as bits**: `CulledData`'s `aInstance` and segment indices (`u32` bits, `bitcast` in WGSL),
  `tInstanceData` (`floatBitsToUint` and `floatBitsToInt` on WebGL2). Why: one buffer of `vec4<f32>` for the culler's
  data (one binding), the bits exact.
- **Picking (step 2)**: the ids as flat integer varyings (`flat out uvec2`, `@interpolate(flat) vec2<u32>`); the
  instance id 32 bits, the RGB and the alpha byte of the `rgba8` target (as now, no new target format); the null id the
  largest value; `PickPass` decodes four bytes.

## Steps

1. **Integer instance ids and data indices** (done): the attribute kinds and their bindings (WebGL2
   `vertexAttribIPointer`, WebGPU `uint32` and `sint32x4`, WebGL1 converted to floats on upload), `aInstance` and
   `aSegment` (the transform data, the instance grid, merged objects), `CulledData` and `tInstanceData` as bits, the
   index arithmetic of every geometry's vertex stage (GLSL, WGSL) and of the raster in integers. Tested by
   `geometry/spheres/instance-ids` (two instances of ids 2^24 + 1 and 2^24 + 3 colored from a layer that long: each its
   own color; it fails on floats): the visual harness now applies `checkImage` to image passes too (it did to canvas
   captures only). Images unchanged (the full suite: WebGL2 and WebGL1 against the commit before, WebGPU against WebGL2).
2. **32-bit pick ids** (done): the instance id packed from the integer (`packPickId`: id + 1, its 24 low bits in rgb, the
   top byte in alpha; volumes and images, few instances, from their floats with alpha 0), WebGPU's pick ids a flat
   `vec2<u32>`; `PickPass` reads four bytes, the white clear of the instance target `PickingId.InstanceNull` (2^32 - 2).
   Tested by `geometry/spheres/instance-ids` picking across the canvas: the ids 2^24 + 1 and 2^24 + 3 (none read on 24
   bits).
3. **Measured** (done): the yeast cell's frame and heap unchanged (an attribute of 4 bytes as before), marking and picking
   instances past 2^24 in the Viewer.
4. **Group and object ids in 32 bits** (built in part, see As built). Steps 1 and 2 made instances exact; a geometry's own group is
   still a float from its source to the pick target, so groups past 2^24 within one geometry (the yeast cell's sites:
   20.3 M) round: their marking, their data by group (`groupIndex = segment.z + group`, the float rounded before the
   integer add) and their picking hit neighbours, and picking writes the group in 24 bits (`packIntToRGB`, alpha 1).
   The object id the same 24 bits (`uObjectId`; never near 2^24 objects, aligned for one encoding).
   - **The group as an integer from its source**: `aGroup` `uint32` (mesh, image, and the geometries with it; WebGL2
     `vertexAttribIPointer`, WebGPU `uint32`, WebGL1 floats as step 1); spheres' and points' `tPositionGroup` w as
     bits (`floatBitsToUint`, `bitcast`: one `rgba32f` texel, the center exact as now); quantized sites' group
     from the item index (integers already); texture meshes' `tGroup` (RGB of a float texture, a GPU-computed surface:
     24 bits there, its own limit, noted).
   - **The vertex stage in integers**: `group` an `int` (`assign_group` both backends, the chunks reading it:
     color, size, marker, clipping, transparency, emissive, substance, the raster and its resolve).
   - **Picking**: the group as the instance (`packPickId`: id + 1, its 24 low bits in rgb, the top byte in alpha), the
     object too; `PickPass` reads four bytes for both (`getInstanceId` for all three). Nulls: the clear (white) reads
     2^32 - 2 as for instances; `PickingId.Null` (2^24 - 2) is the no-group of about twenty representations and the
     object's null: a `GroupNull` 2^32 - 2 for groups (those representations' checks), the object's kept or moved alike.
   - Tests: a geometry of more than 2^24 groups colored by group and picked at 2^24 + 1 and 2^24 + 3 (as
     `geometry/spheres/instance-ids` for instances: it fails on floats); the visual suite on all profiles (images
     unchanged); the yeast cell's sites marked and picked past 2^24 in the Viewer.
   - **As built**: `group` an integer in every geometry's vertex stage (GLSL `GroupId`, `int` on WebGL2 and `float`
     on WebGL1; WGSL `i32`), each source converted where read (`aGroup`, `tPositionGroup`'s w, cylinders', lines'
     and text's w, texture meshes' `tGroup`, marching cubes' packed group), the wiggle reading it as a float; quantized
     sites' group from the item index exact (the yeast cell's spheres, impostors and the raster, its `sphereGroup` a
     `u32` beside the center: an integer passed as float bits reads a denormal, which the GPU flushed: every pixel of
     the raster test wrong). Picking: object, instance and group by `packPickId` (GLSL overloads for `float`, `int`,
     `uint`), volumes and images too; `PickPass` decodes four bytes for all three; one null, `PickingId.Null` 2^32 - 2
     (the white clear; the mesh caps write it), `InstanceNull` folded in (`breaking-v6-changes.md`, Pick ids in 32
     bits). Tested: `pick-buffers.test.ts` (ids past 2^24 and the clear decoded); the visual suite (WebGL2, WebGL1,
     WebGPU) images unchanged; the raster and sites tests on Dawn.
   - **As built, the geometries' group storage**: integers in every geometry, the same memory (the group the fourth
     float of a texel, its bits now). Volume dots (a group per voxel, any volume above about 256³) the one source past
     2^24 today. The data of spheres, points, cylinders, lines and text (`tPositionGroup`, `tCylinder`, `tLine`,
     `tText`) a data format of their own, `rgb32f-i32` (floats but the groups, an item's group the first of its
     `DataSpec.itemTexels` elements' fourth: the CPU arrays `Float32Array`s, the groups written as their bits,
     `int32ToFloatBits`, read through an `Int32Array` view, `getDataWords`): integer bits in a float texel read as
     denormals, which GPUs flush, so the backends read the words as integers: WebGL2 `RGBA32I` pages (`intBitsToFloat`
     for `read_<name>`), WebGPU the storage words; the groups by `readGroup_<name>(i)` (the raster's `sphereGroup` the
     word); WebGL1 float pages, the groups converted where written (`convertGroupsToFloat`, exact to 2^24 as before).
     Images: the group texture's four bytes (id + 1, the top byte in alpha, `packIntToRGBAArray`: the alpha free, the
     plane image's coverage in the color's), read by `unpackGroupId`, their data index an integer (WebGL2, WebGPU).
     One texel for position and group, not two: the hot paths read both, a position alone padded to 16 bytes anyway.
     Tested: visual cases `geometry/spheres/volume-dot-group-ids` and `geometry/points/volume-dot-group-ids` (the
     disc's dots picked past 2^24, odd ones: a float rounds to even); `center-group.test.ts` (groups past 2^24 and -1
     through copies), `number-packing.test.ts` (four bytes); the full visual suite on all profiles unchanged;
     `breaking-v6-changes.md` (Groups as integers in every geometry).
5. **Isosurface groups past 2^24** (done). An isosurface's group is its voxel's index, more than 2^24 for any
   volume above about 256³ (512³: 134 M), and lost on each path:
   - **WebGL, on the GPU** (`isosurface.frag`, texture meshes): the index in float math (`p.z + p.y * dim.z + ...`),
     clamped to white past 2^24 - 1 (which read as the old null, the whole instance; since step 4 as group
     16,777,214), stored in rgb of `tGroup`. Now: the index in integers on WebGL2 (`ivec3`, the dims as `int`),
     its 32 bits as four bytes (rgb and the top byte in alpha, as the pick ids), read by `assign_group` and the color
     smoothing. `tGroup` is `rgba8` on WebGL2 already (4 B a vertex), `rgba32f` on WebGL1 (one bit-plane count for
     all attachments; its index in floats, exact to 2^24 as everything there): no integer texture needed, an `r32ui`
     `tGroup` would save nothing on WebGL2 and WebGL1 has none.
   - **WebGPU, on the GPU** (`marching-cubes-chunks`): the index a `u32`, exact, but the vertex arena's 3 words held
     the 16-bit quantized position, the 24-bit oct normal and 24 bits of the group (`isoPackedGroup`, `packGroup`
     clamped). As built: the group + 1 in the third word whole, the normal octahedral in 16 bits (8 an axis) beside z:
     the arena's layout and size as before. Not a fourth word: the arena is one pool of pages for all surfaces, triangles
     and vertices alike 3 words; a stride per surface splits it. The normal's 16 bits by eye (7QPD's Gaussian surface,
     WebGPU, 1280×720: 4.6 % of the pixels differ, by 0.05 levels on average, single pixels up to 57; the shading the
     same). The CPU readback (`readChunkedVertices`) gives the group's top byte in alpha, the exporter decodes it.
   - **As built (WebGL)**: as above; Gaussian surfaces' packed groups (`dPackedGroup`, 24 bits from the volume) written
     with alpha 0.
   - **On the CPU** (meshes): the index in JS, stored in `aGroup` `float32`. As built: `aGroup` `uint32` for meshes
     (`Mesh`'s group buffer a `Uint32Array`: `Mesh.create`, the mesh and marching cubes builders, the meshes extension,
     the exporter; merged meshes' segmented attributes by their kind; WebGL1 floats as `aInstance`), a breaking change.
     Tested by `geometry/mesh/isosurface-cpu-group-ids` (the GPU case's volume on the CPU: before, no odd group past
     2^24, rounded to even).
   - Tests: `geometry/texture-mesh/isosurface-gpu-group-ids`, a disc in the last slice of a 256 × 256 × 258 grid
     (more than 2^24 voxels; a 2048 × 2048 × 5 grid lays out too large for the GPU isosurface): picked across the canvas,
     every group past 2^24 and some odd, on WebGL2 and WebGPU (before, 16,777,214). The surface faces the empty slices
     and is single-sided: seen from behind it is culled whole (the case blank at first, taken for a limit of 2^24 groups
     before the renderer showed the draws issued). The visual suite unchanged (WebGL2, WebGL1, WebGPU).
