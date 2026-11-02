# Plan: Primitive expansion in the renderer, not the geometry

## Status

- **Current**: every expanded kind is per-item data now: spheres, points, cylinders, lines and text (expanded in the
  vertex shader; levels of detail for spheres and points), with explicit counts (step 5). Open: shared index blocks
  (step 6). Most of the WebGPU section came with other plans (see there).
- Related: [wgsl-strategy.md](wgsl-strategy.md) (points as quads, step 4), [storage-buffers.md](storage-buffers.md)
  (data arrays, `DataSpec`), [attribute-bindings.md](attribute-bindings.md) (vertex buffer limits),
  [sphere-impostors.md](sphere-impostors.md) (indexed quads), [geometry-residency.md](geometry-residency.md) (CPU
  copies), [backend-split.md](backend-split.md), [minimum-capabilities.md](minimum-capabilities.md).
- Done for spheres and points:
  - one array of center xyz and group per item, `positionGroup` (`geo/geometry/center-group.ts`), the geometry's only
    copy and the shaders' `tPositionGroup` (`DataSpec('rgba32f')`): 16 instead of 32 CPU bytes per item
    ([geometry-residency.md](geometry-residency.md), item 9); bounding sphere, group mapping and the position iterator
    from it;
  - the shader chunk `assign_quad_position_group` (GLSL and WGSL: quad corner from the vertex id `% 6`, item from
    `/ 6`, position and group from `tPositionGroup`);
  - drawn indexed by one shared index buffer per context (`gl/quad-indices.ts`, item `i` the indices
    `6i + {0, 1, 2, 1, 4, 2}`): four vertex shader runs per item instead of six, the vertex ids and draw counts as
    unindexed;
  - levels of detail shared (`geo/geometry/lod-levels.ts`: params, presets, `orderLodLevels` reordering the items in
    place; the vertex shaders' `assign_lod_factor`, GLSL and WGSL), merged too.

## Question

Builders currently decide how a primitive is drawn: cylinders, lines and text are expanded into per-vertex data (corner
mappings, duplicated item data, index buffers) when the geometry is built. Spheres and points keep one entry per item
and expand in the vertex shader. Expansion should be a renderer detail for all impostor-like primitives: the geometry
holds per-item data once, in the layout the shaders read.

## Current state

| Kind      | Builder output per item                                                                                       | Draw                                     | Bytes per item (CPU, the same again on the GPU) |
| --------- | ------------------------------------------------------------------------------------------------------------- | ---------------------------------------- | ----------------------------------------------- |
| spheres   | center, group (one texel of `tPositionGroup`)                                                                 | 6 vertex ids, shared quad indices        | 16                                              |
| points    | center, group (one texel of `tPositionGroup`)                                                                 | 6 vertex ids, shared quad indices        | 16                                              |
| cylinders | start, group, end, radius scale, cap, color mode (three texels of `tCylinder`)                                | 6 vertex ids, shared box indices         | 48 (was about 360)                              |
| lines     | start, group, end, whether it continues the line before (two texels of `tLine`)                                | 4 vertex ids, shared line indices        | 32 (was about 170)                              |
| text      | position, group, depth, tether flag, corners and texture coordinates (four texels of `tText`)                  | 4 vertex ids, shared line indices        | 64 (was about 170)                              |

- Cylinders are bonds in ball-and-stick (the default for ligands, common for whole structures): they were the largest
  of the expanded kinds by far, and their 7 per-vertex attributes made WebGPU pack vertex buffers (over
  `maxVertexBuffers`); per-item data since step 2, lines since step 3, text since step 4: no schema exceeds
  `maxVertexBuffers` any more.
- Counts (step 5): `uVertexCount` is the data vertices (`Geometry.getVertexCount`: the values of 'vertex' color,
  size and layers per instance, their vertex-instance stride; one per sphere, point and text item, a start and an end
  per cylinder and line, a mesh's vertices), `drawCount` the drawn indices; the item count follows from `drawCount` and
  the topology (indices per item), e.g. the merged item data (`itemTexels`); merged segments carry `drawBases` for draw
  ranges and `vertexBases` in data vertices.
- Instances (structure transforms) are instanced attributes (`aTransform`, `aInstance`), so per-item data can't use
  instancing too.

## Proposal

1. **Geometry holds per-item data only**, once, in the layout the shaders read (as `positionGroup`). Cylinders: start,
   end, radius scale, caps, color mode, group per cylinder; lines: start, end, group per line; text: per character
   position, group, glyph and offsets. No mappings, no duplicated vertices, no index buffers in the builders.
2. **Item data as `DataSpec`s** (data pages on WebGL, storage pages on WebGPU), read by the vertex id: WebGL1 emulated,
   WebGL2 `gl_VertexID`, WebGPU `vertex_index`.
3. **Shared expansion chunks** per topology (as `assign_quad_position_group`): one quad expansion for spheres, points,
   lines and text, one box expansion for cylinders; shared index patterns per topology (as the quad indices) instead of
   per-geometry index buffers.
4. **Counts made explicit**: item count, data vertex count and drawn vertex count, instead of overloading
   `uVertexCount`/`drawCount`; merged renderables concatenate per-item data and derive draw ranges from the topology.

Expected gains: cylinders about 9× less memory and upload (about 40 instead of 360 bytes), under the WebGPU vertex
buffer limit without packing; simpler merging (per-item data concatenates, no index rewriting); levels of detail for
every kind from `lod-levels.ts`; one expansion path per topology on all backends.

## WebGPU beyond storage buffers

Done by other plans, as far as they apply:

- **GPU-driven culling**: merged items culled by a compute pass that writes indirect draw arguments, frustum, levels and
  Hi-Z occlusion (`merged-culler.ts`; [indirect-draws.md](indirect-draws.md),
  [merged-render-objects.md](merged-render-objects.md)).
- **Render bundles**: draw lists recorded and replayed ([webgpu-draw-submission.md](webgpu-draw-submission.md)).
- **Per-object data in one buffer**: object uniforms in pages at dynamic offsets
  ([uniform-buffers.md](uniform-buffers.md)).
- **Base instance in core**: `firstInstance` in the draws, no base instance extension paths on WebGPU.
- **Compute shaders** for GPU surfaces ([webgpu-compute.md](webgpu-compute.md)).

Open: **quads without index buffers** (4-vertex triangle strips per item with the item from `instance_index`) needs
the instance transforms out of instanced vertex attributes; the shared quad indices took the vertex shader runs to four
per item already, so it pays only with that change for other reasons (see Transforms from storage below).

## Shared index blocks (16-bit, base vertex)

The shared quad and box indices (`quad-indices.ts`) hold absolute vertex ids (`6i + pattern`), so they cover the largest
item drawn, in 32 bits past about 10 k items: 24 bytes a quad, 48 a box. 3J3Q as spacefill (2.44 M spheres) keeps about
59 MB of quad indices, 4V6X's bonds (650 k cylinders) 31 MB of box indices; a PetWorld-sized scene reaches hundreds of
MB, and quads past the buffer limit (`maxQuadItemCount`) draw unindexed (six vertex shader runs per item, not four).

- **Blocks**: one fixed `Uint16` buffer per pattern for a block of K items (K = 10922: `6K - 1` within 16 bits; 64 KB
  for quads, 128 KB for boxes). A draw of a range of items splits per block: `firstIndex` the item's position in its
  block times the indices per item, `baseVertex` the block's first vertex id (`6 × blockStart`), so the vertex ids the
  shaders see stay absolute and the shaders don't change. A range across blocks becomes a draw per block (650 k
  cylinders: about 60 draws instead of one; replayed in render bundles).
- **WebGPU**: `baseVertex` is core (`drawIndexed`, and the indexed indirect records the merged culler writes). The CPU
  draw lists (`SegmentedMdbList`, the multi-draw-base data) and the culler's records split at block boundaries.
- **WebGL2**: no base vertex in practice (Chrome exposes `WEBGL_draw_instanced_base_vertex_base_instance` only behind
  a flag, as with the merged items' `dInstanceData`); it keeps the 32-bit buffer. A per-draw `uVertexBase` added to
  `gl_VertexID` would work for single draws, not for multi-draw lists (`gl_DrawID` would need an array of bases).
- **Gain**: GPU memory only (the context's buffers constant, below 0.2 MB, instead of 24 or 48 bytes per item of the
  largest draw), and no unindexed fallback past the limit; frame time unchanged (the post-transform reuse is within an
  item either way, 16-bit fetches marginally cheaper).

## Transforms from storage (instance strips, WebGPU)

Drawing each item as a 4-vertex strip with the item from `instance_index` needs no index buffer at all, but the
structure transforms use instancing today (`aTransform`, `aInstance` instanced attributes); they would move to storage,
read per vertex by index.

- **Not slower by itself**, measured on its closest analogue: merged items on WebGL2 without base instances read the
  transform and the instance and segment data per vertex from a texture (`dInstanceData`). Read once per vertex, they
  are at parity with or faster than instanced attributes (AMD iGPU, opaque pass: 4V6X merged 19.5 ms against 20.5 ms
  for attributes re-pointed per draw and 19.7 ms unmerged; 1000 cartoon units 7.81 against 8.19 ms). Re-fetched per use
  with an integer division per fetch they were 17-23% slower: the cost is repeated fetches and divisions, not the
  storage read. An instanced attribute is fetched per vertex as well; the vertices of an item read the same transform,
  which stays in cache.
- **Where it could cost**: one `instance_index` carrying both item and transform (`item × transformCount + t`) needs a
  division per vertex (the AMD cost above); the alternatives are a draw per transform (assemblies of 60 copies: 60× the
  draws, cheap in render bundles) or a power-of-two layout (shifts). The merged culler compacts the instanced
  attributes of visible cells today and would write transform index lists instead.
- **Gain**: the same as the blocks (no index memory), no frame time beyond it (four vertex runs per item either way).
  So the blocks first; this only with other reasons to move the transforms (e.g. the culler's compaction). A check
  before it: a spheres variant reading its transform from a storage page by `instance_index`, timed against the
  instanced attributes on both GPUs.

## Steps

1. Spheres and points (done, see Status).
2. **Cylinders** (done): one array per geometry, `Cylinders.cylinderData` (the shaders' `tCylinder`,
   `DataSpec('rgba32f')`), 3 texels per cylinder: start xyz and group, end xyz and radius scale, cap and color mode
   (48 instead of about 360 bytes). The vertex shaders take the cylinder from the vertex id (`/ 6`) and the box corner
   from `% 6` (the mapping a constant table); the box indices are a shared pattern per context like the quads'
   (`quad-indices.ts`: `ItemIndexPattern`, `BoxItems`, 12 indices per item; `ShaderCode.itemIndices`), always indexed
   (boxes have no unindexed form, past `maxBoxItemCount` not drawn). Merged: `tCylinder` an extra data array of the
   descriptor, no elements to rewrite. The counts stay as they were (6 vertex ids, 12 indices per cylinder), so the
   per-vertex theme layers and draw ranges are unchanged; explicit counts are left to step 5. On WebGL2 the three texels
   are read by one division (`read_<name>3`, the accessors use `texelFetch`).
   - Measured (4V6X ball-and-stick bonds, about 650 k cylinders, before → after): CPU geometry 235.7 → 31.4 MB; GPU
     attributes 204 → 0 MB, elements 31.4 MB both (now the shared box indices), data +31 MB; opaque pass NVIDIA
     WebGL2 0.52 → 0.40 ms, WebGPU 0.54 → 0.41 ms, AMD iGPU WebGL2 22.2 → 22.6 ms, WebGPU 29.5 → 14.4 ms; building
     the unmerged representation 257 → 158 ms. Images unchanged on all profiles (the bounding sphere computed as from
     6 vertices per cylinder, as before).
   - Left: levels of detail for cylinders (`lod-levels.ts`); the shared box indices are as large as the largest draw
     (31 MB here): 16-bit patterns per item range would quarter them.
3. **Lines** (done): one array per geometry, `Lines.lineData` (the shaders' `tLine`), 2 texels per line: start xyz
   and group, end xyz and whether the line continues the one before (32 instead of about 170 bytes). Strips
   (`StripLinesBuilder`) are lines too: a line continuing the one before takes that line's end corners as its start
   corners (its start, across the line before, as the shared vertices of a strip did; one more texel read where it
   does), so `stripCount` and `stripOffsets` are gone and the strips' joints are drawn as before (a visual case,
   `geometry/lines/strips`, with a fixed camera). The vertex ids stay 4 per line (`LineItems`, 6 indices over them,
   shared per context; `uVertexCount` and the per-vertex layers as they were); WebGL2 reads the two texels by one
   division (`read_<name>2`).
   - Measured (4V6X as lines, before → after): CPU geometry 110 → 21 MB; GPU attributes 94 → 0 MB, unmerged
     elements 17.7 → 8.2 MB, data +21 MB; opaque pass NVIDIA WebGL2 0.32 → 0.29 ms (merged), 0.39 → 0.36 ms
     (unmerged), WebGPU 0.35 → 0.32 and 0.38 → 0.31 ms; AMD iGPU WebGPU merged 31.2 → 15.3 ms, unmerged 15.3 → 15.6
     ms; AMD iGPU WebGL2 17.4 → 18.8 ms and 18.0 → 18.9 ms (5-8% slower: the data texture reads about 3-4%, the
     strips' branch 2-4%, measured without it; a define for segment-only geometries would take the branch out, but
     has to follow a representation switching builders, e.g. streamlines). Images unchanged on all profiles.
   - The strips' bounding sphere was read over twice the buffers' length (2 vertices per point counted as 4): now that
     of their points.
4. **Text** (done): one array per geometry, `Text.textData` (the shaders' `tText`), 4 texels per item (a glyph,
   a background or a tether): the label's position and group; its depth and whether the item is a tether; then a
   glyph's or background's corners and texture coordinates as rectangles (left, top, right, bottom; the background's
   texture coordinates 10), or a tether's four corners (64 instead of about 170 bytes). No glyph table: a glyph's
   texture rectangle per item is as compact and keeps the atlas lookup in the builder. Drawn by the shared line
   indices (4 vertex ids, `[0, 1, 2, 1, 3, 2]`): a tether's corners ordered base B, tip, base center, base A, so the
   pattern draws its two triangles (tip, base A, base center) and (tip, base center, base B) with their winding as
   before. Images unchanged (`geometry/text/labels`, and a new `geometry/text/tether`, identical to before).
5. **Counts everywhere** (done): `uVertexCount` the data vertices for every kind (before: spheres 6 per sphere,
   though their shaders read one; cylinders 6 and lines 4, though corners only took a start or an end; text 4, all the
   same), the shaders' `vertexId` (what the color, size and layer chunks index 'vertex' data by) the data vertex (a
   cylinder's and a line's corners at the start read the start, the others the end, as the position iterators' strides
   gave them before), the position iterators one location per data vertex. Merged objects size 'vertex' data by
   `uVertexCount` and the impostors' item data by their own granularity, 'item' (`itemTexels`, from `drawCount`): before
   both shared one count, so merged cylinders and lines colored per vertex read their members' colors at the wrong
   offsets (unit test), and spheres colored per vertex and instance read past their colors from the second instance.
   'vertex' color, size and layers of spheres take a sixth, of cylinders a third, of lines a half, of text a quarter.
   Geo-export reads colors and sizes at the data vertices. Images unchanged.
6. **Shared index blocks** on WebGPU (above): a fixed 16-bit block per pattern, draws split per block with base vertex
   (CPU draw lists, the culler's indirect records); WebGL2 keeps the 32-bit buffer. Tests: the split ranges draw what the
   unsplit ones did (images unchanged, draw lists and GPU culling included), index bytes constant; measured on 3J3Q
   spacefill and 4V6X bonds.
7. **Transforms from storage** (above), only with another reason to move them; the storage-read check first.

## Open questions

- Geo-export and other consumers of per-vertex values (cylinders `aStart`/`aEnd`, lines, text) move to the per-item
  arrays with each step.
