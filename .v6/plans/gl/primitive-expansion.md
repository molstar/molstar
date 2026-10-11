# Plan: Primitive expansion in the renderer, not the geometry

## Status

- Planning only. Related: [wgsl-strategy.md](wgsl-strategy.md) (points as quads, step 4), [storage-buffers.md](storage-buffers.md) (data arrays on WebGPU), [attribute-bindings.md](attribute-bindings.md) (vertex buffer limits, vertex pulling), [backend-split.md](backend-split.md), [minimum-capabilities.md](minimum-capabilities.md).
- Done as groundwork: spheres and points share `createCenterGroupData` (`geo/geometry/center-group.ts`: center/group buffers, bounding sphere, group mapping, `tPositionGroup` texture with optional level of detail ordering) and the shader chunk `assign_quad_position_group` (quad corner from `VertexID % 6`, item from `VertexID / 6`, position and group read from `tPositionGroup`).

## Question

Builders currently decide how a primitive is drawn: cylinders, lines and text are expanded into per-vertex data (corner mappings, duplicated item data, index buffers) when the geometry is built. Spheres and points already keep one entry per item and expand in the vertex shader. Should expansion be a renderer/backend detail for all impostor-like primitives, and what does WebGPU allow beyond storage buffers?

## Current state

| Kind | Builder output per item | Draw | Notes |
| --- | --- | --- | --- |
| spheres | center, group (one texel of `tPositionGroup`) | 6 vertices, no elements | level of detail by reordering items into strides |
| points | center, group (one texel of `tPositionGroup`) | 6 vertices, no elements | since WGSL strategy step 4 |
| cylinders | 6 vertices, each with start, end, group, scale, cap, color mode, mapping (13 floats) + 12 indices | 4 triangles, elements | about 360 bytes per cylinder instead of about 40; 9 vertex buffers (over the WebGPU default of 8), all 16 attributes when merged |
| lines | 4 vertices, each with start, end, group, mapping + 6 indices | 2 triangles, elements | |
| text | 4 vertices per character, each with position, group, mapping, depth, texture coordinates + 6 indices | 2 triangles, elements | glyph quads |

- `uVertexCount` means "data vertices" (per-vertex color/size lookups, vertex-instance stride) while the draw count means drawn vertices; points already need both (`drawCount = 6 × uVertexCount`), the emulated vertex id is sized by the larger of the two for renderables without elements, merged segments carry `drawBases` for draw ranges.
- Instances (structure transforms) are instanced attributes (`aTransform`, `aInstance`), so per-item data can't use instancing too.

## Proposal

1. **Geometry holds per-primitive data only.** Cylinders: start, end, radius scale, caps, color mode, group per cylinder; lines: start, end, group per line; text: per character position, group, glyph and offsets. No mappings, no duplicated vertices, no index buffers in the builders.
2. **The renderable declares a topology**, e.g. `quad` (6 vertices or 4 with a shared index pattern) or `box` (cylinder impostor shape) per item, plus the per-item data arrays; the backend chooses how to draw it:
   - WebGL1: emulated vertex id, item data in data textures (as spheres and points today);
   - WebGL2: `gl_VertexID`, data textures;
   - WebGPU: `vertex_index`/`instance_index`, storage buffers ([storage-buffers.md](storage-buffers.md)).
3. **Shared expansion chunks** per topology (like `assign_quad_position_group`), parametrized by the item data reader, so spheres, points, lines and text use one quad expansion and cylinders one box expansion.
4. **Counts made explicit**: separate "item count", "data vertex count" and "drawn vertex count" in the values instead of overloading `uVertexCount`/`drawCount`; merged renderables concatenate per-item data and derive draw ranges from the topology.

Expected gains: cylinders about 9× less memory and upload, under the WebGPU vertex buffer limit without vertex pulling workarounds; simpler merging (per-item data concatenates, no index rewriting); one expansion path per topology on all backends.

## Shared quad index buffer (measure first)

- Quads drawn as 6 non-indexed vertices run the vertex shader 6 times; with indices `4i + {0, 1, 2, 1, 3, 2}` only 4 unique vertices run (post-transform cache). The index pattern depends only on the item count, so one index buffer per context (grown on demand) serves all quad renderables; the vertex id then counts unique vertices (`% 4`, `/ 4`), which works with `gl_VertexID` and the emulated vertex id.
- Relevant for spheres (expensive vertex shader: projected bounding quad), less for points. Measure with the timing mode on large sphere scenes before adopting.

## WebGPU beyond storage buffers

- **Quads without index buffers**: 4-vertex triangle strips per item with the item from `instance_index`; needs instance transforms out of instanced vertex attributes (storage buffer), with item and transform derived from one index (or drawn per transform).
- **GPU-driven culling**: a compute pass for frustum, level of detail and hi-z occlusion culling that writes `drawIndirect` arguments, replacing the per-frame CPU culling and multi-draw lists; multi-draw-indirect is optional, the baseline is one indirect draw per render item.
- **Render bundles**: record static draw sequences once and replay them, lowering the CPU encoding cost for scenes with many objects.
- **Per-object data in one buffer**: all object uniforms in one uniform/storage buffer selected with dynamic offsets, avoiding per-draw uploads (global bindings already map to a bind group).
- **Base instance in core**: `firstInstance` replaces the base instance extension paths of the multi-draw code.
- **Compute shaders** for GPU surfaces and other fragment-shader compute (see the compute section of [wgsl-strategy.md](wgsl-strategy.md)).

## Steps

1. **Counts**: introduce explicit item/data-vertex/drawn-vertex counts (values and merged segments), keep `uVertexCount` as the data vertex count; unit tests for spheres, points and merged layouts.
2. **Lines** (smallest): per-line data in a texture, quad expansion chunk shared with points/spheres; builder without mappings and indices; visual tests byte-identical or within tolerance.
3. **Cylinders**: per-cylinder data (two or three texels per cylinder), box expansion chunk; builder without per-vertex duplication; attribute count drops from 9 buffers to the instance attributes; visual tests including merged and clipped cases.
4. **Text**: per-character data and glyph table; quad expansion.
5. **Shared quad index buffer**: measure (step above), adopt if it pays off.
6. **WebGPU**: storage buffers for the item data, triangle-strip quads with `instance_index`, then GPU-driven culling with indirect draws and render bundles as separate optimizations behind capability checks ([minimum-capabilities.md](minimum-capabilities.md)).

## Open questions

- Per-item data transport on WebGL: one texture per kind (position-group like) vs. the generic data arrays of [storage-buffers.md](storage-buffers.md) (`DataSpec`), which WebGL maps to data textures anyway.
- Level of detail (sphere strides) for other kinds: generalize the reordering in `createCenterGroupData` to all per-item data?
- Geo-export and other consumers of per-vertex values (cylinders `aStart`/`aEnd`, lines, text) move to the per-item buffers.
