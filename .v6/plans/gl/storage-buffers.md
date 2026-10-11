# Plan: Storage buffers for data arrays (WebGPU)

## Status

- Not started; planning only. Part of the WebGPU work in [backend-split.md](backend-split.md); rules in [minimum-capabilities.md](minimum-capabilities.md).

## Question

WebGL has no storage buffers, so per-group/instance/vertex data is pulled from "data textures" (`readFromTexture(tX, index, uXTexDim)`: index to 2D texel, `texture2D`/`texelFetch`). Should WebGPU use storage buffers instead, also for group colors and sizes?

## Answer

Yes for index-addressed data arrays, on WebGPU core. Textures stay for filtered or image-like data, for GPU-produced intermediate results until their producers are ported to compute, and for all of WebGL.

Why:

- Direct indexing (`data[i]`): no 2D layout, padding or `uXTexDim` uniforms, no float index math (precision issues above 2^24 on WebGL1-style float indices).
- Relieves the tightest limit: sampled textures per stage (16, used up to 16 by graphics shaders with all features, see minimum-capabilities). Storage buffers count separately.
- Typed access (`u32`, `f32`, packed `unpack4x8unorm`), partial updates with `writeBuffer` ranges (marker/highlight updates), no sampler state.
- Same mechanism as vertex pulling, which WebGPU needs anyway for cylinders and merged lines/text (8 vertex buffers, [attribute-bindings.md](attribute-bindings.md) step 6), and compute shaders can write the arrays directly later.

Availability (WebGPU core, default limits):

| Limit | Default | Note |
|---|---|---|
| `maxStorageBuffersPerShaderStage` | 8 | core: read-only storage buffers allowed in vertex and fragment stages |
| `maxStorageBufferBindingSize` | 128 MiB | per binding |
| `maxBufferSize` | 256 MiB | per buffer |
| `minStorageBufferOffsetAlignment` | 256 B | for windowed bindings |
| `maxDynamicStorageBuffersPerPipelineLayout` | 4 | |

- Compatibility mode (opt-in `featureLevel: 'compatibility'`, for D3D11/GLES 3.1 devices): `maxStorageBuffersInVertexStage` defaults to 0 and fragment stage to 4; it also restricts flat interpolation (`either` only) and lowers other limits. Decision: WebGPU core only (at least core features and limits); devices without core keep using WebGL2, so no texture-pulling path is needed in WGSL for CPU-provided data.
- Size: one binding holds 33.5 M packed RGBA8 entries (`u32`) or 8.4 M `vec4<f32>`; a default 8192² data texture holds 67 M texels, i.e. textures hold 2x more at default limits. Group-instance arrays grow with groups x instances, vertex pulling with vertex count (10 M vertices x 12 B positions = 120 MB). Desktop adapters usually report much higher binding sizes (requestable via `requiredLimits` up to `adapter.limits`), mobile less.
- Per the minimum-spec rule, larger arrays need a fallback that works at defaults: bind windows of the buffer per draw (instance ranges, offsets aligned to 256 B), or split render items by instance range; requesting higher adapter limits is an optimization on top.

## Design

- Binding budget (8 per stage): one data arena storage buffer per render object (like the data textures today, which are per render object; instances already share one render object) holding all its data arrays (color, size, marker, overpaint, transparency, emissive, substance, clipping, wiggle, palette, sphere `tPositionGroup`, ...), with section offsets and strides in the object uniform block; plus geometry (vertex pulling) and instance data (transforms, instance ids, segments). Target: at most 4 storage bindings per stage.
- Packing into `u32` words: RGBA8/RGB8 as one word (`unpack4x8unorm`), single bytes four per word (shift/mask; marker, transparency, emissive, clipping, wiggle), sizes as `f32`. Memory equal to or below the data textures (no 2D padding).
- Neutral layout (like `uniform-layout.ts`): `getDataLayout(schema, values)` maps `DataSpec`s to sections (format, element count, word offset); unit-tested; used by the WebGPU backend, ignored by WebGL.
- Schema: a new `DataSpec` for index-addressed data arrays that do not benefit from a 2D/3D layout, separate from `TextureSpec`, which stays for genuine 2D/3D textures (filtering, spatial sampling, render targets). WebGL maps a `DataSpec` to a data texture (today's behavior), WebGPU to an arena section.
  - `DataSpec`: `tColor`, `tSize`, `tMarker`, `tOverpaint`, `tTransparency`, `tEmissive`, `tSubstance`, `tClipping`, `tWiggle`, `tPalette`, `tPositionGroup`, image `tGroupTex`.
  - `TextureSpec`: `tColorGrid`/`t*Grid` volumes, `tGridTex`, `tTransferTex`, `tImageTex`, `tValueTex`, `tFont`, depth/dpoit/pass textures, GPU-produced textures.
- Shader access: shader code calls generated accessors (e.g. `readColor(i)`, `readMarker(i)`) instead of `readFromTexture(tColor, i, uColorTexDim)`; GLSL accessors read the data texture (with a backend-owned size uniform), WGSL accessors read the arena. Shared bodies, fewer per-backend differences for the WGSL port.
- Merged render objects: sections concatenate per segment as the data textures do today (`aSegment` offsets unchanged).
- GPU-produced data (texture-mesh `tPosition`/`tGroup`/`tNormal` from GPU isosurface, color smoothing grids, gaussian density): stay textures until the producers are compute shaders writing storage buffers.

## Steps

1. Introduce `DataSpec` with flat typed array values and convert the data arrays listed above; the WebGL backend lays them out into data textures and owns the `uXTexDim` uniforms (removed from the schemas, generated with the accessors); value producers (themes, `createValues`, merged objects) stop computing texture sizes. Unit test that every `DataSpec` is only read through its accessor.
2. Introduce accessor functions in the GLSL chunks (WebGL behavior unchanged, visual tests on all profiles).
3. Neutral data layout (`getDataLayout`) with packing rules and unit tests (offsets, word counts, merged segments).
4. WGSL: arena declaration and accessors generated from the layout; minimum-limits tests for storage buffers per stage and binding sizes.
5. Large arrays: windowed bindings per instance range at default limits; optional higher `requiredLimits`; documented in minimum-capabilities.
6. Vertex pulling on the same mechanism (attribute-bindings step 6).

## Decisions

- `DataSpec` for index-addressed data arrays, `TextureSpec` only for genuine 2D/3D textures.
- One data arena per render object for now (smallest change from the per-render-object data textures); sharing per material/visual can come later.
- WebGPU core only, no compatibility mode.
- `DataSpec` values are flat typed arrays with an element count (and element format); the 2D texture layout (width/height, padding, `uXTexDim`) is a WebGL implementation detail computed by the WebGL backend.

## Open questions

- None currently.
