# Plan: Storage buffers for data arrays (WebGPU)

## Status

- **Current**: done, in part by other work. On WebGPU the index-addressed theme layers (tColor, tSize, tMarker, the
  layers, tPositionGroup) are read from storage since [shared-draw-resources.md](shared-draw-resources.md) step 3.1:
  shared per-context data pages (`webgpu/data-layers.ts`, an item's layers in one range, word offsets in object
  uniforms, range uploads), generated `read_tX(i)` accessors (sizes as rgba8), not an arena per render object. GPU
  surfaces read a vertex arena ([gpu-surface-buffers.md](gpu-surface-buffers.md)), GPU culling writes instance data into
  storage buffers used as vertex buffers ([merged-render-objects.md](merged-render-objects.md) step 6). Palettes and
  images stay sampled textures.
- Steps 1 and 2 done: the layers are `DataSpec`s (`rgba8`, `r8`, `rgba32f`), their values `DataArray`s (`{ array,
  count }`, no 2D layout, no padding; producers allocate whole words), read by accessors `read_<name>(i)` in GLSL and
  WGSL alike. WebGL lays each out as a data texture (`webgl/data-textures.ts`: the schema's data specs become textures
  and object size uniforms `u<Name>TexDim` of the backend, the GLSL accessors generated with the declarations; rows
  uploaded from views of the array, the last row's elements apart, no padded copy, `Texture.loadData`); WebGPU reads
  them from the pages, direct volumes too now (their trilinear vertex colors read element by element). The grids'
  2D layouts have uniforms of their own (`u…GridTexDim`, were the layers'). Deviations from the design below:
  `tPalette` (sampled, filtered) and `tGroupTex` (image) stay textures; no neutral `getDataLayout` (step 3), each
  backend lays out by count. Compute consumers that load a layer into a texture of their own (texture-mesh color
  smoothing) lay it out with `getDataTextureImage` (a copy). Frames, recolor and hover unchanged within the noise on
  WebGL 2 and WebGPU; WebGL's texture count per stage is unchanged (the minimum-limits test counts the accessors).
- WebGL data pages, as WebGPU's: the data arrays in ranges of textures shared per context and format
  (`webgl/data-pages.ts`: 1M elements of the 8-bit formats, 256K of `rgba32f`, 2048 wide; larger arrays a texture of
  their own; the range allocator shared with WebGPU, `range-allocator.ts`), offsets in object uniforms `u<Name>Offset`
  (integer math on WebGL2, float math on WebGL1, whose integers may be 16-bit), writes by element range from any offset.
  Programs skip binding a page where its unit has it (`WebGLState.boundTextures`, pages only: they are never render
  targets). Measured, three interleaved runs: WebGL2 1000 unmerged spacefill units 10.7 to 2.9 ms per frame, hover 21.6
  to 5.4, cartoon 8.4 to 4.2, `boxes-1000` 7.5 to 4.8; WebGL1 spacefill 10.9 to 6.2; merged objects unchanged. Found on
  the way: allocating a texture's storage left it bound at the active unit, where a sampler no pass provided (the DPOIT
  textures before the first DPOIT pass) read it while drawn into, a feedback loop; now unbound, and a pass unbinds the
  global samplers it does not provide.
- Step 5 done in this form: an item's layers together beyond the binding limit (`WebGPUContext.dataBindingSize`, the
  device's `maxStorageBufferBindingSize` or less by option) are packed largest first into up to 4 pages, each within it
  (`sDataPage`, `sDataPage1..3`); the layout is a define (`dDataLayerPages`, 2 bits per layer) present only where not
  0, so items within the limit keep their programs and precompile keys. Items beyond it compile their programs on
  demand, not precompiled (their layout is known only from the data), again when an update changes it; their draws are
  skipped until the pipelines are ready. A single layer beyond the limit is not drawn, with a warning: windows of one
  layer are not done (a group-instance layer of 128 MiB is 32M elements). Updates skip an item whose layers' versions
  are unchanged (the update stage of unmerged recolor of 1000 units 6.0-6.5 to 4.7-5.3 ms; rendering unchanged).
  Multi-page draws measured with the limit lowered (`pnpm benchmark --data-binding 1572864`): the 300k-atom globule's
  sphere items (1.9 MB of layers) read color and size from a second page; `globule-300k-spacefill` 4.39 / 4.39 ms,
  `-large` 14.61 / 14.62, the images the same. Where the design below says an arena per render object, the code shares
  pages per context.

## Question

WebGL has no storage buffers, so per-group/instance/vertex data is pulled from "data textures"
(`readFromTexture(tX, index, uXTexDim)`: index to 2D texel, `texture2D`/`texelFetch`). Should WebGPU use storage buffers
instead, also for group colors and sizes?

## Answer

Yes for index-addressed data arrays, on WebGPU core. Textures stay for filtered or image-like data, for GPU-produced
intermediate results until their producers are ported to compute, and for all of WebGL.

Why:

- Direct indexing (`data[i]`): no 2D layout, padding or `uXTexDim` uniforms, no float index math (precision issues above
  2^24 on WebGL1-style float indices).
- Relieves the tightest limit: sampled textures per stage (16, used up to 16 by graphics shaders with all features, see
  minimum-capabilities). Storage buffers count separately.
- Typed access (`u32`, `f32`, packed `unpack4x8unorm`), partial updates with `writeBuffer` ranges (marker/highlight
  updates), no sampler state.
- Same mechanism as vertex pulling, which WebGPU needs anyway for cylinders and merged lines/text (8 vertex buffers,
  [attribute-bindings.md](attribute-bindings.md) step 6), and compute shaders can write the arrays directly later.

Availability (WebGPU core, default limits):

| Limit                                       | Default | Note                                                                  |
| ------------------------------------------- | ------- | --------------------------------------------------------------------- |
| `maxStorageBuffersPerShaderStage`           | 8       | core: read-only storage buffers allowed in vertex and fragment stages |
| `maxStorageBufferBindingSize`               | 128 MiB | per binding                                                           |
| `maxBufferSize`                             | 256 MiB | per buffer                                                            |
| `minStorageBufferOffsetAlignment`           | 256 B   | for windowed bindings                                                 |
| `maxDynamicStorageBuffersPerPipelineLayout` | 4       |                                                                       |

- Compatibility mode (opt-in `featureLevel: 'compatibility'`, for D3D11/GLES 3.1 devices):
  `maxStorageBuffersInVertexStage` defaults to 0 and fragment stage to 4; it also restricts flat interpolation (`either`
  only) and lowers other limits. Decision: WebGPU core only (at least core features and limits); devices without core
  keep using WebGL2, so no texture-pulling path is needed in WGSL for CPU-provided data.
- Size: one binding holds 33.5 M packed RGBA8 entries (`u32`) or 8.4 M `vec4<f32>`; a default 8192² data texture holds
  67 M texels, i.e. textures hold 2x more at default limits. Group-instance arrays grow with groups x instances, vertex
  pulling with vertex count (10 M vertices x 12 B positions = 120 MB). Desktop adapters usually report much higher
  binding sizes (requestable via `requiredLimits` up to `adapter.limits`), mobile less.
- Per the minimum-spec rule, larger arrays need a fallback that works at defaults: bind windows of the buffer per draw
  (instance ranges, offsets aligned to 256 B), or split render items by instance range; requesting higher adapter limits
  is an optimization on top.

## Design

- Binding budget (8 per stage): one data arena storage buffer per render object (like the data textures today, which are
  per render object; instances already share one render object) holding all its data arrays (color, size, marker,
  overpaint, transparency, emissive, substance, clipping, wiggle, palette, sphere `tPositionGroup`, ...), with section
  offsets and strides in the object uniform block; plus geometry (vertex pulling) and instance data (transforms,
  instance ids, segments). Target: at most 4 storage bindings per stage.
- Packing into `u32` words: RGBA8/RGB8 as one word (`unpack4x8unorm`), single bytes four per word (shift/mask; marker,
  transparency, emissive, clipping, wiggle), sizes as `f32`. Memory equal to or below the data textures (no 2D padding).
- Neutral layout (like `uniform-layout.ts`): `getDataLayout(schema, values)` maps `DataSpec`s to sections (format,
  element count, word offset); unit-tested; used by the WebGPU backend, ignored by WebGL.
- Schema: a new `DataSpec` for index-addressed data arrays that do not benefit from a 2D/3D layout, separate from
  `TextureSpec`, which stays for genuine 2D/3D textures (filtering, spatial sampling, render targets). WebGL maps a
  `DataSpec` to a data texture (today's behavior), WebGPU to an arena section.
  - `DataSpec`: `tColor`, `tSize`, `tMarker`, `tOverpaint`, `tTransparency`, `tEmissive`, `tSubstance`, `tClipping`,
    `tWiggle`, `tPalette`, `tPositionGroup`, image `tGroupTex`.
  - `TextureSpec`: `tColorGrid`/`t*Grid` volumes, `tGridTex`, `tTransferTex`, `tImageTex`, `tValueTex`, `tFont`,
    depth/dpoit/pass textures, GPU-produced textures.
- Shader access: shader code calls generated accessors (e.g. `readColor(i)`, `readMarker(i)`) instead of
  `readFromTexture(tColor, i, uColorTexDim)`; GLSL accessors read the data texture (with a backend-owned size uniform),
  WGSL accessors read the arena. Shared bodies, fewer per-backend differences for the WGSL port.
- Merged render objects: sections concatenate per segment as the data textures do today (`aSegment` offsets unchanged).
- GPU-produced data (texture-mesh `tPosition`/`tGroup`/`tNormal` from GPU isosurface, color smoothing grids, gaussian
  density): stay textures until the producers are compute shaders writing storage buffers.

## Steps

1. (done, but `tPalette` and `tGroupTex`, see Status) Introduce `DataSpec` with flat typed array values and convert the
   data arrays listed above; the WebGL backend lays them out into data textures and owns the `uXTexDim` uniforms
   (removed from the schemas, generated with the accessors); value producers (themes, `createValues`, merged objects)
   stop computing texture sizes. Unit test that every `DataSpec` is only read through its accessor
   (`data-accessors.test.ts`).
2. (done) Introduce accessor functions in the GLSL chunks (WebGL behavior unchanged, visual tests on all profiles).
3. (not needed) Neutral data layout (`getDataLayout`): each backend lays out by count (WebGL rows, WebGPU pages).
4. (done, as pages) WGSL accessors generated from the schema; minimum-limits tests for storage buffers per stage and
   binding sizes.
5. (done, by pages) Large arrays: an item's layers across up to 4 bound pages, each within the binding limit (see
   Status); one layer beyond the limit is not drawn (windows of one layer not done); higher limits with
   `limits: 'adapter'`.
6. Vertex pulling on the same mechanism (attribute-bindings step 6).

## Decisions

- `DataSpec` for index-addressed data arrays, `TextureSpec` only for genuine 2D/3D textures.
- One data arena per render object for now (smallest change from the per-render-object data textures); sharing per
  material/visual can come later: planned in [shared-draw-resources.md](shared-draw-resources.md) (the data arrays as
  ranges of a store shared with the geometry, by representation and material, updated by range).
- WebGPU core only, no compatibility mode.
- `DataSpec` values are flat typed arrays with an element count (and element format); the 2D texture layout
  (width/height, padding, `uXTexDim`) is a WebGL implementation detail computed by the WebGL backend.

## Open questions

- None currently.
