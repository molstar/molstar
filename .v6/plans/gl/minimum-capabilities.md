# Plan: Backend targets and minimum capabilities

Applies to all plans in this folder ([uniform-buffers.md](uniform-buffers.md),
[texture-bindings.md](texture-bindings.md), [attribute-bindings.md](attribute-bindings.md)).

## Targets

- WebGL1 and WebGL2 stay supported for the foreseeable future; a WebGPU backend is added next to them, not instead of
  them.
- WebGPU means core (features and limits); compatibility mode is not targeted, devices without core WebGPU use WebGL2
  ([storage-buffers.md](storage-buffers.md)).
- Optional and upcoming WebGPU features (e.g. bindless textures, (multi-)draw indirect, 64-bit atomics) may be used, but
  only as optimizations: behind a capability check, with the baseline path kept as fallback. The same applies to
  optional WebGL extensions.
- Consequences for the shared, backend-neutral layer (schema, `uniform-layout.ts`, `attribute-layout.ts`,
  `getBindingLayout`):
  - It describes the baseline (per-program bindings, CPU-driven draws, one texture binding per texture); optional paths
    are derived from it in the backend folders, never the other way round.
  - Capabilities are explicit inputs (like `defines`), so a generator can pick the baseline or an optimized path per
    device.

## Rule: minimum specs

Everything must work on WebGL2 minimums and WebGPU default limits without optional features. Anything needing more
capabilities or bigger limits must have a fallback, and both are documented (table below, or the plan of the feature).
WebGL1 is supported with the extensions Mol* requires (visual profile `webgl1-minimal`): `ANGLE_instanced_arrays`,
`OES_element_index_uint`, `OES_standard_derivatives`, `OES_texture_float`, `EXT_blend_minmax` (all core in WebGL2 and
WebGPU, all supported by headless-gl; context creation fails without them, see `createExtensions`). The minimum-limit
guarantee is defined for WebGL2 and WebGPU. A known WebGL1 limit: the emulated vertex id (`aVertex`, a float attribute)
is exact up to 2^24 vertices, e.g. about 2.8M points; larger draws warn and get wrong ids (not worth a wider attribute
for WebGL1).

## Minimum capabilities

| Resource                                                 | WebGL2 minimum                                                                | WebGPU default                                                       | Current use                                                                                                                                                                          | Status / fallback                                                                                                                                                    |
| -------------------------------------------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Vertex attributes                                        | 16 (`MAX_VERTEX_ATTRIBS`)                                                     | 16 (`maxVertexAttributes`)                                           | ≤ 16, merged cylinders exactly 16                                                                                                                                                    | OK (tested), no headroom left for cylinders                                                                                                                          |
| Vertex buffers                                           | -                                                                             | 8 (`maxVertexBuffers`)                                               | one buffer per attribute; exceeded by cylinders (9, merged 12), merged lines and text                                                                                                | Open: interleaving or storage buffers + vertex pulling ([attribute-bindings.md](attribute-bindings.md) step 6); exceedances listed in the test                       |
| Uniform block size                                       | 16 KiB (`MAX_UNIFORM_BLOCK_SIZE`)                                             | 64 KiB (`maxUniformBufferBindingSize`)                               | global 1600 B, material with 16 clip objects ≈ 2.5 KiB                                                                                                                               | OK (tested)                                                                                                                                                          |
| Uniform blocks per stage / bindings                      | 12 / 24                                                                       | 12 per stage                                                         | 3 (global, material, object)                                                                                                                                                         | OK                                                                                                                                                                   |
| Uniform buffer offset alignment                          | implementation-defined (≤ 256 in practice)                                    | 256 (`minUniformBufferOffsetAlignment`)                              | arena aligns to the queried value                                                                                                                                                    | OK                                                                                                                                                                   |
| Bind groups                                              | -                                                                             | 4 (`maxBindGroups`)                                                  | 3 tiers                                                                                                                                                                              | OK, one group left                                                                                                                                                   |
| Sampled textures per stage                               | 16 (`MAX_TEXTURE_IMAGE_UNITS`, `MAX_VERTEX_TEXTURE_IMAGE_UNITS`), 32 combined | 16 (`maxSampledTexturesPerShaderStage`, `maxSamplersPerShaderStage`) | worst-case feature defines: ≤ 16 used per stage for every graphics shader and variant; `getBindingLayout(schema, defines)` with the current `usedIf` annotations binds ≤ 16 in total | OK (tested); keep `usedIf` up to date when adding textures                                                                                                           |
| Varyings / inter-stage variables                         | 15 vectors (`MAX_VARYING_VECTORS`)                                            | 16 (`maxInterStageShaderVariables`)                                  | WebGL2 (packed): ≤ 15 everywhere; unpacked: cylinders color/tracing variants need 18 locations                                                                                       | WebGL2 OK (tested); WebGPU open: pack scalar varyings into `vec4`s in WGSL (exceedances listed in the test)                                                          |
| Draw buffers / color attachments                         | 4                                                                             | 8                                                                    | wboit/dpoit/illumination gated by `isSupported`                                                                                                                                      | Fallbacks exist (blended transparency, no illumination); document per feature                                                                                        |
| 2D / 3D texture size                                     | 2048 / 256                                                                    | 8192 / 2048                                                          | data textures sized from `maxTextureSize`, 3D grids from `max3dTextureSize`                                                                                                          | Check that every 3D path has a 2D fallback above 256                                                                                                                 |
| Float textures as render targets, float linear filtering | optional (`EXT_color_buffer_float`, `OES_texture_float_linear`)               | `float32-filterable`, `float32-blendable` optional                   | GPU isosurface, gaussian density, SSAO, dpoit                                                                                                                                        | Fallbacks: CPU surfaces, blended transparency (tested by `fallback/**` cases); SSAO packs its downsampled/half/quarter depth into RGBA8 without float render targets |
| Float textures (sampling), min/max blending              | core (WebGL1: required `OES_texture_float`, `EXT_blend_minmax`)               | core                                                                 | image `tValueTex` (iso level), grid3d inputs, GPU surfaces, emissive transparency                                                                                                    | OK, required (no fallbacks)                                                                                                                                          |

## Validation

- Unit tests `packages/graphics/src/gl/_test/minimum-limits.test.ts` check every graphics shader against the table:
  attribute locations and WebGPU vertex buffers per schema (plain and merged), uniform block sizes, and per program (all
  render variants, worst-case feature defines, WebGL2 prefix defines) the used samplers per stage, packed and unpacked
  varyings, the `getBindingLayout` texture count and that no used texture is excluded by `usedIf`. Known exceedances are
  listed explicitly with their fallback, so new ones fail and fixed ones must be removed. Shaders are resolved with a
  small GLSL preprocessor (`_test/graphics-shaders.ts`).
- Visual test profile `webgl2-minimum`: WebGL2 without any extension and `getParameter` limits clamped to the WebGL2
  minimums (`UNIFORM_BUFFER_OFFSET_ALIGNMENT` raised to 256). Result (`debug=1`): no errors or failures on any profile
  pair; vs. `webgl2` 59 pass, 12 skipped for missing capabilities.
- Fallback cases (`tests/src/visual/cases/fallback.ts`, `VisualTestCase.fallback`): a profile without
  `fallback.requires` renders the case as is (the feature falls back internally) and is compared against a profile with
  the capability rendering `fallback.overrides` (the documented fallback). Covered: wboit and dpoit → blended,
  illumination → off, GPU isosurface and gaussian surface → CPU; all pass on `webgl2-minimum`, `webgl2-no-float`,
  `webgl1-minimal`.
- The per-program unit tests also check that every `v*` identifier used in a stage is declared, for WebGL2 and for
  WebGL1 without extensions (found: image picking without draw buffers read an undeclared `vColor`; direct-volume with
  `dClipping` read an undeclared `vClipping`, now ignored for direct volumes).
- Clamped limits are reported, not enforced: a shader exceeding them still compiles on a bigger GPU, which is why
  per-stage samplers and varyings are covered by the unit tests.
- With `debug=1`, `checkError` consumes GL errors after each draw and only logs them; the visual test harness captures
  these `WebGL error:` logs as errors (before, such errors were invisible in debug runs).
- Optional-feature paths (once they exist) need a visual test profile with the feature and one without, compared against
  each other.
