# Plan: Uniform Buffer Objects (WebGL2) as a stepping stone to WebGPU

## Status

- Reconstituted from the original `plans/uniform-buffers.md` (lost when the repo was restructured into the v6
  workspace). Paths are updated to the new layout (`src/mol-gl/` → `packages/graphics/src/gl/`, `src/mol-canvas3d/` →
  `packages/graphics/src/canvas3d/`, `src/mol-util/` → `packages/core/src/util/`, `src/tests/visual/` →
  `tests/src/visual/`).
- Phases 0-3 done (uncommitted). Follow-up plans: [texture-bindings.md](texture-bindings.md),
  [attribute-bindings.md](attribute-bindings.md) (Phase 4, incl. an evaluation for WebGL + WebGPU).
- Items lost in the restructure are restored: `webgpu/uniform-declarations.ts`, `getBindingLayout` in
  `uniform-layout.ts`, unit tests `packages/graphics/src/gl/_test/uniform-buffer.test.ts` + `shader-uniforms.test.ts`,
  visual cases `geometry/spheres/lod-levels`, `renderer/clipping/pixel|instance`, and the `debug=1` URL parameter.
- Rule: backend-specific code (GLSL, WGSL, GL/GPU calls) lives only in `webgl/` and `webgpu/`; `uniform-layout.ts` stays
  backend-neutral.
- Rule: backend targets and minimum specs, see [minimum-capabilities.md](minimum-capabilities.md).

## TL;DR

CPU-computed std140 layouts (identical to the WGSL uniform address space layout except `bool`), fixed binding indices
per update-frequency tier (`Global` = 0, `Material` = 1, `Object` = 2, ≈ WebGPU bind groups), version-tracked dirty
writes, sub-allocated buffers bound as ranges (≈ dynamic offsets), a ring buffer for per-pass global snapshots, and
uniform declarations generated from the schema for both WebGL backends (WGSL generator too). The schema is the single
source of truth for non-sampler uniforms.

## 1. What the change does

| Area                                                                                         | Change                                                                                                                                                                                                                                                                                                                                                      |
| -------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `packages/graphics/src/gl/renderable/schema.ts`                                              | `UniformSpec(kind, variant?, block?, length?)`, block defaults to `'material'` (material variant) or `'object'`; `'draw'` = classic per-draw uniforms (`uDrawId`, `uLod`, `uLodLevel`); array uniforms need a max `length` to be in a block. Helpers `GlobalUniformSpec`/`MaterialUniformSpec`/`ObjectUniformSpec`/`DrawUniformSpec`, `MaxLightCount = 16`. |
| `packages/graphics/src/gl/uniform-layout.ts` (new, neutral)                                  | Block info/tiers, std140 layouts (incl. arrays), `getUniformBlock`, `createUniformBlockData`, `createUniformBlockWriter` (packs into any `ArrayBuffer` region).                                                                                                                                                                                             |
| `packages/graphics/src/gl/webgl/uniform-buffer.ts` (new)                                     | `UniformArena` (64 KiB slabs + CPU mirrors, power-of-two size classes aligned to `UNIFORM_BUFFER_OFFSET_ALIGNMENT`, dirty-range uploads, `bindBufferRange`), `createUniformBuffer` (slice), `createUniformRingBuffer` (snapshot per change).                                                                                                                |
| `packages/graphics/src/gl/webgl/shader-declarations.ts` (new, was `uniform-declarations.ts`) | GLSL generation: `getUniformBlockDeclaration`, `getShaderUniformDeclarations(schema, useBlocks)`; also samplers and attributes (`getShaderDeclarations`).                                                                                                                                                                                                   |
| `packages/graphics/src/gl/shader-code.ts`                                                    | `ShaderUniformDeclarations`, `stripUniformDeclarations` (removes hand-written declarations of generated names, throws on type mismatch), injection in `addShaderDefines` for both backends.                                                                                                                                                                 |
| `packages/graphics/src/gl/webgl/program.ts`                                                  | Locations skip block members, `uniformBlockBinding`, debug checks (std140 offsets/strides vs driver, missing-location throw, classic-uniform warning).                                                                                                                                                                                                      |
| `packages/graphics/src/gl/webgl/render-item.ts`                                              | Routing by block, object/material buffers created from the layout, `buffered` uniforms in the object block.                                                                                                                                                                                                                                                 |
| `packages/graphics/src/gl/renderer.ts`                                                       | Global ring buffer, light count clamp, `dispose` destroys the buffer.                                                                                                                                                                                                                                                                                       |
| `packages/graphics/src/gl/webgl/resources.ts`, `state.ts`, `context.ts`                      | Lazy arena per context, `uniformBuffer`/`uniformRingBuffer`/`materialUniformBuffer` (cached per `(materialId, layout)`), program cache key hashes the generated declarations, `bindUniformBufferRange` with cache, `uniform` resource stats.                                                                                                                |
| Shaders (`packages/graphics/src/gl/shader/**`)                                               | Graphics chunks and graphics shaders no longer hand-declare non-sampler, non-array uniforms.                                                                                                                                                                                                                                                                |
| Passes (`packages/graphics/src/canvas3d/passes/`)                                            | Array lengths for bloom/shadow/ssao/tracing (+ clamps).                                                                                                                                                                                                                                                                                                     |
| `packages/core/src/util/clip.ts`                                                             | `Clip.MaxObjectCount = 16` (clamped in `getClip`).                                                                                                                                                                                                                                                                                                          |

## 2. Keep (good decisions for WebGPU)

- CPU-side std140 layout instead of querying the driver: deterministic, backend-neutral, reusable for WGSL.
- Fixed binding indices per tier → `@group(0..2)`.
- Update-frequency tiers match WebGPU's bind-group-per-frequency idiom.
- Version-tracked writers resolve cells once; the per-draw cost is a version compare loop.
- WebGL1 keeps classic uniforms, but declarations come from the same generator.
- Context loss: CPU mirrors are re-uploaded, binding cache reset, programs re-bind blocks.

## 3. Review findings (all fixed in Phase 0 unless noted)

- **F1** per-level `uLod` override silently dropped on WebGL2 (`mdbData.uniforms` → no location → no-op). Fixed by the
  `'draw'` tier + debug throw on `undefined` locations.
- **F2** program cache key ignored block layouts → now hashes the generated declarations.
- **F3** buffers created only when values were present → now created from the layout.
- **F4** `block`/`variant` independent → `UniformSpec` validates combinations, routing by block.
- **S1** material buffers are shared per `(materialId, layout)`; items sharing a material must have identical material
  values (stricter than before). Not changed; candidate for a debug consistency check.
- **S2** in-place global rewrites → fixed in Phase 2 (ring buffer).
- **S3** `Renderer.dispose` destroys the global buffer. **S4** `ivec2/3/4` GLSL types. **S5** type check when stripping
  declarations (always, not only debug).
- **S6** `bool` stored as int; WGSL uses `u32`.
- **S7** one GL buffer per item → fixed in Phase 2 (arena).
- **S8** runtime toggle to disable UBOs: not done (would allow a `webgl2-no-ubo` visual profile).
- **S9** stats: not done (`uniformBufferCount` in `RendererStats`, upload/bind counters in timing mode).

## 4. Phases

### Phase 0 - Make the change mergeable (done)

`uLod` per-draw + missing-location debug throw; layout hash in program key; buffers from layout; spec validation;
layouts cached per schema (`WeakMap`); `Renderer.dispose`; `ivec*` types; unit tests; visual case
`geometry/spheres/lod-levels` (fails with a 28% diff without the F1 fix) and `debug=1` harness parameter.

### Phase 1 - Complete UBO coverage on WebGL2 (done)

- Blocks by default: all renderable-specific and compute-pass uniforms are in blocks without per-file edits, keeping
  update semantics (nothing promoted to `material`).
- Arrays in blocks with max lengths: lights 16, clip objects 16, SSAO samples 256 / levels 16 / blur kernel 25, bloom
  coefficients 11 / factors 5. Objects beyond a max are ignored on all backends (behavior change).
- `buffered` uniforms (`uGeoTexDim`) in the object block via their front copies.
- Debug warning for classic non-`draw` uniforms on WebGL2; driver check covers array strides.
- Visual cases `renderer/clipping/pixel|instance`.
- Not covered by visual tests: picking, overpaint/substance, hi-z, gaussian density GPU, color smoothing, skybox/image
  backgrounds.

### Phase 2 - Backend-neutral uniform layer (done)

- Neutral `uniform-layout.ts` vs. `webgl/uniform-buffer.ts`.
- Arena sub-allocation (object/material slices), ring buffer (64 slots) for the renderer's global block → no in-place
  updates of data used by earlier draws.
- Decision: no global tier split (one ~1.8 KiB snapshot per change is cheap).
- `getBindingLayout(schema)` (uniform block + textures per tier: `GlobalTextureSchema` → global, `variant: 'material'` →
  material, rest → object).
- Open:
  - Slabs are never released (high-water mark).
  - Coalesced uploads only help when several writes precede a `sync`; a write pre-pass would give one upload per slab
    per frame.
  - Object/material slices are still updated in place when values change between draws within a frame (bloom per mip,
    per-iteration illumination uniforms) - WebGPU needs versioned slices.
  - Ring wrap-around relies on WebGL command ordering; WebGPU needs a per-frame reset after submission.
  - Done (backend split step 2): the `UniformBuffer` interface and `UniformKind`/`UniformsList` (and the other neutral
    uniform value types) live in `gl/uniform.ts`.

### Phase 3 - Schema as single source of truth for shader interfaces (done)

- `getShaderUniformDeclarations(schema, useBlocks)` (cached per schema): WebGL2 → blocks + classic `draw` uniforms;
  WebGL1 → classic declarations for all non-array uniforms; preceded by `precision highp float/int` (needed before
  declarations in GLSL ES 1.0 fragment shaders).
- Migrated: `common_vert_params`, `common_frag_params`, `color_vert_params`, `color_frag_params`, `size_vert_params`,
  `light_frag_params`, `common_animation` and the mesh/spheres/cylinders/points/lines/text/image/direct-volume shaders.
  Remaining hand-written: define-sized arrays (lights, clip objects; WebGL1 uniform-vector limits), stripped on WebGL2.
  Samplers and attributes are generated too ([texture-bindings.md](texture-bindings.md) T2,
  [attribute-bindings.md](attribute-bindings.md) step 4).
- Guard test: every `uXxx` identifier in any preprocessor branch of the graphics shaders is in the schema
  (`dGeometryType_*` branches resolved per geometry, allowlist `uBumpFrequency`, `uBumpAmplitude`; extended to
  `tXxx`/`aXxx`).
- WGSL: `getWgslUniformDeclarations(schema)` in `webgpu/uniform-declarations.ts` → `struct Global/Material/Object` +
  `@group(0..2) @binding(0) var<uniform> uGlobal/uMaterial/uObject`, `bool` → `u32`, arrays with ≤ 2 components per
  element padded to `vec4`; tests check WGSL layout rules give the std140 offsets for all graphics schemas; `draw`
  uniforms and arrays without length reported as `unsupported`.
- Access naming decision: GLSL keeps anonymous blocks (bare names); WGSL uses instance access (`uGlobal.uView`); no
  `#define` aliasing.
- Open: migrate compute/postprocessing shaders (still hand-declared, stripped + type-checked); textures/samplers
  ([texture-bindings.md](texture-bindings.md)); a WebGPU `draw` tier (per-draw slice with dynamic offset or
  instance-indexed storage).

### Phase 4 - Full shader interface from the schema: attributes (planned)

Separate plan: [attribute-bindings.md](attribute-bindings.md) (includes the evaluation of the schema as source of truth
for a WebGL + WebGPU backend).

### Phase 5 - WebGPU prerequisites (out of scope, informed by this)

- Pipeline/bind group layout objects derived from schema tiers.
- Shader strategy: hand-ported WGSL vs. GLSL→WGSL translation (naga/tint/glslang WASM); generated declarations are
  needed either way.
- Compute passes as real compute shaders where beneficial.

## 5. WebGL2 ↔ WebGPU mapping

| Concept       | WebGL2 (now)                              | WebGPU                                    |
| ------------- | ----------------------------------------- | ----------------------------------------- |
| Tier binding  | `uniformBlockBinding` fixed index         | `@group(n)`                               |
| Global data   | ring buffer snapshots + `bindBufferRange` | ring + dynamic offset (per-frame reset)   |
| Material data | arena slice per `(materialId, layout)`    | slice, dynamic offset or own bind group   |
| Object data   | arena slice per render item               | dynamic offset                            |
| Per-draw data | classic uniforms (`draw` tier)            | dynamic offset / instance-indexed storage |
| Layout        | std140 on CPU                             | same rules (`bool` → `u32`)               |
| Declarations  | generated from schema (GLSL)              | generated WGSL                            |

## 6. Validation

- Unit tests (to be re-added): std140 layouts incl. arrays and the `GlobalUniformSchema` snapshot (size 1600), packers,
  spec defaults/validation, declaration generation for both backends, arena (alignment, reuse/zeroing, coalescing,
  reset), ring (slots, wrap-around), WGSL layout equality, shader identifier guard.
- Debug mode (`debug=1`): layout/stride check, active uniforms, missing-location throw, classic-uniform warning,
  stripped-type check.
- Visual tests (`tests/src/visual`): `?target=webgl1&reference=webgl2&debug=1` (classic vs. UBO, 61 cases passed before
  the restructure), `webgl1-minimal` and `webgl2-no-float` (only the known SSAO multi-scale and slice/image errors).
- Not done: context loss test (`WEBGL_lose_context`), perf comparison on scenes with many render items.

## 7. Open questions

- Global buffer per context vs. per renderer (per renderer today).
- 256-byte alignment overhead per object slice vs. packing several objects per slice.
- How long WebGL1 stays supported (determines whether classic declarations and define-sized array declarations must
  remain).
