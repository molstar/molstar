# Plan: Generated declarations for compute and pass shaders

## Status

- Done: all steps (the image `tValueTex` follow-up is resolved by requiring `OES_texture_float`, see [minimum-capabilities.md](minimum-capabilities.md)). Follow-up of [texture-bindings.md](texture-bindings.md) T2 and [attribute-bindings.md](attribute-bindings.md) step 4: graphics shaders get all uniform, sampler and attribute declarations from the schema; compute and pass shaders still hand-declare them, and `getShaderDeclarations` strips those at runtime.
- Rules: [minimum-capabilities.md](minimum-capabilities.md); backend-specific code only in `webgl/` and `webgpu/`.
- Why not done together with the graphics shaders: no test can reach these shader/schema pairs, so a wrong spec (dimension, integer-ness, attribute type) or a missing schema entry would only show up when that program is created at runtime, and only some passes are covered by visual tests.

## Current state

- 35 `ShaderCode` definitions outside the graphics shaders (plus grid3d), almost all following one pattern (module-level schema const, `ShaderCode` const, `getXRenderable(ctx, ...)` that spreads the schema and calls `createComputeRenderItem`):
  - `canvas3d/passes`: background, bloom (luminosity, blur, composite), cas, dof, dpoit (blend-back, evaluate), fxaa, illumination (compose, multi-sample compose), marking (edge, dim, overlay), multi-sample (compose), outline, postprocessing, shadow, smaa (edges, weights, blend), ssao (ssao, blur), tracing (trace, accumulate), wboit (evaluate)
  - `gl/compute`: copy, hi-z, histogram-pyramid (reduction, sum), marching-cubes (active-voxels, isosurface), grid3d
  - `geo`: gaussian-density, texture-mesh color-smoothing (accumulate, normalize)
- Deviations from the pattern:
  - `ShaderCode` created inside the factory function (outline, postprocessing, shadow, ssao, color-smoothing, histogram-pyramid, marching-cubes, gaussian-density): new `ShaderCode` per call, so a new id and program cache entry per renderable; can be hoisted.
  - `gl/compute/grid3d.ts`: schema and shader code are generated per call from a template; needs its own handling.
- Visual test coverage (`tests/src/visual`): most passes are covered (postprocessing, transparency, illumination, marking, multi-sample, background, GPU isosurface, gaussian density); not covered: hi-z, color smoothing, grid3d.

## Formalize the pattern

The pattern is already uniform; make the pair addressable from tests without changing behavior:

1. Hoist every schema and `ShaderCode` to module level and export them. Done as plain exports of both consts (`export const BlendSchema`, `export const BlendShaderCode`), paired in the test list; no wrapper object or helper yet. Optionally later a tiny helper (`ComputeShader(name, vert, frag, schema, extensions?, outTypes?)` returning `{ shaderCode, schema }`) and a `createComputeRenderItem` overload taking it; only if it removes boilerplate.
2. Per-shader define variants for tests: string defines enumerate their `DefineSpec` options, booleans both values, number defines (e.g. `dNSamples`, `dLevels`, `dEdgeThreshold`) need example values; keep them next to the export or in the test list.
3. grid3d: the generated schema/code can't be listed statically; test the generator with representative inputs instead.

## Steps

1. **Test list** (done): `_test/compute-shaders.ts` lists every pair (`ComputeShaders`) and the defining modules (`ShaderModules`); `compute-shaders.test.ts` scans `packages/graphics/src` for `ShaderCode(` calls and requires each to be an `export const X = ShaderCode(` that is listed in `GraphicsShaders` or `ComputeShaders` (and vice versa); grid3d is excluded. Define variants come with step 2.
2. **Run the existing checks over the list** (done, no shader changes): define variants from `getComputeDefineVariants` (base: booleans on, first string option, example numbers; plus one variant per other string option / disabled boolean) and prefix defines from `getComputePrefixDefines` (minimum-spec WebGL2 and WebGL1).
   - `shader-uniforms.test.ts`: identifier guard over all branches; `getShaderDeclarations` + the strip functions on the raw sources (as in `addShaderDefines`) for every variant on WebGL2 and WebGL1.
   - `minimum-limits.test.ts`: attribute locations and vertex buffers, uniform block sizes; per variant and version: undeclared varyings, textures per stage, packed varyings / inter-stage locations, WebGPU binding layout, samplers missing from the schema or excluded by `usedIf`.
   - Findings: `dThicknessMode` (tracing) had no options (added); local variables named like uniforms/attributes, listed in `ComputeIdentifierFindings`: `uOffset` (fxaa `sampleLuma` parameter), `aBuff` (illumination compose). No sampler/attribute/uniform type mismatches, undeclared varyings or limit exceedances.
3. **Fix findings** (done): renamed `uOffset`/`vOffset` → `offsetX`/`offsetY` (fxaa) and `aBuff` → `colorBuff` (illumination compose), `ComputeIdentifierFindings` removed.
4. **Hoist** the per-call `ShaderCode` definitions (done; fewer program cache entries, same output).
5. **Remove hand-written declarations** (done): samplers, attributes and non-array uniforms removed from all compute/pass shaders (incl. the shared `quad.vert`); define- or constant-sized arrays stay (`uKernel`, `uSamples`, `uLevelRadius/Bias`, `uLightDirection/Color`, `uGaussianCoefficients`, `uBloomFactors/Tints`). `precision ... sampler*` statements removed unless the sampler type is still used for function parameters (smaa, ssao, illumination compose, marching cubes). The "no hand-written declarations" guard now covers graphics and compute shaders and non-array uniforms. Visual tests pass on webgl1, webgl2-minimum, webgl1-minimal and webgl2-no-float (debug mode).
6. **grid3d** (done): `getGrid3dComputeSchema(spec, isWebGL2)` (loop bounds become number defines on WebGL1) and `getGrid3dComputeFrag(spec)` are exported; the template has no declarations and no `{UNIFORMS}` placeholder any more (also dropped the unused `uQuadShift`), everything is generated from the schema. Tests: 'grid3d compute' in `shader-uniforms.test.ts` (identifiers, no hand-written declarations, generated samplers, loop bound uniform only on WebGL2). Visual cases `geometry/mesh/grid3d-compute` (texture + loop bound) and `geometry/mesh/grid3d-compute-cumulative` (`tCumulativeSum`) compute the blob volume on the GPU; identical to `geometry/mesh/isosurface-cpu`, pass on webgl1 and webgl2-minimum.
7. **Validation** (done): unit tests above; visual tests on all profiles with `debug=1`.
   - New cases: `geometry/texture-mesh/gaussian-surface-gpu-smooth-colors` (GPU color smoothing, `sequence-id` theme) and `renderer/hi-z/occlusion` (new `VisualTestCase.frames`: Canvas3D frames drawn before the capture so interactive-only passes run; the captured image itself is not culled); grid3d cases from step 6.
   - Found and fixed: hi-z did not re-layout its levels when only `minLevel` changed (`copyTexSubImage2D` offset overflow); the harness missed GL errors in debug mode (`checkError` only logs them), now captured.
   - Result vs. `webgl2`: webgl1 75 pass; webgl2-no-float and webgl2-minimum 61 pass, 14 skipped; webgl1-minimal 44 pass, 30 skipped, 1 error: `geometry/image/slice` (pre-existing float `tValueTex` upload without `OES_texture_float`, now visible; resolved by making `OES_texture_float` and `EXT_blend_minmax` required on WebGL1, see [minimum-capabilities.md](minimum-capabilities.md)).
8. **Generated array declarations** (done; graphics and compute): the last hand-written uniform declarations were arrays, which were generated only as block members (WebGL2, sized by the spec's max `length`); WebGL1 has no blocks and used the hand-written declarations, sized by a define to stay within the small classic uniform budget (min. 16 fragment vectors, one per element).
   - Implemented: `UniformSpec(kind, variant, block, length, lengthDefine)`; `getShaderUniformDeclarations` emits classic arrays as `#if dX > 0` / `uniform T name[dX];` / `#endif`, or with the fixed `length`, and adds them to the types map; `addShaderDefines` runs the generated declarations through `preprocess` so ignored count defines (`dLightCount`, see `ignoreDefine`) are replaced like in the shader body. Tests: classic/block output in `uniform-buffer.test.ts`, `array length defines` check (number define of the same schema) and array uniforms in the "no hand-written declarations" guard (`shader-uniforms.test.ts`). Visual tests pass on webgl1, webgl1-minimal and webgl2-minimum (incl. clipping, two lights, ignore-light, SSAO multi-scale, shadows, bloom, direct volume).
   - Spec: optional length define, e.g. `UniformSpec('v3[]', undefined, undefined, MaxSamples, { lengthDefine: 'dNSamples' })`; neutral info, ignored by WebGL2 blocks and WGSL (max `length`). Validate that the define is a number define of the same schema.
   - Generator (`getShaderUniformDeclarations`, classic path): `#if dX > 0` / `uniform T name[dX];` / `#endif` (zero-length arrays don't compile); fixed `length` without a define (`uBloomFactors[5]`, `uBloomTints[5]`); add arrays to the types map so `stripUniformDeclarations` removes the hand-written ones (it handles the `[...]` suffix). Generated declarations follow the `#define` header, so the define is available (`dLightCount`/`dClipObjectCount` are also replaced textually, no conflict).
   - Arrays currently declared under conditions only (`uLevelRadius/Bias` under `dMultiScale`, lights under `dLightCount != 0`, clip objects under `dClipObjectCount != 0`) become unconditional apart from the size guard; unused uniforms are inactive and don't count against limits, the length define must always be set (check with the define-variant tests).
   - Shaders: `uKernel`, `uSamples`, `uLevelRadius/Bias`, `uLightDirection/Color` (light chunk, shadows, trace), `uClipObject*` (common vert/frag params chunks, direct-volume), `uGaussianCoefficients`, `uBloomFactors/Tints`.
   - Tests: extend the "no hand-written declarations" guard to arrays; strip/type checks and the WebGL1 varying/limit checks cover them; visual tests on webgl1 and webgl1-minimal (clipping, lights, SSAO multi-scale, bloom).
   - Result: no uniform, sampler or attribute declarations left in any shader except grid3d (step 6).

## Why it matters for WebGPU

The list is the complete inventory of shader/schema pairs, the input for generating WGSL bindings and for translating or porting the shaders, and it extends the minimum-capabilities checks to every program, not just the graphics ones.

## Open questions

- Helper vs. plain exported object: a helper could also own the `Values` type and the `createComputeRenderItem` call, but changes ~30 call sites; start with plain exports.
- Do any passes rely on hand-written declarations for define-dependent types not expressible in the spec (beyond `{ define }` dimensions and integer sample types)? Step 2 found none.
