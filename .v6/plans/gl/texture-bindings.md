# Plan: Generated texture/sampler declarations and bindings (WebGL now, WebGPU later)

## Status

- **Current**: done. T0–T2 and T4 done, T3.1 done, T3.2 dropped (decision 4), T3.3 not needed. T4 came with the WebGPU
  programs: `getWgslTextureDeclarations` from the binding layout (`webgpu/wgsl-declarations.ts`), per-stage limits in
  `_test/minimum-limits.test.ts`, read conventions in `webgpu/wgsl-conventions.md`.
- T3.1 done: WebGL texture units are a running index over the active samplers in the global, material and object tiers
  of `getBindingLayout` (decision 7, the order of the WebGPU bind groups; inactive samplers get no unit, so the units
  stay as low as the samplers allow), fixed per program; the sampler uniforms are set once after linking (the current
  program restored, the state cache untouched), `bindTextures` binds at the fixed units. WebGL images byte-identical on
  all WebGL profiles (`--compare`; `geometry/direct-volume/blobs` on `webgl1-minimal` varies between full runs before
  and after the change). Render benchmark: no measurable change (`render/boxes-1000` submit within the run-to-run
  variation; ANGLE takes the per-draw `uniform1i` calls cheaply).
- T3.3 not needed: `usedIf` on the remaining feature textures has no consumer. WebGPU declares and binds only the
  textures a program references (`isReferenced` in `getWgslTextureDeclarations`), WebGL skips inactive samplers, and the
  per-stage limits hold (`_test/minimum-limits.test.ts`).
- The survey (Today) and the steps describe the code before them; since: only compute and pass shaders still declare
  samplers by hand (stripped at runtime); units are fixed per program, global → material → object (no `startUnit`);
  T1.2 (`getTextureBinding`, capabilities in `getBindingLayout`) was not built, T4's generator is
  `getWgslTextureDeclarations(layout, isReferenced?)`; the T4 test passes without `usedIf`. Depth `sampleType` is still
  `depth` (`texture_depth_2d` in WGSL); decision 5 open on its own.
- Measured and dropped: a cache of the units' bindings in `WebGLState` (the active unit, the texture per unit; binds
  through it from `bindTextures` and `Texture.bind`, uploads and deletes invalidating it). WebGL binds every texture of
  every object on every draw (with 1000 unmerged objects on WebGL2 about a fifth of the frame), but consecutive objects
  bind their own textures to the same units: few binds repeat, the cache skipped little (`merging/units-1000-*`
  unchanged within noise) and its own bind took the cost. The cube textures binding on image loads without resetting
  the item bindings is harmless without such a cache (another target than the 2D units).
- Rule: must work on WebGL2 minimums and WebGPU default limits; anything needing more has a documented fallback (see
  [minimum-capabilities.md](minimum-capabilities.md)). Relevant here: 16 sampled textures per stage on both. Checked by
  `_test/minimum-limits.test.ts`: with worst-case feature defines every graphics program uses ≤ 16 per stage, and
  `getBindingLayout(schema, defines)` with the current `usedIf` annotations binds ≤ 16 in total, so T4.2 already holds
  for those define sets.
- T0 done: guard tests in `_test/shader-uniforms.test.ts` (`t*` identifiers, and hand-written sampler declarations of
  all graphics shaders vs. spec `dimension`/`sampleType`), bind-time `checkTextureBinding` in `webgl/program.ts` (debug
  mode), audit of `texture` kind specs (results in section 1a, decisions 4. and 5. revised).
- T2 done: `getShaderDeclarations` generates `uniform highp <sampler type> tX;` for every schema texture in both stages,
  with define-dependent dimensions resolved per program and integer sample types as `isampler*`/`usampler*` on WebGL2
  (`sampler*` on WebGL1, 3D omitted on WebGL1); hand-written declarations are stripped when their type is one of the
  spec's possible types and fail otherwise (`stripSamplerDeclarations`); program cache key hashes the generated code.
  Graphics shaders no longer hand-declare samplers (guard test); compute/pass shaders keep theirs (stripped at runtime,
  see [compute-shader-declarations.md](compute-shader-declarations.md)). Precision: `highp` everywhere (background and
  CAS used `mediump`). Depth `sampleType` is still `depth`, only integer-ness matters for GLSL; change it with T4.
- T1 done except capabilities: `TextureSpec` has `dimension` (incl. `{ define }` resolver, used by `tGridTex`,
  `tMinDistanceTex`), `sampleType` (explicit `sint` for the histogram pyramid / marching cubes textures, `cube` for
  `tSkybox`) and `usedIf` (`UsedIf`, `isUsed`, `DefineLookup` in `renderable/schema.ts`);
  `getBindingLayout(schema, defines?)` returns `TextureBinding`s with binding numbers (0 = uniform block). Debug mode
  checks dimension, integer sample type and `usedIf` of every active sampler against the linked program.
  Float-filterability downgrade (capabilities) not done yet. Tests: `_test/binding-layout.test.ts`.
- Follow-up of [uniform-buffers.md](uniform-buffers.md) (Phase 3 open item: "Textures/samplers are still
  hand-declared"). Builds on `getBindingLayout` and the WGSL generator (both present). Attributes follow the same
  approach ([attribute-bindings.md](attribute-bindings.md)); `usedIf` is shared between textures and attributes.
- Goal: the schema is the single source of truth for texture bindings too, so that WebGL declarations, WebGL texture
  units and WebGPU bind group layouts + WGSL declarations are all derived from it. Backend-specific code stays in
  `webgl/` and `webgpu/`; `uniform-layout.ts` (or a sibling `binding-layout.ts`) stays neutral.

## TL;DR

`TextureSpec` doesn't carry what a declaration needs: `'texture'` kind specs (externally created textures, 30+ of them)
don't say 2D/3D/cube, and their format/type/filter are only nominal. Add `dimension` and `sampleType` to `TextureSpec`
(defaults derived from `kind`/`format`/`dataType`, explicit for the exceptions), generate GLSL sampler declarations like
uniforms (strip + type-check hand-written ones, then migrate shaders), assign WebGL texture units from
`getBindingLayout` tiers, optionally use WebGL2 sampler objects so the spec's filter is authoritative, and emit WGSL
`texture_*` + a few shared `sampler`s. The hard WebGPU constraints are binding-count limits (graphics schemas have 18-21
textures vs. a default limit of 16 per stage) and WGSL's lack of a preprocessor, so bindings must be generated per
program (defines), not per schema.

## 1. Current state (facts)

- `TextureSpec(kind, format, dataType, filter, variant?)` (`packages/graphics/src/gl/renderable/schema.ts`); kinds
  (`webgl/texture.ts`): `image-uint8|float32|float16|int32|depth` (2D, created/uploaded by the render item),
  `volume-uint8|float32|float16` (3D, WebGL2), `texture` (a `Texture` object created elsewhere; spec format/type/filter
  are not enforced).
- Shaders hand-declare samplers (all shaders): 131x `sampler2D`, 3x `isampler2D`, 2x `sampler3D` (`tGridTex`,
  `tMinDistanceTex`), 1x `samplerCube` (`tSkybox`); precision via statements: 31x `precision highp sampler2D`, 2x
  `highp sampler3D`, 3x `highp isampler2D`, `mediump sampler2D` in background and CAS, `mediump samplerCube` in
  background; no inline precision on sampler declarations.
- Type depends on defines/backend in a few places:
  - `tGridTex` (direct-volume), `tMinDistanceTex` (gaussian density): `sampler2D` if `dGridTexType_2d`, `sampler3D` if
    `dGridTexType_3d`.
  - `tPreviousLevel`, `tTexture` (histogram pyramid), `tActiveVoxelsPyramid` (marching cubes): `isampler2D` on GLSL 300,
    `sampler2D` on GLSL 100.
  - `tSkybox` (background): `samplerCube`, `mediump`.
- Binding: `program.bindTextures(list, startUnit)` binds in list order; units = shared/global textures (renderer
  `sharedTexturesList`) first, then object textures, then material textures (`variant: 'material'`). Unused (null
  location) samplers are skipped.
- Sampler state: per texture object (`texParameteri`), wrap always `CLAMP_TO_EDGE`, min/mag from the filter passed at
  creation, mipmaps only for cube/skybox; no `gl.createSampler`.
- Reads: `texture2D`/`texture`, `textureCube`, `texture2DLodEXT`/`textureLod`, `texelFetch` (direct-volume, gaussian
  density), many in vertex shaders (`readFromTexture`) and inside non-uniform branches.
- Counts: `BaseSchema` has 14 textures + 4 in `GlobalTextureSchema`; spheres 20, texture-mesh 21. WebGL only counts
  active samplers per stage; WebGPU's default `maxSampledTexturesPerShaderStage`/`maxSamplersPerShaderStage` is 16 and
  counts every layout entry visible to the stage.
- Float linear filtering depends on `OES_texture_float_linear` (WebGL) / `float32-filterable` (WebGPU); code already
  falls back to `nearest` at texture creation (e.g. SSAO).

## 1a. Audit of `texture` kind specs (T0.4)

Bound texture vs. spec, collected with a temporary debug warning over all visual tests (incl. new skybox, image
background, gaussian surface GPU/CPU cases) on `webgl1`, `webgl2`, `webgl2-no-float`, `webgl1-minimal`. Dimension and
integer-ness never mismatched (checked by the debug checks). Everything else is only nominal:

- **Storage depends on capabilities**, so one spec can't describe it:
  - Depth (`tDepth` spec `depth`/`ushort`, pass `tDepth`/`tDepthOpaque` spec `rgba`/`ubyte`): `DEPTH_STENCIL` with
    `FLOAT_32_UNSIGNED_INT_24_8_REV` or `UNSIGNED_INT_24_8`, or packed `RGBA`/`UNSIGNED_BYTE` without depth textures
    (`depthTextureSupport` define).
  - `float` specs bound with `HALF_FLOAT` or `UNSIGNED_BYTE` storage on lower profiles: dpoit textures,
    `tMinDistanceTex`, histogram pyramid and marching cubes textures, `tNormal`, `tVolumeData` (spec `ubyte`, actual
    half float).
- **Format differs**: integer pyramid textures are `RED_INTEGER`/`INT` (spec `rgba`/`float`), `tGroup` `RGBA` (spec
  `alpha`), `tPosition`/`tNormal` `RGBA` (spec `rgb`), `tDpoitDepth` `RGBA` on WebGL1 (spec `rg`), SSAO
  `tDepthHalf`/`tDepthQuarter` `RED`/`FLOAT` (spec `rgba`/`ubyte`).
- **Filter differs**: `LINEAR` textures on `nearest` specs: `tColor` (copy, postprocessing, dof, illumination compose),
  `tTransparentColor`, `tSsaoDepth`, `tVolumeData`; `NEAREST` texture on a `linear` spec: SSAO `tDepth`. Generic shaders
  (copy) sample textures with different filters per use.

Consequences: `format`/`dataType`/`filter` of `texture` kind specs are documentation only and must not drive bindings;
only `dimension` and `sampleType` (resolved with capabilities) do. See decisions 4. and 5.

## 2. Design decisions

1. **Extend `TextureSpec`** with
   - `dimension: '2d' | '3d' | 'cube'` - default from kind (`image-*` → 2d, `volume-*` → 3d); for `texture` kind default
     `2d`, explicit otherwise (`tSkybox`: cube).
   - `sampleType: 'float' | 'unfilterable-float' | 'sint' | 'uint' | 'depth'` - default from `dataType`/`format` (`int`
     → sint, `depth` → depth, `float` + `linear` → float but downgraded at program creation without float filtering, see
     5.).
   - Keep `filter`; it becomes the sampler state for the binding (see 4.).
   - Signature: add an options object as 6th parameter
     (`TextureSpec(kind, format, dataType, filter, variant?, { dimension?, sampleType? })`) to stay compatible with
     existing callers.
2. **Define-dependent dimension** (`tGridTex`): allow `dimension` to be resolved per program, e.g.
   `dimension: { define: 'dGridTexType', map: { '2d': '2d', '3d': '3d' } }`. Programs are already per define set, and
   WebGPU pipelines/layouts are per program, so this needs no schema split. Alternative (simpler, more churn): split
   into `tGridTex2d`/`tGridTex3d`. Recommend the resolver since it's the only case and keeps values unchanged.
3. **Backend-dependent type** (`isampler2D` vs `sampler2D`): not a schema concern; the GLSL generator maps `sint`/`uint`
   to `sampler2D` on WebGL1 (current behavior) and to `isampler2D`/`usampler2D` on WebGL2.
4. **Sampler state from the bound texture, not the spec** (revised after the audit, 1a): making the spec filter
   authoritative (WebGL2 sampler objects) would change rendering for at least 7 textures and can't serve generic shaders
   like copy. WebGL keeps per-texture parameters. WebGPU picks one of the shared samplers (decision 8.) when creating
   the bind group, from the bound texture's filter, forced to a non-filtering sampler for
   `unfilterable-float`/`sint`/`uint`. The spec `filter` stays meaningful only for render-item created kinds (`image-*`,
   `volume-*`).
5. **Capabilities and storage at program creation**: `sampleType` is resolved with capabilities, as the audit shows
   storage varies per profile: 32-bit float storage + `linear` without float-linear filtering → `unfilterable-float` +
   nearest (WebGPU: non-filtering sampler); depth reads → `unfilterable-float` for both depth-stencil and packed RGBA
   storage (WebGPU allows binding the depth aspect as `unfilterable-float`, which keeps `texture_2d<f32>` and the GLSL
   `sampler2D` semantics; the current `depth` default for `depth` format specs should change accordingly when T2/T4
   consume it). Generator inputs are `(schema, defines, capabilities)`.
6. **Declarations per program, not per schema**: only textures the program variant uses are declared/bound. Derive usage
   from defines via an optional spec field (e.g. `usedIf: 'dOverpaint'` / predicate on define values), falling back to
   "always". Needed to stay under WebGPU binding limits and because WGSL has no `#ifdef`; for WebGL it also reduces unit
   usage. Until `usedIf` exists, GLSL declares all schema textures (unused samplers are inactive and cost nothing in
   WebGL).
7. **Tier binding numbers** (neutral, from `getBindingLayout`): per tier, binding 0 = uniform block, then textures in
   schema key order. WebGL texture units = running index over global → material → object tiers (same order as WebGPU
   groups; today object comes before material, harmless to change).
8. **WGSL conventions** (for whoever writes/translates WGSL shaders):
   - Types: `texture_2d<f32|i32|u32>`, `texture_3d<f32>`, `texture_cube<f32>`, `texture_depth_2d`.
   - Samplers: shared immutable samplers in group 0 (`sNearest`, `sLinear`, `sLinearMipmap`, all clamp-to-edge) instead
     of one per texture; `nearest` data textures don't need a sampler at all.
   - Reads: `textureLoad` for nearest/data/integer/depth textures (replaces `readFromTexture`'s UV math with integer
     coords), `textureSampleLevel` for filtered reads; never `textureSample` (implicit derivatives are invalid in
     non-uniform control flow and in vertex shaders, both common here).
   - Formats: WebGPU has no `alpha`/`rgb` 8-bit formats → upload `alpha` as `r8unorm` and `rgb` as `rgba8unorm`; shader
     code reading `.a` of alpha textures must read `.r` (texture resource concern, listed here because it affects shader
     code).

## 3. Phases

### T0 - Inventory and guards (no behavior change)

1. One-off scan (like the uniform scan): every hand-written sampler declaration vs. the spec it binds to; list
   define/version-dependent types (expected: the four above).
2. Extend the shader identifier guard test (re-add from the uniform plan): every `t[A-Z]\w*` identifier referenced in
   graphics shaders is a texture in the schema.
3. Debug-mode bind-time check in `program.bindTextures`: the bound `Texture`'s target/kind matches the spec's
   dimension/sampleType (catches `texture` kind mismatches, e.g. a 3D texture bound to a 2D spec).
4. Audit `texture` kind specs: spec `filter`/`format`/`dataType` vs. how the texture is actually created (prerequisite
   for 4. and 5.).

### T1 - Schema

1. Add `dimension`/`sampleType` (+ defaults) to `TextureSpec`, annotate `tSkybox`, `tGridTex` (define resolver), integer
   pyramid textures, `tDepth` (`depth`), float linear textures (`tValueTex`, grid textures if float).
2. Add neutral `getTextureBinding(spec, defines, capabilities)` → `{ dimension, sampleType, filter }` and extend
   `getBindingLayout(schema, defines?, capabilities?)` to return resolved texture bindings with binding numbers per
   tier.
3. Unit tests for defaults, resolver and binding numbering.

### T2 - Generated GLSL sampler declarations (`webgl/shader-declarations.ts`)

1. Generate `uniform highp sampler2D|sampler3D|samplerCube|isampler2D|usampler2D tX;` (inline precision; `sampler3D` has
   no default precision in GLSL ES 3.00 fragment shaders) for both backends; WebGL1 maps integer to `sampler2D` and
   rejects 3D/cube where unsupported.
2. Extend `stripUniformDeclarations` to samplers with the same type check (currently it keeps sampler lines untouched).
3. Program cache key already hashes generated code; include the resolved dimension defines once T1.2 lands.
4. Migrate graphics shaders/chunks, then compute passes (remove hand-written sampler declarations and
   `precision ... sampler*` statements that become redundant).
5. Validate: visual tests `webgl1` vs `webgl2` with `debug=1`, plus `webgl1-minimal`, `webgl2-no-float`.

### T3 - WebGL binding by tier

1. Texture units from the T1 binding layout instead of list order (renderer shared textures = global tier; render item
   material/object tiers).
2. ~~WebGL2 sampler objects (decision 4.)~~ dropped, see revised decision 4.
3. `usedIf` (decision 6.) for BaseSchema feature textures (overpaint, transparency, emissive, substance, clipping,
   wiggle, marker, color variants), with a debug check that a texture referenced by an active sampler is not excluded.

### T4 - WGSL declarations (`webgpu/`)

1. `getWgslTextureDeclarations(schema, defines, capabilities)`: `@group(tier) @binding(n) var tX: texture_*;` + shared
   samplers in group 0; report `unsupported` (e.g. cube without support, formats without WebGPU equivalent).
2. Tests: binding numbers match `getBindingLayout`, per-stage counts ≤ 16 for every graphics schema with representative
   define sets (fails today without `usedIf` → drives T3.3).
3. Document the WGSL read conventions (decision 8.) next to the generator.

## 4. Validation

- Unit tests: spec defaults/resolvers, binding layouts, GLSL/WGSL generator output, identifier guard for `t*` names,
  per-stage binding counts.
- Debug mode: bind-time texture/spec check, existing active-uniform checks (sampler types come from generated
  declarations).
- Visual tests on all profiles; cases exist for skybox and image backgrounds, direct-volume (2D grid texture on WebGL1,
  3D on WebGL2), GPU isosurface (integer pyramid textures), gaussian surface on GPU (`tMinDistanceTex`) and CPU.
  Missing: image with float values (unfilterable path on `webgl2-no-float`).

## 5. Open questions

- `usedIf` as define-name strings vs. predicates on define values (strings are serializable and testable; predicates are
  more flexible).
- Sampler precision: generate `highp` everywhere or keep `mediump` where shaders use it today (background)?
- Should `nearest` data textures move to `texelFetch` in GLSL 300 already (closer to WGSL `textureLoad`), or wait for
  WGSL?
- Combining small per-group data textures (marker, transparency, clipping, emissive: all `alpha`) into one multi-channel
  texture or a texture array to cut binding counts - bigger change to the data model, only if `usedIf` isn't enough.
