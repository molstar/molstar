# Plan: Generated attribute declarations and locations (WebGL now, WebGPU later)

## Status

- Rule: must work on WebGL2 minimums and WebGPU default limits; anything needing more has a documented fallback (see
  [minimum-capabilities.md](minimum-capabilities.md)). Relevant here: 16 vertex attributes on both, 8 vertex buffers on
  WebGPU.
- Step 1 (guards) done: `a*` identifier guard in `_test/shader-uniforms.test.ts`, `checkAttributeBinding` (debug mode,
  now part of `bindAttributeBuffers` in `webgl/buffer.ts`; `AttributeBuffer` exposes `itemSize`). Findings:
  `aColor`/`aSize` are only referenced in the never-used `dColorType_attribute`/`dSizeType_attribute` branches (schema
  entries are commented TODOs), `aColorMode` is cylinders-only (`dDualColor`); all allowlisted.
- Step 2 (schema) done: `aVertex` is a `BaseSchema` entry with `usedIf: ['!noNonInstancedActiveAttribs']` (value from
  `BaseGeometry.createValues`/merged values, filled by the render item when used; no schema mutation); `AttributeSpec`
  has `usedIf` (`aSegment*`: `dSegmented`). The render item skips attributes that are unused for the context
  (`getContextDefineLookup`); debug mode checks `usedIf` of active attributes.
- Steps 3 and 5 done: neutral `getAttributeLayout(schema)` in `gl/attribute-layout.ts` (schema order, `mat4` = 4
  locations, location 0 = first per-vertex attribute without `usedIf`, else first without `usedIf`; cached per schema;
  hash in the program cache key). Programs call `gl.bindAttribLocation` for all layout entries before linking (works for
  GLSL 100 and 300, no shader changes); debug mode checks active attributes are at their layout location.
  `bindAttributeBuffers` binds every buffer at its fixed location, so one VAO per render item is shared by all variants
  (was one per variant) and no longer recreated on define changes; `getLocations` no longer queries attributes. Max size
  over all graphics schemas incl. `SegmentSchema` is 16 (merged cylinders), tested. Step 4 (generated declarations,
  `layout(location = N)`) open; not needed for fixed locations.
- Step 4 done: `getShaderDeclarations` (`webgl/shader-declarations.ts`) generates attribute declarations for the vertex
  shader from the attribute layout (`layout(location = N) in` on WebGL2, `attribute` on WebGL1, `bindAttribLocation`
  kept); hand-written declarations are stripped with a type check (`stripAttributeDeclarations`) and removed from all
  graphics shaders (guard test in `_test/shader-uniforms.test.ts`). The unused `aColor`/`aSize` branches were removed,
  see [attribute-color-size.md](attribute-color-size.md). Remaining hand-written: compute shaders (stripped at runtime,
  see [compute-shader-declarations.md](compute-shader-declarations.md)).
- Fixed: without base-instance draw extensions, `offsetAttributes` changed instance attribute offsets in the VAO during
  multi-draw and didn't reset them, so a later draw without `mdbDataList` (culling turned off) used stale offsets; the
  render item now resets them to 0 after the multi-draw loop.
- Split out of [uniform-buffers.md](uniform-buffers.md) (formerly its Phase 4). Sibling of
  [texture-bindings.md](texture-bindings.md); `usedIf` is shared between textures and attributes.
- Goal: like uniforms (uniform plan, Phase 3) and textures, attribute declarations and locations come from the schema
  instead of hand-written `attribute ...` lines + `getAttribLocation` after link (see `getLocations` in
  `packages/graphics/src/gl/webgl/program.ts`). Backend-specific code stays in `webgl/` and `webgpu/`.

## 1. Current state

- `AttributeSpec(kind, itemSize, divisor)`: `kind` is always `'float32'`, `itemSize` 1/2/3/4/16 (16 = `mat4`, bound as 4
  consecutive `vec4` locations), `divisor` 0 (per vertex) or 1 (per instance). ~57 hand-written `attribute`
  declarations; some conditional (`aColor` for `dColorType_attribute`, `aSize` for `dSizeType_attribute`, `aSegment*`
  for `dSegmented`, `aVertex` without `gl_VertexID`).
- Locations are assigned by the driver and queried after link; VAOs are created per program variant because locations
  may differ between variants. `aVertex` is added to the schema object at render item creation (schema mutation).
- One GL buffer per attribute (no interleaving). Per-instance data: `aTransform` (mat4), `aInstance`, `aSegment*`;
  everything else per vertex.
- Location counts (worst case per graphics schema, e.g. cylinders): `aTransform` 4 + `aInstance` + `aGroup` + 6 cylinder
  attributes + `aColor` + `aSize` + `aVertex` + 3 `aSegment*` ≈ 18; never all active at once. Buffers for plain
  cylinders: 9.

## 2. Evaluation: is the schema the right source of truth for a WebGL + WebGPU backend?

Yes for the CPU↔GPU interface (attributes, uniforms, textures), no for shader-internal interfaces (varyings, fragment
outputs).

- **Required for WebGPU anyway**: WebGPU has no reflection (`getAttribLocation`/`getUniformLocation`/active-uniform
  queries). Pipelines need explicit `GPUVertexBufferLayout`s (format, offset, stride, step mode, shader location) and
  bind group layouts at creation time, so this information must exist on the CPU side. The schema already has it (kind,
  item size, divisor; uniform kinds/blocks; texture specs) and is the only place that knows about the data, so deriving
  both backends from it avoids a second, parallel description.
- **Separate logical from physical**: the schema describes _what_ (name, element type, components,
  per-vertex/per-instance, optional `usedIf`); each backend decides _how_ (WebGL: attribute + location + VAO; WebGPU:
  vertex buffer slot/interleaving or storage buffer + vertex pulling). This keeps the schema backend-neutral and lets
  the backends differ where their limits differ.
- **WebGPU limits force a different transport**: default `maxVertexBuffers` is 8 and `maxVertexAttributes` 16; plain
  cylinders already use 9 buffers, and `mat4` vertex inputs don't exist in WGSL (4 × `vec4<f32>`). So WebGPU needs
  either interleaving into few buffers or, more robustly, per-instance data (transforms, instance ids, segment data) in
  storage buffers read by `instance_index` (vertex pulling; similar to what Mol* already does with data textures). A
  schema-driven generator can make that choice per backend without touching the shared description; hand-written
  declarations can't.
- **WebGL benefits on its own**: explicit, stable locations (`layout(location = N)` on GLSL 300, `bindAttribLocation`
  before link on GLSL 100) make locations identical across program variants → one VAO per render item instead of per
  variant, no post-link attribute queries (helps `KHR_parallel_shader_compile`), deterministic location 0 (should be an
  always-active per-vertex attribute; some drivers penalize a disabled attribute 0).
- **Costs/risks**: shaders no longer show their inputs (mitigate: debug dump of the generated header, guard tests, keep
  names unchanged); conditional attributes need `usedIf` (same mechanism as textures) or are declared unconditionally
  (fine in WebGL, inactive attributes don't count); the `aVertex` schema mutation conflicts with per-schema caches and
  should become a regular schema entry with a `usedIf`-style condition.
- **Not in the schema**: varyings (~113 hand-written declarations) are internal to a shader pair; they belong to the
  shader code. For WGSL they become an inter-stage struct (watch `maxInterStageShaderVariables`, default 16). Fragment
  outputs already live in `ShaderCode.outTypes`/`drawBuffers`. Defines stay `DefineSpec`s in the schema; WGSL has no
  preprocessor, so WGSL sources need a preprocessing step (or override constants) either way.

## 3. Steps

1. **Guards**: extend the identifier guard test to `a[A-Z]\w*` (every referenced attribute is an attribute in the
   schema); debug check that bound attribute buffers match the spec (item size, divisor).
2. **Schema**: make `aVertex` a regular entry (no schema mutation); add an optional `usedIf` for conditional attributes
   (shared with textures); document that `kind` + `itemSize` + `divisor` is the complete logical description (add
   `normalized`/integer kinds only when needed).
3. **Neutral attribute layout** (`uniform-layout.ts` sibling, e.g. `attribute-layout.ts`): ordered list of
   `{ name, kind, itemSize, divisor, locations }` with deterministic locations (schema order, `mat4` → 4 locations, an
   always-used per-vertex attribute at location 0), cached per schema.
4. **WebGL generation** (`webgl/`): generate `attribute`/`in` declarations (`layout(location = N)` on WebGL2,
   `gl.bindAttribLocation` before `linkProgram` on WebGL1), strip + type-check hand-written ones (extend
   `stripUniformDeclarations`), migrate graphics shaders, then compute shaders (`aPosition` of the quad).
5. **WebGL VAOs**: use layout locations in `bindAttributes`/`offsetAttributes` and share one VAO per render item across
   variants; keep `getLocations` only for classic uniforms/textures.
6. **WebGPU generation** (`webgpu/`): vertex input struct with `@location(n)`, `GPUVertexBufferLayout`s from the
   attribute layout (`mat4` split into 4 `vec4`); decide the per-instance transport (interleaved instance buffer vs.
   storage buffer + `instance_index`) and enforce `maxVertexBuffers`/`maxVertexAttributes` in tests for every graphics
   schema.
7. **Validation**: visual tests on all profiles (`debug=1`), unit tests for location assignment and the generated
   declarations/layouts.

## 4. Open questions

- Per-instance transport on WebGPU: interleaved instance vertex buffer vs. storage buffer + vertex pulling (the latter
  also lifts the attribute limits, but changes how instancing and `MultiDrawBaseData` base instances are handled).
- Should WebGL2 also move per-instance data to textures (vertex pulling via `texelFetch`) to keep both backends
  structurally closer, or keep attributes there?
- Location 0 choice per geometry (`aPosition`/`aGroup`/`aVertex`) when the obvious candidate is conditional.
