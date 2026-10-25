# WGSL conventions

Rules for WGSL shader sources and the WGSL declaration generators in this folder. WGSL sources are written for WebGPU
and evolve independently of the GLSL ones; the shared contract is the renderable schema (values, defines, uniforms,
textures, attributes) and rendering equivalent images where both backends support a feature.

## Sources and composition

- One module per program, both stages in one source: `shader/<name>.wgsl.ts` exporting `<name>_wgsl` (e.g. `mesh_wgsl`),
  chunks in `shader/chunks/<name>.wgsl.ts` exporting `<name>_wgsl`.
- Chunks are composed in TS: a module imports the chunks it uses and interpolates them (`${common_wgsl}`) or calls chunk
  functions for structural choices; there is no `#include`.
- Template literals with WGSL source are tagged `/* wgsl */` (`export const mesh_wgsl = /* wgsl */ `...``) for editor
  highlighting.
- Entry points: `vertexMain`, `fragmentMain`, `computeMain`.
- No hand-written declarations of uniforms, textures, samplers, storage buffers or vertex inputs: they are generated
  from the schema (see Bindings, Vertex inputs). Names stay the schema names (`uView`, `tColor`, `aTransform`), so the
  same schema drives both backends.

## Preprocessing and defines

- Sources go through the shared shader preprocessor with macro expansion
  (`preprocessShader(source, defines, { expandMacros: true })`); its strict subset applies (`#if`/`#elif` with
  `defined`, `!`, `&&`, `||`, parentheses and integer comparisons, `#ifdef`, `#ifndef`, `#else`, `#endif`, `#define`,
  `#undef`).
- Defines have the flattened schema names, as in GLSL: string defines as `dX_value`, booleans as `dX`, numbers as `dX`
  with their value.
- Prefer WGSL functions to function-like macros; keep `#if` for small switches inside chunks and TS composition for
  structural choices (which chunks, outputs, entry points).
- A program binds every texture its source names (layouts come from the referenced names), so a feature a variant turns
  off is left out of that variant's source, not only out of its calls: `WgslProgram.source` may be a function of the
  variant's define lookup (`resolveWgslSource`) that includes the feature's helpers and `read_<name>` functions only
  when the feature is on (the compose passes, the marking overlay, SSAO). Checked on Dawn (`webgpu-unused-bindings`): a
  pipeline without any one of its bindings must fail validation.
- No `__VERSION__`, extension or precision switches: WGSL has one language version and `f32` everywhere.

## Override constants

- Number defines that don't change declarations (loop bounds over arrays sized by their maximum length, sample counts,
  thresholds) are read through `override` constants, so one shader module serves several pipelines. Defines that select
  code paths, bindings or array sizes stay preprocessor defines.
- Override names replace the `d` prefix with `o` (`dLightCount` → `override oLightCount: i32;`), so they never collide
  with the preprocessor macro of the same define.
- The backend sets the values of the overrides a module declares (`constants` of the pipeline stage) from the define
  values; overrides without a define keep their default.

## Bindings

Binding numbers of the uniform blocks and textures come from `getBindingLayout` (`uniform-layout.ts`), the same layout
the WebGL backend uses for texture units (block at binding 0, textures from binding 1 in schema key order). Samplers and
storage buffers follow the textures of their group; they are not part of that layout yet.

| Group | Contents                                                                                                          |
| ----- | ----------------------------------------------------------------------------------------------------------------- |
| 0     | global: binding 0 `uGlobal: Global` (uniform), then global textures in schema key order, then the shared samplers |
| 1     | material: binding 0 `uMaterial: Material`, then material textures                                                 |
| 2     | object: binding 0 `uObject: Object`, then object textures, then storage buffers of the data arrays                |
| 3     | draw: binding 0 `uDraw: Draw` (uniform with dynamic offset, written per draw)                                     |

- Uniform blocks: structs and `var<uniform>` from `getWgslUniformDeclarations`; member layouts equal std140. `bool`
  members are `u32` (non-zero is true); arrays of scalars and 2-component vectors are padded to 4-component vectors
  (`array<vec4<f32>, N>`, element `i` read as `a[i].x`).
- Per-draw values (WebGL classic uniforms with the `draw` scope: `uDrawId`, `uSolidInteriorPass`, `uSolidInteriorClip`)
  live in the `Draw` struct of group 3. Every other uniform must belong to a block (`global`, `material`, `object`);
  WebGPU has no classic uniforms, the generator reports the rest as unsupported.
- Samplers: shared immutable samplers in group 0, `sNearest`, `sLinear` and `sLinearMipmap`, all clamp-to-edge; textures
  read with `textureLoad` need none.
- Four groups is the WebGPU default `maxBindGroups`; nothing goes beyond group 3.

## Vertex inputs and builtins

- A generated `VertexInput` struct with `@location(n)` per attribute, using the same locations as the WebGL attribute
  layout; `mat4` attributes are split into four `vec4<f32>` inputs (`aTransform0` to `aTransform3`) and rebuilt in the
  shader.
- `@builtin(vertex_index)` and `@builtin(instance_index)` replace the vertex id emulation and `gl_InstanceID`; convert
  to `i32` where the GLSL code used `int`.
- Data that doesn't fit the vertex buffer limits (8 buffers, 16 attributes) is pulled from storage buffers by index
  instead of being a vertex input.

## Inter-stage variables

- A `VertexOutput` struct: `@builtin(position) position: vec4<f32>`, then the varyings with `@location(n)`, named like
  the GLSL varyings (`vColor`, `vGroup`); flat ones with `@interpolate(flat)`. The fragment stage takes the same struct.
- At most 16 locations (`maxInterStageShaderVariables`); pack scalar varyings into `vec4`s where a program would exceed
  it.

## Fragment outputs

- A `FragmentOutput` struct with `@location(i)` per color attachment (`vec4<f32>`, or `vec4<i32>` for integer targets)
  and `@builtin(frag_depth)` when the program writes depth (impostors).
- Multiple render targets are just more locations; there is no `gl_FragData`.

## Coordinates and depth

- Clip space depth is `[0, 1]`: the WebGPU backend builds the projection for it, so `uProjection` already maps to
  `[0, 1]`; code that linearizes or reconstructs depth uses the `[0, 1]` formulas.
- `@builtin(position)` in the fragment stage is in pixels with the origin top-left (y down) and pixel centers at `.5`.
  Render targets also store row 0 at the top, so screen-space texture coordinates computed from it
  (`position.xy / size`) need no flip; only values compared with GL-style bottom-left coordinates do.

## Textures

- Types: `texture_2d<f32|i32|u32>`, `texture_3d<f32>`, `texture_cube<f32>`, `texture_depth_2d`.
- Reads: `textureLoad` for nearest, data, integer and depth textures (integer texel coordinates, e.g. index `i` of a
  data texture at `(i % width, i / width)`); `textureSampleLevel` for filtered reads; never `textureSample` (its
  implicit derivatives are invalid in vertex shaders and in non-uniform control flow, both common here).
- Bind group layouts declare textures only read with `textureLoad` as `unfilterable-float`, which accepts float, unorm,
  unfilterable float and depth formats (e.g. a depth texture bound where WebGL samples `.r`); textures sampled with a
  filtering sampler (`textureSample*`) are `float`.
- Filtering is a property of the bound texture, as on WebGL: the render item defines `dLinear_<texture>` for bound
  textures that are linear and of a filterable format (part of the program key), and `texture_read_wgsl(name)` generates
  `read_<name>(uv)`, a filtered read with it and the nearest texel (`textureLoad`) without.
- Quad (compute) programs that port GLSL window-space math keep it: `glCoords` gives the GL fragment uv (bottom-up),
  `read_<name>` takes such a uv and flips it to the texture rows, offsets and bounds stay in GL space; their own GL
  projection uniforms keep the `[-1, 1]` formulas (only the global block holds `[0, 1]` matrices).
- Formats: there are no 8-bit `alpha`/`rgb` formats; `alpha` textures are stored as `r8unorm` and read from `.r`, `rgb`
  as `rgba8unorm`.

## Control flow and loops

- Loops over light and clip object arrays are plain loops bounded by override constants
  (`for (var i = 0; i < oLightCount; i++)`); there is no unroll pragma (WGSL allows dynamic bounds and dynamic indexing
  of uniform arrays, and constant bounds can be unrolled by the compiler).
- Derivatives (`dpdx`, `dpdy`, `fwidth`) only in uniform control flow; `discard` is fine.
- No implicit conversions: write `f32(i)`, `i32(x)`; use `select(f, t, cond)` where GLSL used a ternary on values.
