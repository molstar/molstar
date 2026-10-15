# Plan: Backend split (WebGL / WebGPU)

## Status

- Steps 1, 2 and 3 (3a-3f) done, step 4 done; steps 5-6 not started. Goal: code outside `gl/webgl/` (and later
  `gl/webgpu/`) depends only on backend-neutral interfaces, so a WebGPU backend can be added next to WebGL1/2 (see
  [minimum-capabilities.md](minimum-capabilities.md) for targets and the minimum-spec rule).
- Already neutral: renderable schema, `uniform-layout.ts`, `attribute-layout.ts`, `getBindingLayout`, the shader/schema
  inventory (`_test/graphics-shaders.ts`, `_test/compute-shaders.ts`); all GLSL declarations are generated from the
  schema ([compute-shader-declarations.md](compute-shader-declarations.md)).
- Still WebGL-specific: readback and fences (pick, hi-z, image, compute: step 5), `shader-code.ts` and GLSL 100
  workarounds (WGSL strategy plan); guard test allowlists.

## Step 1: backend-neutral capability object (done)

Done as planned: `gl/capabilities.ts` (`GraphicsFeatures`, `GraphicsCapabilities`), `gl/webgl/capabilities.ts`
(`getWebGLCapabilities(gl, extensions)`, half-float features also require `textureHalfFloat`), `gl/context.ts`
(`GraphicsContext`, `GraphicsResources`, `GraphicsLimits`); the `max*` getters were removed from `WebGLContext`
(breaking, recorded). All feature checks in the table below use capabilities; pass texture type choices
(`halfFloatRenderTarget ? 'fp16' : floatRenderTarget ? 'float32' : 'uint8'`) too. Guard test:
`gl/_test/backend-neutral.test.ts` (graphics, plugin, extensions; allowlist with reasons, fails on stale entries).
Visual tests: `CapabilityKeys` = `GraphicsFeatures` + derived, `getCapabilities(graphics)`.

Original plan: `GraphicsCapabilities` in `gl/capabilities.ts` (neutral), mapped per backend (`gl/webgl/capabilities.ts`:
`getWebGLCapabilities(gl, extensions, parameters)`), exposed as `capabilities` on a new, minimal `GraphicsContext`
interface (`gl/context.ts`, neutral) that `WebGLContext` implements and step 3 extends; recomputed on context restore
(after `resetExtensions`). Checks outside the backend folders take or read the `GraphicsContext` and use its
capabilities instead of extensions, `isWebGL2` or WebGL parameters.

- Features (named by what they enable, not by extension; same names as the visual test capability keys):
  - `texture3d` (sampling and rendering into slices; WebGL2), `flatInterpolation` (`flat` varyings; WebGL2),
    `integerTextures` and `mixedAttachmentFormats` (WebGL2; added in 3f)
  - `floatTextureLinear`, `halfFloatTexture`, `halfFloatTextureLinear`, `floatRenderTarget`, `halfFloatRenderTarget`,
    `depthTexture`
  - `fragDepth`, `multipleRenderTargets`, `textureLod`, `multiDraw`, `clipControl`, `depthClamp`
  - `backend: 'webgl' | 'webgl2' | 'webgpu'` for reports only, never for checks
- Limits are not capabilities: texture and render target size limits (`maxTextureSize`, `max3dTextureSize`,
  `maxRenderbufferSize`, `maxDrawBuffers`, `maxTextureImageUnits`) go to the resource interface as `resources.limits`
  (minimal neutral `GraphicsResources` with `limits` in this step, extended in step 3); `max3dTextureSize` is 0 without
  `texture3d`, `maxDrawBuffers` 0 without MRT.
- Not listed: required features (instancing, uint indices, derivatives, float textures, min/max blending), they are
  always available.
- Derived rendering features (wboit, dpoit, illumination) stay `isSupported` functions of the passes, computed from
  capabilities.

Checks to migrate (outside `gl/webgl/`):

| Check                                     | Today                                                              | Capability                                                         | Where                                                                                                                                                                            |
| ----------------------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Impostors                                 | `extensions.fragDepth`                                             | `fragDepth`                                                        | `checkSphereImpostorSupport`, `checkCylinderImpostorSupport` (`repr/structure/visual/util/common.ts`), element-sphere, bond-intra-unit-cylinder, particles spacefill, volume dot |
| Per-triangle groups without flat varyings | `!isWebGL2` → `Mesh.uniformTriangleGroup`                          | `flatInterpolation`                                                | blob-surface (particles), blob-surface-mesh (2x), gaussian-surface-mesh (2x), molecular-surface-mesh (2x), isosurface, segment; extension `sb-ncbr` tunnels                      |
| GPU surfaces                              | `colorBufferFloat`, `textureFloatLinear`, `drawBuffers`            | `floatRenderTarget`, `floatTextureLinear`, `multipleRenderTargets` | gaussian-surface-mesh, isosurface, segment `gpuSupport`                                                                                                                          |
| 3D vs. 2D textures                        | `isWebGL2`                                                         | `texture3d`                                                        | direct-volume (3d/2d), `geo/gaussian-density.ts`, `geo/gaussian-density/gpu.ts`                                                                                                  |
| Half-float data                           | `extensions.textureHalfFloat`                                      | `halfFloatTexture`                                                 | direct-volume (2x), isosurface                                                                                                                                                   |
| Texture size limits                       | `webgl.maxTextureSize`, `max3dTextureSize`                         | `resources.limits.*`                                               | gaussian-surface-mesh, `util/gaussian.ts`, direct-volume, isosurface, segment                                                                                                    |
| Passes                                    | `drawBuffers`, `depthTexture`, `isWebGL2`, `colorBufferFloat`      | `multipleRenderTargets`, `depthTexture`, ...                       | shader-manager, pick (MRT; async pick needs fences, a backend method), hi-z, dpoit/wboit/illumination `isSupported`                                                              |
| Plugin                                    | `extensions.colorBufferFloat`, `textureFloat`, `gl.drawingBuffer*` | `floatRenderTarget`; drawing buffer size from the context          | `packages/plugin/core/src/util/viewport-screenshot.ts`                                                                                                                           |

- Format choices inside passes (`isWebGL2 ? 'depth32f-stencil8' : ...`, `'alpha'` vs. `'rgba'` targets, integer textures
  in compute) are not feature checks; they move behind backend resource/render-target creation in step 3, not to
  capabilities.
- Backend-specific readback (`extensions/geo-export`: `isWebGL2 ? Uint8Array : Float32Array`) needs a neutral
  `readPixels` contract (step 3).
- Visual tests: `tests/src/visual/framework/capabilities.ts` derives its `Capabilities` from `capabilities` (plus the
  derived features) and its report limits from `resources.limits`; `floatTexture` and `blendMinMax` (now required) are
  already dropped from its keys and the cases' `requires`.
- Guard test: no `.extensions.`, `isWebGL2` or WebGL parameter access outside `gl/webgl/` except an explicit, shrinking
  allowlist (passes and compute until step 3).

## Step 2: neutral uniform buffer interface (done)

Done: `gl/uniform.ts` holds the `UniformBuffer` interface and the neutral uniform value types and helpers
(`UniformKind`, `UniformKindValue`, `UniformType`, `UniformValues`, `UniformsList`, `isArrayUniform`,
`isUniformValueScalar`, `cloneUniformValues`); `webgl/uniform.ts` keeps the GL setters and GLSL types,
`webgl/uniform-buffer.ts` the arena/ring buffer implementations.

Original plan: move the `UniformBuffer` interface (and writer types) out of `webgl/` next to `uniform-layout.ts`; the
WebGL arena/ring buffer implementations stay in `webgl/uniform-buffer.ts` ([uniform-buffers.md](uniform-buffers.md),
user TODO).

## Step 3: neutral context and resource interfaces

Split into sub-steps, each validated by type checks, unit tests, the guard tests and the visual tests:

- **3a neutral resource types (done):** `gl/texture.ts` (texture kinds/types/formats/filters/attachments, neutral
  `Texture` without GL fields and texture-unit binding, `Textures`/`TextureValues`, `getBytesPerElement`,
  `createNullTexture()`/`isNullTexture` (WeakSet marker, bound as an empty 2D texture by the WebGL program),
  `loadImageTexture`, `CubeFaces`), `gl/framebuffer.ts`, `gl/renderbuffer.ts`, `gl/render-target.ts`, `gl/buffer.ts`
  (array/attribute/elements types), `gl/timer.ts` (`GraphicsTimer`, `TimerResult`, `printTimerResults`),
  `gl/render-item.ts` (`DrawMode`, `Transparency`, render variants). `webgl/texture.ts` has `GLTexture extends Texture`
  (target/format/internalFormat/type/filter, `bind`/`unbind`); `WebGLResources.texture` returns it, WebGL-only code
  narrows neutral textures with a cast. Guard test `backend-neutral imports`: no `gl/webgl/*` imports outside the
  backend except pending modules (`context`, `render-item`, `program`, `compat`) and a per-file allowlist (`GLTexture`
  in compute and hi-z, `PixelPackBuffer`, `shader-code.ts`).
- **3b extended `GraphicsContext`/`GraphicsResources` (done):** `gl/context.ts` now has resource creation (`texture`,
  `cubeTexture`, `framebuffer`, `renderbuffer`, uniform buffers, `getByteCounts`, `linkPrograms`/`finalizePrograms`),
  render targets (`createRenderTarget`, `createDrawTarget`), named texture/framebuffer caches, drawing buffer
  (`bindDrawingBuffer`, `getDrawingBufferSize`), `pixelRatio`/`setPixelScale`, `timer`, `stats` (`GraphicsStats` in
  `gl/stats.ts`, was `WebGLStats`), context loss/restore, `xr` (`GraphicsXR`) and `destroy`.
  `WebGLContext`/`WebGLResources` only add WebGL members (`gl`, `extensions`, `state`, buffers, programs, shaders,
  vertex arrays, `namedComputeRenderables`, sync `readPixels`/fences, `clear`, `checkError`).
  `isWebGLContext(ctx): ctx is WebGLContext` in `webgl/context.ts`. Not in the neutral context yet: readback and fences
  (step 5), `clear` (step 4), compute renderable cache (done in 3d).
- **3c `VisualContext.graphics` (done):** `VisualContext.graphics` and `RepresentationContext.graphics` (renamed from
  `webgl`) are `GraphicsContext`, as are the visual callbacks (`mustRecreate`, `processValues`, `setOverpaint`/...
  smoothing), visual constructors and mesh color smoothing. WebGL-only GPU paths narrow with `isWebGLContext`: `useGpu`
  of gaussian surface, isosurface and segment also requires it (otherwise the CPU mesh is used), texture-mesh creation
  and texture-mesh color smoothing throw/skip without it, gaussian density volume requires it. Passes and canvas3d
  helpers stay on `WebGLContext`: they create scenes, renderers and render items (3d) and set GL state (step 4).
- **3d render items, programs and compute (done):** `gl/render-item.ts` has the neutral `GraphicsProgram` (`id`,
  `variant`, `isReady`, `link`, `finalize`), `RenderItem`/`GraphicsRenderItem`/`ComputeRenderItem` and
  `MultiDrawBaseData`; the WebGL `Program` extends `GraphicsProgram`, `WebGLRenderItem` narrows `getProgram`.
  `GraphicsContext` has `createGraphicsRenderItem`/`createComputeRenderItem` and `namedComputeRenderables`; all
  renderables, passes and compute create render items through it. Renderable constructors, `createRenderable`,
  `Scene.create` and the camera, handle, pointer, shader-manager and XR helpers take a `GraphicsContext`; direct-volume
  narrows with `isWebGLContext` for the GLSL 100 loop bound. Still on `WebGLContext` until step 4: renderer (casts to
  `Program` for uniforms and texture binding), passes, compute (GL state, `GLTexture` bind/format, readback), pick/ray
  helpers (create passes) and the debug registry/`Helper` (done in 3e).
- **3e external consumers (done):** debug-helpers (helper constructors, `DebugRegistry.ctx`, `Helper`), sb-ncbr tunnels
  (`graphics` param, was `webgl`), alpha-orbitals (`createSphericalCollocation(Density)Grid`;
  `canComputeGrid3dOnGPU(graphics?: GraphicsContext): graphics is WebGLContext` narrows to the WebGL GPU path) and
  geo-export (`RenderObjectExporter.add`, `AddMeshInput.graphics`, `getInterpolated*`) take a `GraphicsContext`;
  geo-export narrows at texture readback and throws without WebGL (step 5). `viewport-screenshot` uses
  `getDrawingBufferSize()` instead of `gl.drawingBuffer*`. Headless `screenshot` (`createContext`, `readPixels`) and the
  internal alpha-orbitals GPU compute stay WebGL-only.
- **3f format selection (done):** storage-only choices moved into the WebGL backend: depth texture types and depth
  renderbuffer formats are precision preferences (`float`/`float-stencil`, `depth32f`/`depth32f-stencil8`; WebGL1 falls
  back to `ushort`/`uint24-8`, `depth16`/`depth-stencil`), `alpha` render targets fall back to `rgba` on WebGL1 (value
  in `.r` either way), `Texture.mipmap()` is a no-op for non-power-of-two sizes on WebGL1,
  `resources.getValueChannel(texture)` tells where `alpha` format textures keep their value (`red` for WebGL2 float/int
  storage). Choices that change the data layout seen by shaders and readers stay with the consumers behind two new
  capabilities: `integerTextures` (int histogram pyramid and sum vs. packed RGBA) and `mixedAttachmentFormats` (MRT
  attachments with different formats and `rg`, vs. all `rgba` float32 on WebGL1: dpoit, tracing, texture-mesh color
  smoothing, GPU isosurface). Draw buffers follow with step 4, readback types with step 5, GLSL 100 workarounds (grid3d,
  direct-volume, shader-code) with the WGSL strategy plan.

Original scope:

- `GraphicsContext` interface (capabilities, resources, state-free render/compute submission, render targets, readback,
  timers, context loss), extended from the minimal step 1 interface and implemented by `WebGLContext`;
  `VisualContext.graphics` (renamed from `webgl`) and pass constructors take it.
- Neutral interfaces for program/pipeline, render item, texture, render target, framebuffer; format and filter selection
  moves into the backends (decision 4 of [texture-bindings.md](texture-bindings.md)).
- Compute renderables (`createComputeRenderItem`, `ComputeRenderable`) behind the same interface; the compute shader
  list is the inventory.
- Consumers outside the graphics package that import `WebGLContext` or `gl/webgl/*`: extensions alpha-orbitals,
  debug-helpers, geo-export (framebuffer + `readPixels`), sb-ncbr; plugin (`animation-loop` timer (done in 3a),
  `viewport-screenshot`, headless `screenshot` with `createContext`/`readPixels`); mesoscale-explorer (`Transparency`
  type, done in 3a). They move to the neutral types; headless (`gl` npm package) stays WebGL-only.

## Step 4: render state and pass model

Passes and compute set GL state imperatively (`state.enable/disable(gl.BLEND|DEPTH_TEST|SCISSOR_TEST|CULL_FACE)`,
`blendFunc`, `blendEquation`, `colorMask`, `depthMask`, `clearColor` + `gl.clear`; about 240 call sites in 24 files
under `canvas3d/passes`, `geo`, `gl/compute`). WebGPU bakes blend, depth, cull and color-write state into pipelines and
clears through render pass load ops.

Done:

- `gl/render-state.ts` (neutral, WebGPU names): `RenderState` with `blend` (`BlendState`: color/alpha `operation`,
  `srcFactor`, `dstFactor`; `min`/`max` with `one` factors), `colorWrite`, `depth` (`compare`, `write`; unset = no depth
  test or write), `stencil` (front/back compare and ops, masks, reference), `cull`, `frontFace`; unset fields reset to
  the defaults, so no state leaks between passes. Presets `RenderState.Default` (the collapsed "fullscreen quad, no
  blend, no depth" setup), `Additive`, `Premultiplied`, `Alpha`, `Constant`; `BlendState.Max` replaces
  `extensions.blendMinMax`. `ClearValues` (`color`, `depth`, `stencil`).
- `GraphicsContext`: `setRenderState`, `setBlendConstant` (dynamic, `constant` factors), `setViewport`,
  `setScissor`/`disableScissor` (dynamic; disabled = whole target), `clear(values)` (within the scissor rectangle,
  ignoring write masks, all color attachments get the same value, integer attachments included), `flush`.
  `WebGLContext.clear(r, g, b, a)` (unused) was removed.
- WebGL backend: `webgl/render-state.ts` maps the state onto the cached `WebGLState` setters (constant blend factors as
  `CONSTANT_COLOR`), restores the write masks of the current render state after a clear, and clears integer color
  attachments with `clearBufferiv` (attachment recorded per GL framebuffer by `GLTexture.attachFramebuffer`).
- `Framebuffer.setDrawBuffers(count)` replaces `extensions.drawBuffers` and is re-applied on context restore;
  `Texture.copyFromFramebuffer(framebuffer, x, y, width, height, dstX, dstY)` replaces `GLTexture.bind` +
  `copyTexSubImage2D` (hi-z, histogram pyramid).
- Renderer: per-pass state constants (opaque, back depth, emissive, blended, volume, WBOIT, DPOIT) instead of
  `enable`/`blendFunc`/`depthFunc` sequences; `drawObject` derives cull mode and front face per object (cached variants
  per pass state); solid-interior caps use stencil states and `clear({ stencil })`/`clear({ depth })`. The WBOIT/DPOIT
  blend states moved from the passes' `bind()` into the renderer, the background blend from `clear()` into `render()`.
  Behavior fixes from explicit state: the depth-back pass sets the front face per object (was leaked), temporal
  multi-sample accumulation always blends additively (was the leaked blend function), GPU gaussian density clears all
  channels (was alpha only, through the color mask).
- Migrated: all passes, `gl/compute` (copy, hi-z, histogram pyramid, marching cubes, grid3d), GPU gaussian density,
  texture-mesh color smoothing. Passes without readback take a `GraphicsContext` (draw, postprocessing, antialiasing,
  FXAA, SMAA, CAS, DOF, bloom, outline, shadow, SSAO, marking, multi-sample, illumination, tracing, WBOIT, DPOIT,
  background); pick, hi-z, image and compute keep `WebGLContext` for readback/fences (step 5).
- Guard tests: new `backend-neutral render state` (no `state.*` setters, `gl.clear*`/`copyTexSubImage2D`/draw buffer
  calls outside `gl/webgl/`, no allowlist); feature-check and import allowlists shrunk (dpoit, wboit, tracing, gaussian
  density, color smoothing, histogram reduction, isosurface removed; `GLTexture` imports gone from hi-z and reduction).
- Validation: type checks (full workspace), unit tests, guard tests; visual tests in debug mode on all five profiles are
  byte-identical to renders of the previous code (two webgl1 cases vary by one pixel between runs of the same code). New
  visual cases `renderer/clipping/solid-interior` and `transparency/{blended,wboit,dpoit,dpoit-4}/solid-interior` cover
  the stencil cap path, also identical to the previous code.

- Render pass grouping (done, follow-up): `GraphicsContext.beginRenderPass(target, { clear })` and `endRenderPass()`
  (`RenderPassTarget`: render target, framebuffer or `'drawing-buffer'`) replace
  `target.bind()`/`framebuffer.bind()`/`bindDrawingBuffer()` in neutral code; `bindDrawingBuffer` is WebGL-only now.
  Beginning a pass ends the open one; copies (`copyFromFramebuffer`), readback (`readPixels`), attachment changes
  (`attachFramebuffer`/`detachFramebuffer`, `setDrawBuffers`) and `flush` end it too, as WebGPU runs them outside of
  passes (on WebGL they also bind another framebuffer). Draws and `clear` need an open pass: `WebGLState.inRenderPass`,
  checked in debug mode by render items and `clear` (throws). Passes that changed attachments mid-flow now attach first
  and begin (again) afterwards, e.g. the emissive bloom target is cleared with its own depth, then drawn against the
  opaque depth texture in a second pass; volumes and helpers begin the target they used to inherit from the last bind.
  The begin clear and `clear` stay limited to the scissor rectangle (stereo and XR render several viewports into one
  target); a WebGPU backend turns a begin clear that covers the whole target into a load op, otherwise into a clear
  draw. Guard test `no target binds outside gl/webgl/` with an allowlist of WebGL readback binds (image, pick, density
  field, texture-mesh, `readTexture`, headless screenshot, geo-export; step 5). Visual tests in debug mode on all
  profiles: byte-identical, no render pass errors (verified to fire by removing a `beginRenderPass`).

- Global bindings (done, follow-up): neutral `GlobalBindings` (`gl/render-item.ts`: shared `textures`, `invalidate`,
  `destroy`) created with `GraphicsContext.createGlobalBindings(schema, uniforms, textures)`;
  `RenderItem.render(variant, globals?, mdb?)` and `GraphicsRenderable.render(variant, globals)` take it instead of a
  shared texture count. The WebGL implementation (`webgl/global-bindings.ts`) applies it in the render item before
  drawing: uses the program, writes the global uniform ring buffer on WebGL2 (classic uniforms otherwise, and for arrays
  and draw uniforms), binds the shared textures from unit 0, only when the program changed or after `invalidate` (which
  also resets the render item binding cache). The renderer calls `invalidate` per render method and for the
  solid-interior pass uniforms (instead of setting them on the current `Program`);
  `GraphicsContext.hasRenderPassStencil()` replaces its framebuffer binding query. The renderer takes a
  `GraphicsContext` and has no WebGL imports or checks left (feature and import allowlist entries removed). Visual tests
  in debug mode on all profiles byte-identical.

Not done here, follow-ups:

- Render item binding cache (done, follow-up): the manual `WebGLState.currentRenderItemId = -1` resets in hi-z,
  histogram pyramid, marching cubes, GPU gaussian density, texture-mesh color smoothing and `readAlphaTexture` are gone.
  The fast path for repeated draws of the same render item skips attributes/VAO, elements, front-buffer and material
  uniforms and material textures, so the WebGL backend now invalidates the cache where these can change:
  `RenderItem.update()` drops the item's cache (each reset preceded a render after `update()` or of a new item), texture
  `define`/`load`/`mipmap`/`reset` and the WebGL1 copy path reset the shared cache (they replace the texture of the
  active unit), `GlobalBindings.invalidate` as before. Guard test: no `state.current(Program|Material|RenderItem)Id`
  outside `gl/webgl/`. Visual tests in debug mode on all profiles byte-identical.
- Backend downcasts (done, follow-up): replace casts of neutral objects to WebGL implementations with type guards in the
  style of `isWebGLContext`, backed by a module-level `WeakSet` the WebGL factories register their objects in; the WebGL
  code throws on objects of another backend instead of misbehaving silently.
  - Guards: `isWebGLGlobalBindings` (`webgl/global-bindings.ts`), `isGLTexture` (`webgl/texture.ts`).
  - Casts to replace: `globals as WebGLGlobalBindings` (`webgl/render-item.ts`, render), `texture as GLTexture`
    (`webgl/program.ts` texture binding, `webgl/resources.ts` `getValueChannel`, `gl/compute/util.ts` readback type
    checks; the latter moves with the neutral readback, step 5).
  - Considered and rejected: narrowing the parameter types in the WebGL implementations (fails `strictFunctionTypes` for
    function-typed properties; method syntax only hides the cast), a backend type parameter on
    `RenderItem`/`GraphicsContext` (spreads through renderables and the renderer), opaque handles with a backend
    `WeakMap` to implementations (same safety, more restructuring).
  - Validation: type checks, unit tests (a foreign `GlobalBindings`/`Texture` throws), guard test pattern for the
    remaining casts (`as (GLTexture|WebGLGlobalBindings)` outside `gl/webgl/`).
  - Done: `registerGLTexture`/`isGLTexture`/`asGLTexture` in `webgl/texture.ts` (the resources register the objects they
    return, `wrap` returns a copy), `isWebGLGlobalBindings`/`asWebGLGlobalBindings` in `webgl/global-bindings.ts`; the
    `as*` helpers throw for objects of another backend. All listed casts replaced, no
    `as GLTexture`/`as WebGLGlobalBindings` left anywhere (the guard test pattern checks outside `gl/webgl/`). Unit
    tests in `context.test.ts` (foreign objects throw, WebGL objects narrow; the WebGL halves need a headless GL
    context). Visual tests in debug mode on all profiles byte-identical, every texture bind and draw goes through the
    guards.

## Step 5: synchronization and readback

WebGPU readback is async only (`mapAsync`), WebGL uses sync `readPixels` and fences in many places.

- Inventory: pick (fence + `clientWaitSync`), hi-z (fence + pixel pack buffer), histogram pyramid sum, grid3d and
  gaussian density (`readPixels`), `waitForGpuCommandsCompleteSync`, `ImagePass.getImageData`, geo-export, plugin
  screenshots.
- Neutral async readback (`read(target, rect) => Promise<TypedArray>`), sync variants WebGL-only; callers that need sync
  results (e.g. histogram sum in GPU isosurface) get an async path or a documented WebGL-only status.
- Goal: compute and GPU paths take a `GraphicsContext`; `WebGLContext` remains only where a WebGL-only path is intended
  (narrowed with `isWebGLContext`). Users of `WebGLContext` outside `gl/webgl/` and what ties them to it:
  - `gl.finish()` (sync): `gl/compute/histogram-pyramid/reduction.ts` (`createHistogramPyramid`, also calls the sum),
    `gl/compute/marching-cubes/active-voxels.ts`, `gl/compute/marching-cubes/isosurface.ts`.
  - sync `readPixels`/`waitForGpuCommandsCompleteSync`: `gl/compute/histogram-pyramid/sum.ts`, `gl/compute/grid3d.ts`
    (also the GLSL 100 schema, WGSL strategy plan), `geo/gaussian-density/gpu.ts`, `gl/compute/util.ts`,
    `geo/geometry/texture-mesh/texture-mesh.ts`.
  - fences, pixel pack buffers, readback: `canvas3d/passes/pick.ts`, `canvas3d/passes/hi-z.ts`,
    `canvas3d/passes/image.ts`, `canvas3d/canvas3d.ts`.
  - pass-through only (move with the functions they call): `geo/gaussian-density.ts`,
    `repr/structure/visual/util/gaussian.ts`, `canvas3d/passes/passes.ts`, `canvas3d/helper/pick-helper.ts`,
    `canvas3d/helper/ray-helper.ts`.
  - done: `geo/geometry/texture-mesh/color-smoothing.ts` takes a `GraphicsContext` (its `readPixels` calls were
    commented-out debug code); the color smoothing callers in `repr/visual.ts` and `gaussian-surface-mesh.ts` no longer
    narrow with `isWebGLContext`.
  - not step 5: `gl/renderable/direct-volume.ts` (`isWebGL2` for GLSL 100 loop bounds, WGSL strategy plan).
- Guard test (done): the blanket `context` entry of `PendingModules` in `backend-neutral.test.ts` is replaced by
  per-file entries in `ImportAllowlist` (reason and removing step per file, grouped by readback/sync, pass-through,
  `isWebGLContext` GPU paths, context creation; stale-entry check), so new `gl/webgl/context` imports fail and the list
  shrinks with this step.

- Passes on `GraphicsContext` (done, second WebGPU milestone sub-step 2): `Passes`, `ImagePass`, `PickPass`,
  `PickHelper`, `RayHelper`, `HiZPass`, `Canvas3D` and `Canvas3DContext` take a `GraphicsContext`;
  `GraphicsContext.readRenderTarget` (async, rgba8, bottom-left rectangle, rows top to bottom; WebGL flips) and
  `waitForGpuCommandsComplete`; `PickBuffers` is an interface with the WebGL implementation (`readPixels`, pixel pack
  buffers, fences) and an unsupported one for other backends, where picking finds nothing; `HiZPass` is unsupported on
  non-WebGL contexts; `Canvas3D.graphics`/`Canvas3DContext.graphics`, their `webgl` throws on other backends. WebGL
  visual tests byte-identical; on WebGPU `Canvas3DContext.create` and `Canvas3D.create` succeed (rendering needs the
  copy quad and passes).
- Instead of lazy pass internals, WebGPU render items (graphics and compute) resolve their WGSL program at the first
  render of a variant, so renderables of unported programs can be created and only rendering them fails; this removes
  the construction-time failures that motivated the lazy internals. Lazy internals remain an optional startup and memory
  optimization for both backends:
- Lazy pass internals (optional): `Passes` creates the draw, pick, multi-sample and illumination passes eagerly, and
  each pass its render targets and renderables (postprocessing, SSAO, outlines, bloom, OIT, ...) in its constructor.
  Passes create them on first use instead, so a backend only needs the programs of enabled features (WebGPU without
  postprocessing ports can run the draw pass), and startup time and memory drop on WebGL too. Visual tests
  byte-identical.

## Step 6: context creation and loss

- Detailed in [backend-selection.md](backend-selection.md): backend descriptors, async creation, runtime switching,
  builds with a subset of backends.

- Backend selection in `Canvas3DContext.fromCanvas` (WebGPU needs async `requestAdapter`/`requestDevice`, so creation
  becomes async or gets an async variant); WebGL stays the default.
- WebGPU device loss vs. WebGL context loss/restore (`handleContextRestored`, `resetExtensions`, resource `reset`)
  behind one neutral event.

## Related

- Shaders: [wgsl-strategy.md](wgsl-strategy.md) (hand-written WGSL evolving independently of GLSL, shared TS
  preprocessor; preprocessing of defines/includes/unrolled loops, hand-port vs. transpile vs. generate); inputs are
  `GraphicsShaders` and `ComputeShaders` with generated declarations. Program variants and the program cache map to
  pipeline permutations; async pipeline creation fits the existing parallel shader compile in `shader-manager`.
- WebGPU limit gaps: vertex buffers (vertex pulling), inter-stage variables for cylinders, WGSL texture declarations
  (T4) - see [minimum-capabilities.md](minimum-capabilities.md).
- Data arrays (group colors, sizes, markers, ...) in storage buffers instead of data textures on WebGPU core:
  [storage-buffers.md](storage-buffers.md).
- Primitive expansion (cylinders, lines, text quads/boxes) moved from the builders to the renderer, and WebGPU drawing
  optimizations (instanced triangle-strip quads, GPU-driven culling with indirect draws, render bundles):
  [primitive-expansion.md](primitive-expansion.md).
- Later: a type-only GPU compute interface in core so GPU-accelerated algorithms can live in `core/math` without
  importing graphics: [core-gpu-compute.md](core-gpu-compute.md). The neutral `GraphicsContext`/capabilities/async
  readback from this plan are its prerequisites.

## Order and validation

- Steps 1 and 2 are small and independent; 3-6 are incremental on WebGL, each keeping all visual test profiles green
  (`debug=1`) and extending the guard test allowlist downwards.
- First WebGPU milestone after step 4: a `webgpu` visual test profile rendering `geometry/mesh/primitives` (one graphics
  shader, no passes), compared against `webgl2`. Reached ([wgsl-strategy.md](wgsl-strategy.md) step 6): 14 mesh cases
  pass with 0 differing pixels.
- Second WebGPU milestone: `Canvas3D` and `ImagePass` on WebGPU with postprocessing off, the `webgpu` visual profile
  rendering through them (no Scene + Renderer shortcut in the harness), the 26 cases of the first milestone unchanged.
  Progress: sub-steps 1 (descriptors, partly), 2 (passes on `GraphicsContext`), 3 (`copy.wgsl` with a shared quad vertex
  chunk, verified by reading every `webgpu` case through a copy) and 4 (`ImagePass` readback through `readRenderTarget`,
  done with sub-step 2) are done; sub-step 5 done too: the `webgpu` profile renders through `Canvas3D` and `ImagePass`
  like the WebGL profiles (the Scene + Renderer shortcut and the harness-side list of unsupported props are gone,
  unported programs are skipped when a pass reaches them), 27 cases pass (the 26 plus `renderer/hi-z/occlusion`, hi-z
  being unsupported only disables occlusion culling). Milestone reached. Found on the way: depth textures bound where
  WebGL samples `.r` (SSAO depth copies) need `unfilterable-float` bindings, now the rule for textures read with
  `textureLoad`. Order: backend descriptors and async `Canvas3DContext.create`
  ([backend-selection.md](backend-selection.md) step 1); passes on `GraphicsContext` with lazy internals (step 5); the
  copy quad as the first WGSL compute program ([wgsl-strategy.md](wgsl-strategy.md) step 7); async readback in
  `ImagePass` via `readTexture` (step 5; picking stays WebGL-only until its readback is async); the harness on the real
  path. WebGL stays byte-identical at each step.
- Breaking API changes (`VisualContext.webgl`, `canvas3d.webgl`, exported WebGL types) go to
  [breaking-v6-changes.md](../breaking-v6-changes.md).

## Decisions

- The capability object lives on a new neutral `GraphicsContext` (not first on `WebGLContext`).
- Texture and render target size limits belong to the resource interface (`resources.limits`), not to the capabilities.
- Name: `GraphicsCapabilities`.
- `VisualContext.webgl` becomes `VisualContext.graphics` (typed `GraphicsContext`; breaking, recorded in
  [breaking-v6-changes.md](../breaking-v6-changes.md), done in 3c).
- Step 3 is split into sub-steps 3a-3f (format selection split off from 3e).
- WebGL-only paths (GPU compute until step 4) get a `WebGLContext` through a narrowing helper `isWebGLContext(ctx)`;
  other backends use the CPU path.
- The neutral `Texture` has no GL fields and no texture-unit `bind`/`unbind`; `GLTexture` in `webgl/texture.ts` extends
  it.
- Renderer globals are a backend object (`GlobalBindings`) applied by the render items, not uniforms set on the current
  program; the backend chooses uniform buffers, classic uniforms or bind groups.
- Render state is set immediately before draws (`setRenderState`, full state, unset = default) rather than stored on
  render items; WebGPU folds the current state into the pipeline key at draw time. Viewport, scissor, blend constant and
  stencil reference stay dynamic.
- Render passes are explicit begin/end calls on the context (not callbacks or pass objects): beginning a pass and
  operations that can't run inside one end the open pass, draws and clears require one. Begin clears are scissored like
  `clear`; the backend decides about load ops.

## Open questions

- None currently.
