# v6 compatibility changes

Record public API and behavior changes introduced by the refactor here. The current prototype audit is complete for the
scope below. Update this ledger as implementation continues; release/downstream validation remains in
[checklist.md](checklist.md).

## Audit checkpoint

Compared v5 master `a1bf00ba55d737461e7e149ccc1996d979523b34` with prototype source at
`094c412a2f40692440434b5cd6eea8f1d0b77821` on 2026-10-04. Using the migration map, compared 1,547 production
TypeScript/TSX modules across packages, extensions, apps, CLI tools, and servers. No mapped targets or upstream
production source mappings were missing. Tests/performance fixtures were outside this source scan.

A TypeScript syntax/token scan removed top-level imports, normalized export module paths, and ignored
comments/formatting. It identified 63 remaining source differences, reviewed for exported names, namespace/class
members, signatures, and behavior. Reviewed the 17 new/extracted source modules separately. Internal shape-factory call
changes, ESM inline import suffixes, alias-only type changes, and CLI path rewrites account for much of the remaining
diff.

Checked compiled public modules and actual classic Viewer globals in local Chrome. Node probes confirmed the MVS
cache/file-URL changes and the standalone/runtime MolQL validation distinction. This checkpoint records the API/behavior
diff; it does not certify arbitrary downstream applications, native GL, or browser platforms. Changes already present in
the merged v5 master, including rendering/marking and BinaryCIF/PDB fixes, are not attributed to this refactor.

## Preserved classic globals

The classic Viewer still exposes `molstar.Viewer`, `ViewerAutoPreset`, `PluginExtensions`, `ExtensionMap`, `lib`,
`version`, `consoleStats`, `isDebugMode`, `isProductionMode`, `isTimingMode`, `setDebugMode`, `setProductionMode`, and
`setTimingMode`.

`molstar.lib` retains `structure`, `volume`, `shape`, `loci`, `math`, `plugin`, and `extensions`. The Viewer and MVS
Stories entrypoint implementations retain their public exports. MVS Stories still registers its custom elements and
exports its context/load/download functions and `MVSData`. Retaining these namespaces does not preserve every nested
API, particularly the shape changes below.

## Accepted: classic `molstar.lib.shape`

The classic bundle still exposes `molstar.lib`, but its shape objects now expose the model-side API following the
separation of model data and graphics helpers.

- `Shape.createRenderObject`, `Shape.createTransform`, `Shape.getTheme`, and `Shape.groupIterator` are no longer on the
  classic model-side `Shape` object.
- `ShapeGroup.getBoundingSphere` is no longer on its model-side object.
- Model-side `Shape.create` requires `groupCount` as argument 7, with optional `transforms` as argument 8. The old
  seventh-argument transforms call must change. JavaScript callers passing the old signature can silently store the
  transforms array as `groupCount`.

Graphics integrations belong to `@molstar/graphics/geo/shape/shape`. These classic API changes are accepted for v6; do
not restore them through a compatibility shim.

The graphics module exports functions directly, rather than a replacement `Shape`/`ShapeGroup` namespace:

| v5 call                                       | v6 graphics call                        |
| --------------------------------------------- | --------------------------------------- |
| `Shape.create(..., transforms?, groupCount?)` | `create(..., transforms?, groupCount?)` |
| `Shape.createRenderObject(shape, props)`      | `createRenderObject(shape, props)`      |
| `Shape.createTransform(...)`                  | `createTransform(...)`                  |
| `Shape.getTheme(shape)`                       | `getTheme(shape)`                       |
| `Shape.groupIterator(shape)`                  | `groupIterator(shape)`                  |
| `ShapeGroup.getBoundingSphere(loci, out?)`    | `getGroupBoundingSphere(loci, out?)`    |

For rendered shapes, import that graphics module as a namespace and use its `create` function; it preserves the old
optional-argument order and computes the initial group count from geometry when omitted. Model-side `Shape` defaults its
generic geometry type to `ShapeGeometry`, not the full graphics `Geometry` union. Supply a concrete generic type when
accessing geometry-specific fields.

`shape.groupCount` is now a stored value rather than a getter deriving the count from geometry on each access. Consumers
changing geometry group counts after creation must recreate the shape with the correct count instead of relying on the
old getter. The model contract still marks `groupCount` readonly.

## Standalone MolQL validation

The standalone builder and `mvs-validate` validate MolQL syntax, callable names, and argument names against
`@molstar/query-language` symbol tables without loading the molecular query runtime. Unknown callable names fail
standalone validation by default; custom definitions can be supplied to `MVSData.validationIssues`/`isValid` through
`getMolQLSymbol`. Runtime validation also checks executable compilation.

`@molstar/mvs-builder/molql` exposes `MolScriptBuilder` and `compileScript` for MolScript, PyMOL, VMD, and Jmol text.
Language-only imports moved to `@molstar/query-language`; executable compiler/runtime functionality remains under
`@molstar/model/script/...`. See [the detailed migration notes](#standalone-query-language-and-mvs-molql-validation) and
the [validation design](../designs/architecture.md#72-molql-builder-and-validation-design).

The `@molstar/mvs/behavior-id` constants keep the existing `molviewspec` name and `ms-plugin.molviewspec` transformer
ID; the loader checks that stable ID without importing the behavior module. Serialized snapshot identity is unchanged.

## Video export (breaking)

The MP4 extension is renamed from `@molstar/mp4-export-extension` to `@molstar/video-export-extension`
(`packages/extensions/video-export`). Viewer and Mesoscale Explorer extension options use `video-export` instead of `mp4-export`.
Public names change from `Mp4Export`, `Mp4Controls`, `Mp4EncoderUI`, `Mp4EncoderParams`, and `Mp4HeadlessPluginContext`
to their `Video` equivalents. `Mp4AnimationParams` is replaced by backend-generated `PD.Params` exposed through
`VideoControls.behaviors.paramsSchema`. The serialized behavior ID `ms-plugin.extension-mp4-export` is retained so
existing snapshots remain loadable.

`encodeMp4Animation` becomes `encodeVideoAnimation`; it and headless `getAnimation` return
`{ data, mimeType, extension }` rather than raw bytes. `fps` remains a shared option; codec settings live in
`backendParams`. The renderer produces timestamped frames at the requested FPS without the former extra terminal frame.
`saveAnimation` validates the filename extension and propagates filesystem errors.

Mediabunny is the default backend. Codec/container choices and codec-specific parameters are generated for the current
runtime, resolution and FPS. Discovery is shared across formats and includes quality fallbacks and explicitly selected
codec settings. Quality/FPS presets preserve the selected format and codec; Automatic hardware selection remains the
default. Programmatic exports reject invalid or unknown settings before normalization can substitute defaults.
Unsupported settings fail before animation playback. ProRes encoding requires a registered custom/server encoder; it is
not provided by browser WebCodecs. No native encoder dependencies are installed automatically.

`HeadlessPluginContext` owns image rendering only. Use `VideoHeadlessPluginContext` from
`@molstar/video-export-extension/headless` and register `VideoExport` for animation rendering. For Node environments
without WebCodecs, explicitly import `H264WasmBackend` from `@molstar/video-export-extension/backends/h264-wasm` and
pass `{ backend: H264WasmBackend, backendParams: { quantization: 18 } }` to `getAnimation`/`saveAnimation`. This
replaces the former top-level `quantization`/`quantizationParameter` options. The MVS rendering CLI explicitly uses this
backend for MP4 output. Browser entry points do not import the WASM encoder. Native GL remains opt-in.

## Shape bounds and graphics helpers

Model-side shapes use a structural `ShapeGeometry` contract. Geometry-aware shape creation and group bounds belong to
`@molstar/graphics/geo/shape/shape`. `Loci.getBoundingSphere` for shape groups delegates to the optional
`shape.getGroupBoundingSphere` callback; the graphics factory supplies it. Shapes created directly through the model
factory need that callback for group bounds.

Streamline picking and location-iterator helpers moved from model props to `@molstar/graphics/geo/streamlines`:
`getStreamlinesVisualLoci`, `getStreamlinesLoci`, `eachStreamlines`, and `createStreamlinesLocationIterator`.

`PositionLocation` and `isPositionLocation` are owned by the model package and re-exported by the graphics location
iterator. `Viewport` is owned by core and re-exported by the graphics camera utility. Those re-exports preserve existing
behavior at the mapped module paths.

## Other helper moves requiring new imports

These are symbol-level moves within formerly larger modules; the whole-file migration map alone is insufficient for
these imports. Symbol moves from the plugin composition split are recorded in `migration-symbols.json` and described in
_Plugin composition step 1_ below.

| Public API                                                                                                                          | v5 module                            | v6 import                                              |
| ----------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ | ------------------------------------------------------ |
| `CollapsableControls`, `CollapsableProps`, `CollapsableState`                                                                       | `mol-plugin-ui/base`                 | `@molstar/plugin-ui/controls/collapsable`              |
| `GaussianDensityTextureData`, `computeGaussianDensityTexture`, `computeGaussianDensityTexture2d`, `computeGaussianDensityTexture3d` | `mol-math/geometry/gaussian-density` | `@molstar/graphics/geo/gaussian-density`               |
| `fillIdentityTransform`                                                                                                             | `mol-geo/geometry/transform-data`    | `@molstar/core/math/linear-algebra/3d/transform-array` |

Collapsible controls are no longer exported from `@molstar/plugin-ui/base`. Gaussian-density CPU APIs remain in
`@molstar/core/math/geometry/gaussian-density`; their GPU functions/data types moved to graphics with the existing
signatures. The GPU helper behavior is carried over from v5; this audit does not redesign it. `fillIdentityTransform` is
no longer exported by the graphics transform-data module.

`Tokens` now lives in `@molstar/core/util/token-ranges`, with the same `data`, `count`, and `indices` fields. The IO
tokenizer retains its type re-export, so the mapped tokenizer import remains valid. `openRead` now lives in
`@molstar/common-server/open-read`; the volume server's mapped common/file module retains its re-export and behavior.

## WebGL context capabilities and limits

`WebGLContext` implements the backend-neutral `GraphicsContext` (`@molstar/graphics/gl/context`). Feature checks use
`webgl.capabilities` (`GraphicsCapabilities`, `@molstar/graphics/gl/capabilities`) instead of `webgl.extensions` or
`webgl.isWebGL2`, e.g. `capabilities.fragDepth`, `capabilities.floatRenderTarget` (was `extensions.colorBufferFloat`),
`capabilities.multipleRenderTargets` (was `extensions.drawBuffers`), `capabilities.texture3d` and
`capabilities.flatInterpolation` (were `isWebGL2`). `extensions` and `isWebGL2` remain for WebGL-specific code. New
capabilities `integerTextures` and `mixedAttachmentFormats` replace `isWebGL2` checks for `int` textures and for
framebuffers with attachments of different formats.

Depth texture types and depth renderbuffer formats are precision preferences: on WebGL1 `float`/`float-stencil` fall
back to `ushort`/`uint24-8` and `depth32f`/`depth32f-stencil8` to `depth16`/`depth-stencil` instead of throwing. `alpha`
render targets fall back to `rgba` on WebGL1, `Texture.mipmap()` is a no-op for non-power-of-two textures on WebGL1, and
`resources.getValueChannel(texture)` returns the channel holding the value of `alpha` format textures.

The limit getters `maxTextureSize`, `max3dTextureSize`, `maxRenderbufferSize`, `maxDrawBuffers` and
`maxTextureImageUnits` were removed from `WebGLContext`; use `webgl.resources.limits.*`.

Backend-neutral resource types moved out of `@molstar/graphics/gl/webgl/*`:

| Symbols                                                                                                                                                                                              | v6 import                            |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| `Texture`, `Textures`, `TextureValues`, `TextureKind`, `TextureType`, `TextureFormat`, `TextureFilter`, `TextureAttachment`, `createNullTexture`, `isNullTexture`, `decodeImage`, `CubeFaces`        | `@molstar/graphics/gl/texture`       |
| `Framebuffer`                                                                                                                                                                                        | `@molstar/graphics/gl/framebuffer`   |
| `Renderbuffer`, `RenderbufferFormat`, `RenderbufferAttachment`                                                                                                                                       | `@molstar/graphics/gl/renderbuffer`  |
| `RenderTarget`                                                                                                                                                                                       | `@molstar/graphics/gl/render-target` |
| `ArrayType`, `AttributeItemSize`, `AttributeKind`, `AttributeValues`, `ElementsType`, `ElementsKind`, `UsageHint`, `DataType`                                                                        | `@molstar/graphics/gl/buffer`        |
| `TimerResult`, `printTimerResults`, `GraphicsTimer` (was `WebGLTimer`)                                                                                                                               | `@molstar/graphics/gl/timer`         |
| `GraphicsStats` (was `WebGLStats`)                                                                                                                                                                   | `@molstar/graphics/gl/stats`         |
| `DrawMode`, `Transparency`, `GraphicsRenderVariant(s)`, `ComputeRenderVariant(s)`, `ProgramVariant`, `RenderItem`, `GraphicsRenderItem`, `ComputeRenderItem`, `MultiDrawBaseData`, `GraphicsProgram` | `@molstar/graphics/gl/render-item`   |
| `UniformBuffer`, `UniformKind`, `UniformsList`, `UniformValues`, ...                                                                                                                                 | `@molstar/graphics/gl/uniform`       |

The neutral `Texture` has no WebGL fields (`target`, `format`, `internalFormat`, `type`, `filter`) and no
`bind`/`unbind`; WebGL code uses `GLTexture` from `@molstar/graphics/gl/webgl/texture`. `createNullTexture` takes no
arguments.

`RepresentationContext.webgl` and `VisualContext.webgl` are renamed to `graphics` and typed `GraphicsContext`, e.g.
`provider.factory({ graphics: plugin.canvas3d?.graphics, ...themes }, getParams)`. Visual callbacks (`mustRecreate`,
`processValues`, `setOverpaint`/`setTransparency`/`setEmissive`/`setSubstance`/`setWiggle`), visual constructors
(`visualCtor(materialId, data, props, graphics)`) and `applyMesh*Smoothing` take a `GraphicsContext`; use
`isWebGLContext` from `@molstar/graphics/gl/webgl/context` for WebGL-only code paths.

The visual constructors of `UnitsRepresentation`, `ComplexRepresentation`, `VolumeRepresentation` and
`ParticleRepresentation` are a `VisualFactory` (`@molstar/graphics/repr/visual-factory`): either
`(materialId) => visual`, creating the same visual whatever the data, props and context, or a factory made by
`selectVisual` (choosing by the data, props and context) or `chooseVisual` (by the props and context), which names the
visuals it may create for signatures. A selecting factory gets the data its visuals' `mustRecreate` gets: the symmetry
group (`StructureGroup`) of a units representation, `{ volume, key }` (`VolumeKey`) and `{ particles }`
(`ParticleKey`) in place of the separate arguments. `withSignatureVisuals` is removed; impostor and mesh visuals no
longer recreate themselves (`mustRecreate`), their factory does when its choice changes.

Render items are created through the context: `createGraphicsRenderItem(ctx, ...)` and
`createComputeRenderItem(ctx, ...)` from `@molstar/graphics/gl/webgl/render-item` become
`ctx.createGraphicsRenderItem(...)` and `ctx.createComputeRenderItem(...)`. `RenderItem.getProgram` and
`Renderable.getProgram` return the neutral `GraphicsProgram` (`id`, `variant`, `isReady`, `link`, `finalize`); the WebGL
`Program` extends it. `namedComputeRenderables` is part of `GraphicsContext`. The renderable constructors
(`MeshRenderable`, ...), `createRenderable`, `Scene.create` and the camera, handle, pointer, shader and XR helpers take
a `GraphicsContext`.

Extensions take a `GraphicsContext` named `graphics` instead of a `WebGLContext` named `webgl`: debug-helpers
constructors (and `DebugRegistry.ctx`), sb-ncbr `TunnelShapeParams.graphics` (was `webgl`), alpha-orbitals
`createSphericalCollocationGrid`/`createSphericalCollocationDensityGrid` and geo-export `RenderObjectExporter.add` and
`AddMeshInput.graphics` (was `webgl`). `canComputeGrid3dOnGPU(graphics, computable)` accepts a `GraphicsContext` and the
computable of `createGrid3dComputeRenderable` (true on WebGPU for specs with WGSL code, `wgsl`), and computables take a
`GraphicsContext`; geo-export texture readback throws without a WebGL context.

Texture readback is asynchronous and backend-neutral: `GraphicsContext.readTextureData(texture)` reads a whole 2d color
texture (data order, four values per texel; ubyte, float and half float). `readTexture` and `readSingleChannelTexture`
(`@molstar/graphics/gl/compute/util`, WebGL, synchronous) were removed; the geo-export `MeshExporter.getInterpolated*`
helpers return promises. `TextureMesh.meta.webgl` was removed: the position iterator of texture meshes uses
`meta.vertexData`, read back with `TextureMesh.readVertexData` (visuals do it before colors per vertex).

The fragment-shader GPU compute is WebGL backend code: `gl/compute/histogram-pyramid/{reduction,sum}` and
`gl/compute/marching-cubes/{active-voxels,isosurface,tables}` moved to `@molstar/graphics/gl/webgl/compute/…`, and
`getGrid3dComputeFrag`/`getGrid3dComputeSchema` to `@molstar/graphics/gl/webgl/compute/grid3d`. Neutral code reaches
them through `getAtlasIsosurfaceExtraction(graphics)` (`@molstar/graphics/gl/compute/marching-cubes/atlas`) and
`createGrid3dComputeRenderable`, which use what the backend of the context registered. `canComputeGrid3dOnGPU(graphics)`
without a computable reports whether fragment shaders are available (a backend registered grid3d). `KernelCode` takes
optional `sources` (`wgsl`) for kernels generated at runtime.

Picking and Hi-Z read their textures through `GraphicsContext.createTextureReader({ format })`
(`@molstar/graphics/gl/texture-reader`, `rgba8` or `r32f`, a new member every context implements). `PickBuffers`
(`canvas3d/passes/pick`) is a class over a `TextureReader` (with an `asynchronous` getter) instead of an interface with
WebGL and WebGPU implementations; `checkAsyncPickingSupport` takes the `PickBuffers` instead of the context;
`PickPass.bindObject`/`bindInstance`/`bindGroup`/`bindDepth` were removed (use `PickPass.getTextures()`). `HiZPass`
works on WebGPU too; `HiZSchema` has a `uTexSize` uniform (the size of the level rendered).

`GraphicsCapabilities.floatBlend` (new) is rendering into float32 color targets with blending (WebGL: `EXT_float_blend`,
WebGPU: `float32-blendable`); `floatRenderTarget` is rendering without blending, always on WebGPU (it meant
`float32-blendable` there). Wboit, dpoit and illumination need `floatBlend`, as do the float32 fallbacks of blended
targets without half-float ones (multi-sample, marking, the Gaussian density). `requestWebGPUDevice` and
`createWebGPUBackend` take `features` (`'optional'`, the default, or `'none'`).

`GaussianDensityGPU` (`@molstar/graphics/geo/gaussian-density/gpu`, the WebGL density read back into a CPU field) was
removed; it had no callers. Use `computeGaussianDensity` (`@molstar/core/math/geometry/gaussian-density`, CPU) for a
field, or `computeGaussianDensityTexture` (`@molstar/graphics/geo/gaussian-density`) for a GPU texture.

Fixed-function state goes through the context instead of `webgl.state` and `gl` calls: `ctx.setRenderState(state)` with
a neutral `RenderState` (`@molstar/graphics/gl/render-state`: blend, color write mask, depth, stencil, cull, front face;
unset fields reset to the defaults, presets `RenderState.Default`/`Additive`/`Premultiplied`/`Alpha`/`Constant`),
`ctx.setViewport`, `ctx.setScissor`/`ctx.disableScissor`, `ctx.setBlendConstant`, `ctx.clear({ color, depth, stencil })`
(within the scissor rectangle, ignoring write masks, integer attachments included) and `ctx.flush()`.
`WebGLContext.clear(red, green, blue, alpha)` (bind and clear the drawing buffer) was removed. Draw buffers are set with
`framebuffer.setDrawBuffers(count)` (kept on context restore) instead of the `drawBuffers` extension, and
`texture.copyFromFramebuffer(framebuffer, x, y, width, height, dstX, dstY)` replaces `GLTexture.bind` +
`copyTexSubImage2D`. Most passes (`DrawPass`, postprocessing, antialiasing, bloom, DOF, marking, multi-sample,
illumination, tracing, WBOIT, DPOIT, background) take a `GraphicsContext`. `ShadowPass.webgl` is renamed to `graphics`.

Drawing happens in render passes: `ctx.beginRenderPass(target, { clear })` (a render target, framebuffer or
`'drawing-buffer'`) and `ctx.endRenderPass()` replace binding targets before drawing. `bindDrawingBuffer` moved from
`GraphicsContext` to `WebGLContext` (readback only). Copies, readback, attachment changes and `flush` end the open pass;
in debug mode draws and clears outside of a pass throw.

`Renderer.create` takes a `GraphicsContext`. `GraphicsRenderable.render(variant, globals)` and
`RenderItem.render(variant, globals?, mdbDataList?)` take `GlobalBindings`
(`ctx.createGlobalBindings(schema, uniforms, textures)`) instead of the number of shared textures; global uniforms and
shared textures are applied by the render items.

Shader sources compose their chunks in TS instead of `#include`: import the chunk modules from
`@molstar/graphics/gl/webgl/shader/chunks/*` and interpolate them (`` `${common}` ``); registering a program throws on
`#include`. Conditionals are resolved by a strict TS preprocessor (`@molstar/graphics/gl/shader-preprocessor`) before
compilation: `#if`/`#elif` support `defined`, `!`, `&&`, `||`, parentheses and integer comparisons of integer literals
and defines with integer values, plus `#ifdef`/`#ifndef`/`#else`/`#endif`, object- and function-like `#define` and
`#undef`; anything else (arithmetic in expressions, undefined names in evaluated expressions, other directives such as
`#extension` or `#version` in sources) throws a `ShaderPreprocessorError`.

Points are drawn as quads (two triangles per point) with the point data read from a texture: `PointsValues` have
`tPositionGroup`, `uTexDim`, `centerBuffer` and `groupBuffer` instead of the `aPosition` and `aGroup` attributes,
`drawCount` is `pointCount * 6` (`uVertexCount` stays the point count), and `Points` has `shaderData`. The point size is
no longer limited by the GPU point size range. Texture formats follow WebGPU: `rgb` is removed from `TextureFormat`, and
`alpha` is renamed `r` (also for `RenderTargetFormat` and `pixelPack`). `r` textures are `R8`/`R16F`/`R32F`/`R32I` on
WebGL2 and `LUMINANCE` on WebGL1; shaders read the value from `.r` (previously `.a` for 8-bit data). Color, palette and
size textures (`tColor`, `tPalette`, `tSize`) have 4 channels per item (colors with an opaque alpha, sizes packed into
`rgb` as before), so code reading their arrays uses a stride of 4. `resources.getValueChannel` is removed:
`extractIsosurface`, `createIsosurfaceBuffers` and `calcActiveVoxels` take the `valueChannel` (`red` or `alpha`) of the
volume texture after it. `readAlphaTexture` is renamed `readSingleChannelTexture` and returns the values in `.r`.

The GLSL define for the `aVertex` vertex id emulation is `emulateVertexId` (set where `gl_VertexID` is unusable),
replacing `noNonInstancedActiveAttribs`; `usedIf` conditions referring to it change accordingly.

The global texture `tDepth` is only the opaque scene depth; depth with alpha packed into `rgba` (front transparent depth
for the emissive variant, marking depth) is bound as `tDepthAlpha`, and custom shaders reading packed depth with alpha
from `tDepth` switch to `tDepthAlpha`.

`Canvas3DContext` has a `graphics: GraphicsContext` field (custom implementations of the interface add it) and an async
`Canvas3DContext.create(canvas, assetManager, backends, attribs, props)` that takes the backends to try in order
(`WebGL2Backend`, `WebGLBackend` from `@molstar/graphics/gl/webgl/backend`, `WebGPUBackend` from
`@molstar/graphics/gl/webgpu/backend`). `Canvas3DContext.fromCanvas` moved to
`createWebGLCanvas3DContext(canvas, assetManager, attribs, props)` in `@molstar/graphics/canvas3d/webgl`, so that
importing `Canvas3D` doesn't bundle the WebGL backend;
`Canvas3DContext.fromGraphics(canvas, graphics, assetManager, attribs, props)` wraps a context created otherwise. The
plugin's default backends (without `PluginConfig.General.GraphicsBackends`) come from
`@molstar/graphics/gl/backends/default` (the `all` set: WebGPU where capable, else WebGL2, then WebGL); builds with one backend alias that module
to `@molstar/graphics/gl/backends/webgl` or `@molstar/graphics/gl/backends/webgpu`.

`Canvas3D` and `Canvas3DContext` have a `graphics: GraphicsContext` field; `Canvas3D.webgl` and `Canvas3DContext.webgl`
were removed (narrow `graphics` with `isWebGLContext` from `@molstar/graphics/gl/webgl/context` where WebGL is
required). The WebGL context handles the loss and restore events of its canvas itself; `GraphicsContext.recreate` is
implemented by every backend (a new context replacing a lost one; on WebGL once the browser restored the context),
`simulateContextLoss()` loses the context for debugging, `release()` loses it on purpose to free its GPU resources
(WebGL `WEBGL_lose_context`, WebGPU the device destroyed; custom implementations add it), and `destroy` takes `GraphicsDestroyOptions` (`doNotForceWebGLContextLoss`). The passes and
pick helpers take a `GraphicsContext`; `PickBuffers` (a class, see the picking paragraph above) is created with
`createPickBuffers(graphics, pickPass, maxAsyncReadLag?)`. `GraphicsContext` has `readRenderTarget` (async readback of
render targets, rows top to bottom) and `waitForGpuCommandsComplete`.

## Image textures and pending work

`loadImageTexture(src, cell, texture)` is removed; use `resources.loadImage(texture, source, cell?, label?)` of the
graphics context, which returns a promise and is tracked in `resources.pending` (`isReady`, `whenReady`, `settled`).
`ImagePass.render` waits for pending work, `Canvas3D` redraws when it settles.

## Program descriptors and GLSL sources

GLSL sources moved from `@molstar/graphics/gl/shader/*` to `@molstar/graphics/gl/webgl/shader/*`. `ShaderCode` is a
descriptor without sources: `ShaderCode(name, extensions?, ignoreDefine?)` (no `vert`, `frag` and `outTypes`); each
backend looks the sources up by the program name, which must be unique. Custom GLSL programs register their sources with
`registerGlslProgram(name, vert, frag, outTypes?)` from `@molstar/graphics/gl/webgl/programs`. The GLSL assembly
(`addShaderDefines`, `getContextDefineLookup`, `stripUniformDeclarations`, `stripSamplerDeclarations`,
`stripAttributeDeclarations`, `ShaderUniformDeclarations`) moved to `@molstar/graphics/gl/webgl/shader-code`;
`addShaderDefines` returns the GLSL sources (`{ vert, frag }`). The program name `accumulate` is now
`color-smoothing-accumulate` (texture mesh color smoothing); tracing's is gone (`trace` blends into its average).
Marking's `edge` program is gone with `EdgeSchema` and `EdgeShaderCode` (`canvas3d/passes/marking`): the `overlay`
program computes the edges from the mask.

## Color type define

The `dColorType` define of renderables has the values `uniform`, `texture`, `grid` and `direct`: the color types read
from `tColor` by index (instance, group, groupInstance, vertex, vertexInstance) share `texture`, the grid types (volume,
volumeInstance) share `grid`, one program each, and the object uniform `uColorIndexType` selects the index or the grid
position (`ColorTypeIndex` in `@molstar/graphics/geo/geometry/color-data`). The color type itself is the value `colorType` of `ColorData`. Code that
set `dColorType` directly (custom geometry or renderables) sets the color type with `setColorType(values, type)`, which
keeps the three in step; code that read the type from `dColorType` reads `colorType`. Custom shaders test
`dColorType_texture` or `dColorType_grid` and branch on `uColorIndexType`.

## Size type define

As the color type: `dSizeType` has the values `uniform` and `texture`; the size types read from `tSize` by index
(instance, group, groupInstance, vertex, vertexInstance) share `texture`, and the object uniform `uSizeIndexType`
selects the index (`SizeTypeIndex` in `@molstar/graphics/geo/geometry/size-data`). The size type itself is the value
`sizeType` of `SizeData`, set with `setSizeType(values, type)`; code that read the type from `dSizeType` reads
`sizeType`. Custom shaders test `dSizeType_texture` and branch on `uSizeIndexType`.

## Clip object uniforms

The six clip object uniforms of renderables (`uClipObjectType`, `uClipObjectInvert`, `uClipObjectPosition`,
`uClipObjectRotation`, `uClipObjectScale`, `uClipObjectTransform`) are two: `uClipObjectFrame` (per object the transform
of a point into its frame) and `uClipObjectShape` (half the scale, the type in `w`, negative where inverted), derived by
`Clip.getClip` (`Clip.Objects.frame`, `shape`). The objects as given are the value `clipObjects`; code that read the
uniforms reads it. Custom shaders test shapes in the frame (`localSD`, `clipObjectSD` in GLSL, `localSD` and
`clipLocal` in WGSL); `getSignedDistance` with position, rotation, scale and transform remains for image trims. Whether
there are clip objects is the define `dClipObjects`, their count also the material uniform `uClipObjectCount` (custom
GLSL tests `defined(dClipObjects)` and loops to `uClipObjectCount`; WGSL to the override `oClipObjectCount`, as
`dClipObjectCount` sets it).

## Draw lists on by default

`RendererProps.drawLists` defaults to true: on WebGPU the passes replay recordings (render bundles), objects changing
every frame are drawn directly; WebGL ignores it. Set it to false for direct drawing (e.g. to compare).

## Renderable draw keys

`Renderable` has `programVersion` (incremented by `setTransparency`) and `drawVersion` (advanced when the instance and
level of detail ranges it draws change, see `createDrawVersion`), which the renderer's draw lists compare; custom
implementations add them. `GraphicsContext.createDrawList` is optional (WebGPU render bundles).

## Instance granularity

The instance granularity is the object uniform `uInstanceGranularity` (was the value `instanceGranularity`): the marker
and layer data are read per instance or per group instance by it at runtime, one program for both. The defines
`dMarkerType`, `dClippingType` and `dWiggleType` are removed; `dOverpaintType`, `dTransparencyType`, `dEmissiveType`
and `dSubstanceType` have the values `texture` (instance and groupInstance), `vertexInstance` and `volumeInstance`
(`getLayerTypeDefine` in `@molstar/graphics/geo/geometry/layer-type`). The types themselves are the values
`markerType`, `overpaintType`, `transparencyType`, `emissiveType`, `substanceType`, `clippingType` and `wiggleType`;
code that set a layer's type define sets the type with `setOverpaintType(values, type)` (and the other layers' setters),
code that read the type from the define reads the value. Custom shaders test `d…Type_texture` and select the index
with `uInstanceGranularity`.

## Data arrays (`DataSpec`)

The data read by element index (`tColor`, `tSize`, `tMarker`, `tOverpaint`, `tTransparency`, `tEmissive`,
`tSubstance`, `tClipping`, `tWiggle`, and `tPositionGroup` of spheres and points) are `DataSpec`s (`rgba8`, `r8` or
`rgba32f`, `@molstar/graphics/gl/renderable/schema`), not textures: their values are `DataArray`s, `{ array, count }`
(`createDataArray`, `@molstar/graphics/gl/renderable/util`), a flat array of `count` elements (it may be longer) without
a 2D layout, and the size uniforms `uColorTexDim`, `uSizeTexDim`, `uMarkerTexDim`, `uOverpaintTexDim`,
`uTransparencyTexDim`, `uEmissiveTexDim`, `uSubstanceTexDim`, `uClippingTexDim`, `uWiggleTexDim` and `uTexDim` are
removed from the values (the WebGL backend lays the arrays out and owns them). The 2D layouts of the grid layers
(`tColorGrid`, `tOverpaintGrid`, `tTransparencyGrid`, `tEmissiveGrid`, `tSubstanceGrid`) are the new
`u…GridTexDim` uniforms (were the layers' `u…TexDim`). `createTextureColor` and `createTextureSize` take a `DataArray`,
`getTrilinearlyInterpolated` (mesh color smoothing) returns one, `CenterGroupData` has no `texDim`. Shaders read the
arrays with generated accessors, `read_tColor(i)` (GLSL and WGSL), instead of `readFromTexture(tColor, i,
uColorTexDim)`. `Texture.loadRows` and `TextureChanges.getRows` are removed: data textures load with
`Texture.loadData`; a consumer loading a layer into a texture of its own lays it out with `getDataTextureImage`.
`getValueVersions` and `Versions` (`@molstar/graphics/gl/renderable/schema`) are removed: the render items keep their
versions in arrays per list.

## Context loss and headless WebGPU

`GraphicsContext` has a `contextLost` subject (emitted by `setContextLost`), `recreate()` for a new context replacing a
lost one (WebGPU: a new device) and `release()`; custom implementations of the interface add them.
`Canvas3DContext.contextLost` forwards the context's losses, also on WebGL; a lost context is replaced and everything
built again on it (`contextRestored`). `Canvas3DContext.release()` frees the GPU resources of a canvas on purpose (e.g.
out of view), reported as a loss but not restored until `restore()`, on both backends (`isReleased`). `requestWebGPUDevice(powerPreference, gpu?)` takes the `GPU` to use (e.g. Dawn's in Node.js) and
`WebGPUContextOptions.requestDevice` replaces lost devices. Headless: `ExternalModules.gl` is optional, `webgpu` (e.g.
the `webgpu` package) enables the WebGPU backend with `HeadlessScreenshotHelperOptions.backend` and the async
`HeadlessScreenshotHelper.create` / `HeadlessPluginContext.create`; images are read with `readRenderTarget` on both
backends.

## Graphics backend option of the Viewer and the Mesoscale Explorer

The `preferWebgl1` option and the `prefer-webgl1` URL parameter of the Viewer and the Mesoscale Explorer are removed:
`graphicsBackend: 'webgl1'` (`graphics-backend=webgl1`) asks for WebGL1. The default is `graphicsBackend: 'auto'`
(WebGPU where it has the features of WebGL2, else WebGL2, then WebGL1). The Prism viewer's `renderer.preferWebgl1`
option and `prefer-webgl1` URL parameter are removed as well (`graphics-backend=webgl1`).

WebGL1 is asked for by listing only its backend, nothing else: `PluginFeatureDetection.preferWebGl1` (iOS and Safari
15.1–15.3 no longer default to WebGL1), `PluginConfig.General.PreferWebGl1` (use
`PluginConfig.General.GraphicsBackends: [WebGLBackend]`), `Canvas3DContext.Attribs.preferWebGl1`, the `preferWebGl1`
attribute of `getGLContext` and the `preferWebGl1` option of the backend sets' `getDefaultGraphicsBackends` are
removed.

## Picking results

`Canvas3D.identify(target)` returns `PickResult`, i.e. `PickData | 'pending' | undefined` (was `PickData | undefined`):
`'pending'` where the backend reads pick data asynchronously (WebGPU) and the read for the current scene and camera
hasn't returned yet, so that it is no longer mistaken for "nothing at the target". Callers that use the result directly
check for `'pending'`, or await `canvas3d.asyncIdentify(target)?.wait()`, which resolves with the pick data on every
backend (`AsyncPickData` has a `wait()` next to `tryGet()`). On WebGL, which reads synchronously, `identify` never
returns `'pending'`.

## GPU timer groups

`GraphicsTimer.mark`/`markEnd` are replaced by `group`/`groupEnd` (no aliases). Only render passes are timed: a group
sums the passes that begin inside it and lists them as leaves (`pass <label>`), named by the new `label` option of
`beginRenderPass` (derived from the target otherwise). Groups must begin and end outside of passes (checked in debug
mode); a group without passes reports its GPU time as `NaN` (printed as `n/a`). The `Renderer.*` timings and other
timings inside passes were removed, and timings of a single pass became the label of that pass.

## Reversed depth

Depth is reversed on both backends: the near plane at depth 1, the far plane (and the background, the cleared depth) at
0, so nearer fragments have greater depth. Camera projections stay GL projections; the backends convert the projection
uniforms, which are tagged with `UniformSpec.depthRange` (`'projection'` or `'inverse'`, `@molstar/graphics/gl/depth-range`):
custom renderable schemas holding a projection or its inverse add the tag. Custom shaders follow the new convention:
background tests compare with 0 (was 1), depth compares and min/max of depths flip, packed depth of the far plane is
zeros (was ones), and `perspectiveDepthToViewZ`/`orthographicDepthToViewZ` (GLSL and WGSL chunks) take reversed depth.
Depth render states of custom passes use `greater`/`greater-equal` and clear depth to 0. `PickBuffers.getDepth` still
returns GL window depth (near 0, far 1), for `cameraUnproject`. With `EXT_clip_control` WebGL clips depth to `[0, 1]`
(GLSL defines `clipZeroToOne`): custom GLSL converts between clip z / w and window depth with `clipZToDepth` and
`depthToClipZ` (common chunk) instead of `* 0.5 + 0.5` and `* 2.0 - 1.0`, and vertex shaders writing clip z keep it
inside the context's range.

## Top-left window origin

Window positions and rectangles have a top-left origin (y down) everywhere outside the WebGL backend, as in the
browser (CSS, pointer events, `ImageData`):
- `Viewport` and `Camera.viewport`; `cameraProject`, `cameraUnproject`, `Camera.project`, `Camera.unproject` and
  `getRay` take and return window coordinates with y down (consumers that converted with `height - y` drop it).
- `GraphicsContext.setViewport`, `setScissor`, `readRenderTarget` and `readTextureAsync` (rows top to bottom, as
  before), texture readers (rows now top to bottom, were bottom to top) and `Texture.copyFromFramebuffer`.
- The `static-frame` and `relative-frame` viewport params of `Canvas3D` count `y` from the top.
- Uniforms holding a window rectangle or uv offset are tagged with `UniformSpec.windowRect` (`bounds`, `rect`, `pixels`,
  `offset`, `@molstar/graphics/gl/window-rect`) and flipped by WebGL where it binds them; custom passes tag theirs.
  GLSL keeps GL's bottom-left `gl_FragCoord`; WGSL works top-left (`fragmentUv`, `read_<name>` without a flip,
  `screenSpaceToViewSpace` with y down). The Hi-Z uniform `uOffset` became `uRect` (corner and extent).
- `GraphicsXR` has `isPresenting()` and `getView(view, out)` (the view's viewport, top-left, and projection with GL's
  clip depth); `StereoCamera.update` takes `{ pose, graphics }` (the context's `GraphicsXR`) instead of
  `{ pose, layer }` (an `XRWebGLLayer`); custom implementations of the interface add both.
- Data rendered into tiles of a texture by rectangles (Gaussian density, color smoothing atlases) is placed by texel
  rows from the bottom, converted explicitly; the WebGL backend's own compute converts its texel rectangles internally.

## Lazy texture storage

`Texture.define(width, height)` records the size; the GPU storage is allocated when the texture is first written (a
render pass on a framebuffer it is attached to, a copy, an upload) or read back, not at `define`. Sampled before, a
texture reads zero (WebGPU binds a placeholder, WebGL keeps a 1×1 image). Renderbuffers are allocated on first use in
the same way. `getByteCount` of textures, renderbuffers, render targets and passes counts allocated storage only, so it
varies with what frames draw.

## Pass targets

The passes that draw into targets take a `PassTargets` (the registry their targets are created through, see
`canvas3d/passes/pass-targets.ts`) instead of the `GraphicsContext`: `DrawPass(targets, assetManager, width, height,
transparency)`, `WboitPass`, `DpoitPass`, `MarkingPass`, `BloomPass`, `DofPass`, `OutlinePass`, `ShadowPass`,
`SsaoPass`, `SmaaPass`, `AntialiasingPass`; `Passes.targets` and `DrawPass.targets` hold it. The storage of targets
the current props don't draw into is released (see Lazy texture storage), also of an `ImagePass` a second after its
last render. `Texture` and `Renderbuffer` have `releaseStorage`, `Framebuffer` has `getAttachments`; custom
implementations add them.

## SSAO at half resolution

`SsaoParams.resolutionScale` defaults to 0.5 (was 1): the occlusion is computed at half the resolution of the drawing
buffer and upsampled with the depth (postprocessing, illumination compose and the marking overlay). The scale is now at
most that of CSS pixels (`min(resolutionScale, 1 / pixelRatio)`, was `resolutionScale / pixelRatio`), so that the
default stays at half resolution on high pixel ratios. Set 1 for the previous full-resolution occlusion.

## External color themes

`ColorTheme.createRegistry()` in graphics no longer includes `external-structure` or `external-volume`, and the special
registration helper in `PluginContext` is gone. The implementations live at `@molstar/plugin/themes/external-structure`
and `@molstar/plugin/themes/external-volume` because they depend on plugin state objects and selection queries. Plugin
composition registers them again through the `ExternalColorThemes` entry (`@molstar/plugin/themes/external`), which
`DefaultRegistry` includes; specs with their own registry add the entry to use them.

## Plugin composition step 1

Step 1 of [plugin-composition.md](plugin-composition.md) split facades and large modules so that applications compose
the plugin from catalogs instead of importing everything. Whole-file ownership stays in `migration-map.json`; where a v5
module was split, its entry maps to an array of the files that received its contents, or to `null` when nothing succeeds
it. [`migration-symbols.json`](migration-symbols.json) maps each exported symbol of an affected v5 module to its current
defining module (`null` when removed). `scripts/workspace/check.mjs` verifies that every target and symbol in both files
still exists.

### `StateTransforms` removed

`StateTransforms` and `@molstar/plugin/state/transforms` no longer exist. Import each transformer from its defining
module; the `StateTransforms.<Namespace>.<Name>` keys in `migration-symbols.json` list all 120 members.

| v5                                                         | v6 import                                                   |
| ---------------------------------------------------------- | ----------------------------------------------------------- |
| `StateTransforms.Data.Download`                            | `@molstar/plugin/state/transforms/data/fetch`               |
| `StateTransforms.Data.ParseCif`                            | `@molstar/plugin/state/formats/cif`                         |
| `StateTransforms.Model.TrajectoryFromMmCif`                | `@molstar/plugin/state/formats/trajectory/mmcif`            |
| `StateTransforms.Model.StructureFromModel`                 | `@molstar/plugin/state/transforms/structure/hierarchy`      |
| `StateTransforms.Model.StructureComponent`                 | `@molstar/plugin/state/transforms/structure/selection`      |
| `StateTransforms.Representation.StructureRepresentation3D` | `@molstar/plugin/state/transforms/structure/representation` |
| `StateTransforms.Volume.VolumeFromCcp4`                    | `@molstar/plugin/state/formats/volume/ccp4`                 |

The classic Viewer global `molstar.lib.plugin.StateTransforms` keeps its shape (the same seven namespaces and 120
members). It is now an app-level object in `@molstar/viewer/state-transforms` for classic-script pages that cannot
import leaf modules. The `StateTransforms` type alias (`typeof` the facade) is gone with the facade. Neither
`@molstar/viewer/lib` nor `@molstar/viewer/state-transforms` exports a `StateTransforms` type by name; where the type is
needed, use `typeof StateTransforms` imported from `@molstar/viewer/state-transforms`.

### Transform and format modules split

The seven transform modules and the six format modules are split by functionality:
`state/transforms/{data,misc, particles,shape,structure,volume}/...` (structure effects live in `structure/effects/`)
and `state/formats/{trajectory,coordinates,topology,volume,shape,particles}/<format>.ts`. Format-specific transformers
(`ParseCcp4`, `VolumeFromCcp4`, `TrajectoryFromPDB`, ...) live next to their provider in the format module. Each format
family also has `provider.ts`, `category.ts`, and `catalog.ts` modules for shared helpers, the category, and the
built-in provider list.

`@molstar/plugin/state/transforms/catalog` imports every built-in transform and format module for its registration side
effects and exports nothing. The default spec imports it, so a plugin built from `DefaultPluginSpec` registers all
built-in transformer ids. A custom spec must import the modules it uses, or the catalog, before it restores snapshots
(see _Snapshot loading_).

### Default specs

`DefaultPluginSpec` moved to `@molstar/plugin/default-spec` and `DefaultPluginUISpec` to
`@molstar/plugin-ui/default-spec`. `@molstar/plugin/spec` and `@molstar/plugin-ui/spec` keep only types and helpers
(`PluginSpec`, `PluginUISpec`). The classic Viewer globals `molstar.lib.plugin.DefaultPluginSpec` and
`DefaultPluginUISpec` are unchanged.

Both default-spec modules load `DefaultRegistry` and every built-in catalog. Apps that compose their own registry and
only want the default UI parts import `DefaultPluginUIComponents()` and `DefaultPluginUICustomParamEditors()` from
`@molstar/plugin-ui/default-ui` instead; `DefaultPluginUISpec()` is composed from them.

### Script languages

PyMOL, VMD, and Jmol are no longer enabled implicitly. Import `@molstar/query-language/transpilers/pymol`, `.../vmd`,
`.../jmol`, or `@molstar/query-language/transpilers/all` to enable a language; the default plugin spec and the Viewer
import `all`. Using a language that was not imported throws `Script language '<x>' is not available in this build`.
`mol-script` is always available. `Script.getAvailableLanguages()` returns `mol-script` plus the registered languages,
and the script-language select in the plugin UI lists only those. The `_transpiler` object exported by `transpilers/all`
is removed; after enabling the language, call `parse(language, text)` from `@molstar/query-language/transpile` (or
`Script.toExpression`), or import `transpiler` from `.../transpilers/<lang>/parser` directly.

### Data format providers

`DataFormatProvider` requires a `name` and gains an `Id extends string` type parameter so the name keeps its literal
type; the `DataFormatProvider(...)` helper takes a `const` type parameter for the same reason. Custom providers must add
`name`. `BuiltInTrajectoryFormats`, `BuiltInVolumeFormats`, `BuiltInShapeFormats`, `BuiltInTopologyFormats`,
`BuiltInCoordinatesFormats`, and `BuiltInParticlesFormats` are now `as const` provider arrays instead of
`[name, provider]` tuples; read `provider.name` instead of the first tuple element. `PluginSpec.customFormats` and
`DataFormatRegistry.add(name, provider)` keep their `[name, provider]` shape in this step.

The format name types (`BuiltInTrajectoryFormat`, `BuiltInVolumeFormat`, ...) are derived from the catalogs and live in
`@molstar/plugin/state/formats/<family>/catalog`. `BuiltInVolumeFormat` and `BuiltInShapeFormat` are new; the misspelled
`BuildInVolumeFormat` and `BuildInShapeFormat` remain as deprecated aliases. `MmcifProvider.parse` and
`CifCoreProvider.parse` take an optional `params` argument. The volume formats' parameter type, formerly an unexported
`Params`, is exported as `VolumeFormatParams` from `@molstar/plugin/state/formats/volume/provider`.

### Graphics catalogs

The `BuiltIn` values were removed from the `ColorTheme` and `SizeTheme` namespaces and from
`StructureRepresentationRegistry`, `VolumeRepresentationRegistry`, and `ParticleRepresentationRegistry`. The `BuiltIn`
types, `BuiltInParams`, and `createRegistry()` keep their names; `ColorTheme.createRegistry()` and
`SizeTheme.createRegistry()` return empty registries
([Registries start empty](#plugin-composition-step-3-registries-start-empty)). Read the values from the catalogs:

| v5 value                                  | v6 catalog export                                                                 |
| ----------------------------------------- | --------------------------------------------------------------------------------- |
| `ColorTheme.BuiltIn`                      | `BuiltInColorThemes` from `@molstar/graphics/theme/color/catalog`                 |
| `SizeTheme.BuiltIn`                       | `BuiltInSizeThemes` from `@molstar/graphics/theme/size/catalog`                   |
| `StructureRepresentationRegistry.BuiltIn` | `BuiltInStructureRepresentations` from `@molstar/graphics/repr/structure/catalog` |
| `VolumeRepresentationRegistry.BuiltIn`    | `BuiltInVolumeRepresentations` from `@molstar/graphics/repr/volume/catalog`       |
| `ParticleRepresentationRegistry.BuiltIn`  | `BuiltInParticleRepresentations` from `@molstar/graphics/repr/particles/catalog`  |

Catalogs are built with the new `namedCatalog` helper (`@molstar/graphics/util/named-catalog`), which checks at compile
time that each key equals its provider's `name`; the runtime key/name check in the registry constructors is removed.

### Selection queries

`@molstar/plugin/state/helpers/structure-selection-query` is split into `@molstar/plugin/state/queries/structure/*`:
`catalog` (`StructureSelectionQueries`), `registry` (`StructureSelectionQueryRegistry`), `query`
(`StructureSelectionQuery`, `StructureSelectionCategory`), `dynamic` (`ResidueQuery`, `ElementSymbolQuery`,
`EntityDescriptionQuery`, and the `get*Queries` helpers), and the category modules `basic`, `bond`, `common`,
`manipulate`, `residue`, `structure`, and `type`. `applyBuiltInSelection` was unused and is removed; apply
`StructureSelectionFromExpression` with `StructureSelectionQueries[name].expression` instead.

### Structure presets

Representation presets moved from `state/builder/structure/representation-preset` to
`state/builder/structure/representation-presets/`, and hierarchy presets from `.../hierarchy-preset` to
`.../hierarchy-presets/`: one module per preset, plus `types` (`StructureRepresentationPresetProvider`,
`TrajectoryHierarchyPresetProvider`, `presetStaticComponent`, `presetSelectionComponent`) and `catalog`
(`PresetStructureRepresentations`, `PresetTrajectoryHierarchy`). Presets gain an optional `alias` (for example `auto`,
`default`, `all-models`) next to `id`; `PresetProvider` and both provider helpers take `Id` and `Alias` type parameters.
The catalogs export the name types `BuiltInStructureRepresentationPresetId`,
`BuiltInStructureRepresentationPresetAlias`, `BuiltInTrajectoryHierarchyPresetId`, and
`BuiltInTrajectoryHierarchyPresetAlias`. The `auto` preset is its own module, `representation-presets/auto`.

### Behaviors and markdown extensions

`BuiltInPluginBehaviors` moved to `@molstar/plugin/behavior/built-in`; `PluginBehaviors` stays in
`@molstar/plugin/behavior`. `@molstar/plugin/behavior` no longer re-exports `PluginBehavior`; import it from
`@molstar/plugin/behavior/behavior`. The built-in markdown extensions moved out of
`@molstar/plugin/state/manager/markdown-extensions` into `@molstar/plugin/state/markdown/*`: `BuiltInMarkdownExtension`
is exported from `.../markdown/catalog`, and the individual extensions are in the sibling modules (`audio`, `camera`,
`highlight`, `query`, `snapshots`). `MarkdownExtension`, `MarkdownExtensionManager`, and
`defaultParseMarkdownCommandArgs` stay in the manager module.

### Snapshot loading

Restoring a snapshot (`PluginState.setSnapshot`, and every entry of `PluginStateSnapshotManager.setStateSnapshot`) first
checks that each transformer id in the behavior tree, data tree, and transition frames is registered. A snapshot that
names unregistered transformers throws one error listing all missing ids
(`Snapshot uses transformers that are not available in this plugin: <ids>. Import the modules that define them.`) before
the plugin or the snapshot manager changes. `StateTransformer.has(id)` is a non-throwing registration check;
`StateTransformer.get` still throws for unknown ids.

## Plugin composition step 2

Step 2 of [plugin-composition.md](plugin-composition.md) replaces three spec fields with registry entries. Step 3
([Registries start empty](#plugin-composition-step-3-registries-start-empty)) then removes the constructor preloads
(formats, representations, themes, presets, selection queries, markdown extensions).

### `PluginSpec.actions`, `animations`, and `customFormats` removed

`PluginSpec` no longer has `actions`, `animations`, or `customFormats`. The `PluginContext` constructor throws
`PluginSpec.<key> was removed in 6.0; use registry entries (see the migration guide)` when any of the three keys has a
value other than `undefined`; a key present with `undefined` (for example `customFormats: o?.customFormats`) is ignored.
This is a check, not an alias. The replacements are entries of `PluginSpec.registry`:

| 5.x                                 | 6.0                                                                                                                                     |
| ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `actions: [PluginSpec.Action(a)]`   | `registry: [..., { actions: [a] }]`                                                                                                     |
| `animations: [...]`                 | `registry: [..., { animations: [...] }]`                                                                                                |
| `customFormats: [[name, provider]]` | `registry: [..., { formats: [DataFormatProvider.withName(provider, name)] }]`                                                           |
| `actions: defaultSpec.actions`      | `registry: DefaultRegistry` (or `defaultSpec.registry`), or the `DefaultActions` entry                                                  |
| `animations: [AnimateModelIndex]`   | replace the default entry: `registry: [...DefaultRegistry.filter((e) => e !== DefaultAnimations), { animations: [AnimateModelIndex] }]` |

`DefaultPluginSpec()` and `DefaultPluginUISpec()` return `{ registry: DefaultRegistry, behaviors, ... }` and no longer
carry `actions` or `animations`, so code that read `defaultSpec.actions`/`defaultSpec.animations` (or
`plugin.spec.actions`) must read the registry entries (`DefaultActions.actions`, `DefaultAnimations.animations`) or the
plugin (`plugin.state.data.actions`, `plugin.managers.animation.animations`). A literal spec that did not spread the
default spec gets no actions and no animations until it lists them (for example `registry: DefaultRegistry`); a spec
that listed `actions: defaultSpec.actions` and `animations: defaultSpec.animations` becomes
`registry: defaultSpec.registry`.

What changes now: actions and animations come only from `registry` (and behaviors); `PluginContext` no longer registers
`AnimateStateSnapshotTransition` when no animations are listed, so `plugin.managers.animation.current` is `undefined`
(it is typed `Current | undefined`) until an animation is registered. `AnimateStateSnapshotTransition` stays base
functionality: snapshot transitions and the UI snapshot controls play it directly and `play()` registers it on demand.
`DefaultAnimations` keeps it in its current position, so the default animation remains `AnimateModelIndex`.
`PluginAnimationManager.setSnapshot` skips `snapshot.current` and logs a warning when the snapshot names an unregistered
animation, instead of writing `paramValues` into the previously selected animation. `PluginContext.initCustomFormats`,
`initDataActions`, and `initAnimations` are removed (they were private).

What step 3 changes: see [Registries start empty](#plugin-composition-step-3-registries-start-empty); a literal spec
gets no formats, representations, themes, presets, selection queries, or markdown extensions either, and must list
`DefaultRegistry` (or its own entries).

### Viewer and mesoscale-explorer `customFormats`

`Viewer.create` options are unchanged, including `customFormats: [name, provider][]` with `DataFormatProvider.Unnamed`
providers. The Viewer builds
`registry: [...DefaultRegistry, { formats: customFormats.map(([n, p]) => withName(p, n)) }]` with the custom-formats
entry last. A custom name that equals a built-in format name (for example `['pdb', MyPdb]`) still overrides it:
`DefaultFormats` is replaced by a copy without that provider, so `dataFormats.get(name)` returns the custom provider and
no name is registered twice. The overriding provider is listed after the built-ins, so for that format `auto()`
tie-breaking by registration order sees it last, as in 5.x. mesoscale-explorer keeps its own `customFormats` option
shape and applies the same rule.

## Plugin composition step 3: registries start empty

`new PluginContext(spec)` (and `PluginUIContext`, `HeadlessPluginContext`) starts with empty registries. After `init()`
a plugin whose spec lists nothing has no representations, color or size themes (in the structure, volume, and particles
scopes), data formats, trajectory-hierarchy or representation presets, selection queries, markdown extensions,
drag-and-drop handlers, actions, or animations; it has only what `spec.registry` and the behaviors register. A literal
spec that relied on the implicit catalog adds `DefaultRegistry` (or its own entries):

```ts
// 5.x: a literal spec received every built-in provider
const spec = { behaviors: [...] };
// 6.0
import { DefaultRegistry } from '@molstar/plugin/default-registry';
const spec = { registry: DefaultRegistry, behaviors: [...] };
```

Script-tag code using the classic Viewer global does the same with `molstar.lib.plugin.DefaultRegistry`, which the
Viewer now exposes next to `DefaultPluginSpec` and `DefaultPluginUISpec`.

Specs built from `DefaultPluginSpec()` or `DefaultPluginUISpec()` list `DefaultRegistry` already and are unchanged.

- `ColorTheme.createRegistry()` and `SizeTheme.createRegistry()` return empty registries. Code that creates a
  `RepresentationContext` by hand (for example `colorThemeRegistry: ColorTheme.createRegistry()`) must add the themes it
  uses (`BuiltInColorThemes` from `@molstar/graphics/theme/color/catalog`, or single providers).
- `StructureRepresentationRegistry`, `VolumeRepresentationRegistry`, and `ParticleRepresentationRegistry` constructors
  add nothing; `DataFormatRegistry`, `TrajectoryHierarchyBuilder`, `StructureRepresentationBuilder`,
  `StructureSelectionQueryRegistry`, and `MarkdownExtensionManager` constructors add nothing either. The modules that
  define these classes no longer import any catalog.
- `PluginContext.dropOverriddenPreloadedFormats` (the transitional handling that let a `customFormats` entry replace a
  preloaded built-in format) is removed. The Viewer and mesoscale-explorer replace the matching default entry as
  described above, so a custom format with a built-in name still overrides it; a different provider under an existing
  format name in any other spec is an error, as for every registry.
- mesoscale-explorer no longer clears registries after `init()`; it lists spacefill with its themes, the `uniform` and
  `illustrative` color themes, the mmCIF format, its three animations, and its custom formats.
- The import-graph check requires `@molstar/plugin/context`, `@molstar/plugin/spec`, `@molstar/plugin-ui`, and the base
  UI entry points to reach no catalog module by value.

## Plugin composition step 3: UI

### `createPluginUI` requires a spec

`createPluginUI({ target, render, spec })` no longer falls back to `DefaultPluginUISpec()`; `spec` is a required option
and `@molstar/plugin-ui` (`index.ts`) no longer imports `@molstar/plugin-ui/default-spec`. Pass
`spec: DefaultPluginUISpec()` (import it from `@molstar/plugin-ui/default-spec`) for the previous behavior.

### Base UI renders minimal structure tools

The base `Plugin` layout (`ControlsWrapper`) no longer falls back to `DefaultStructureTools`. Without
`spec.components.structureTools` it renders `MinimalStructureTools`: structure source, measurements, components, and the
behavior-registered custom structure controls. The full tools component (also adding superposition, quick styles,
procedural animation, volume streaming, volume source, and particle source) moved from `@molstar/plugin-ui/controls` to
`@molstar/plugin-ui/default-structure-tools` and is set as `DefaultPluginUISpec().components.structureTools`.

A spec that spreads `DefaultPluginUISpec()` but replaces `components` wholesale must spread the default components to
keep the full tools, or it silently gets the minimal list:

```ts
const defaultSpec = DefaultPluginUISpec();
const spec = { ...defaultSpec, components: { ...defaultSpec.components, remoteState: 'none' } };
```

### Quick styles and volume controls

Quick Styles resolves its presets by id through the representation preset registry (`DefaultRepresentationPreset`
config, then `preset-structure-representation-auto` for Default; `-polymer-and-ligand`, `-illustrative`, and
`-molecular-surface` for Cartoon, Spacefill, and Surface) and hides the buttons whose preset is not registered.
Lazy-volume loading in the volume source controls uses the isosurface representation and uniform color theme only when
they are registered, and otherwise falls back to the registry default.

## Plugin composition step 3

Step 3 of [plugin-composition.md](plugin-composition.md) removes implicit defaults from the plugin context.

### View models and hooks require a spec

`PluginViewModel` and `PluginUIViewModel` (`@molstar/plugin-extension`) take `{ spec }` as a required constructor option
and no longer fall back to `DefaultPluginSpec()` / `DefaultPluginUISpec()`. `useCreatePluginViewModel` and
`useCreatePluginUIViewModel` require `options.spec` (a spec object) and no longer accept the
`spec: (defaultSpec) => spec` callback form. Callers import the default spec themselves. The classic Viewer global
(`molstar.PluginExtensions.plugin.models`) exposes app-level subclasses that keep the 5.x default when no spec is given
(`apps/viewer/src/view-models.ts`).

```ts
// 5.x
useCreatePluginUIViewModel({ spec: (spec) => ({ ...spec, behaviors: [...spec.behaviors, MyBehavior] }) });
// 6.0
const spec = DefaultPluginUISpec();
useCreatePluginUIViewModel({ spec: { ...spec, behaviors: [...spec.behaviors, MyBehavior] } });
```

### Drag-and-drop open-files fallback is a registry entry

`DragAndDropManager.handle()` no longer opens unrecognized files itself. The open-anything handler is the
`{ name: 'open-files', handle, fallback: true }` entry of `DefaultDragAndDrop` (part of `DefaultRegistry`), so it still
runs after the handlers and the built-in session handling (`.molx`/`.molj`). A plugin without that entry ignores drops
it does not recognize; add `DefaultDragAndDrop` (or your own fallback handler) to open arbitrary files.

### Snapshot name warnings

`PluginState.setSnapshot` calls the new `PluginState.reportUnregisteredNames(snapshot)` after the behavior tree is
applied and before the data tree. It logs a warning for each representation type, color theme, and size theme name of
the structure, volume, and particles representation transformers (data tree and transition frames) that is not
registered in the matching scope; the registry default replaces it when the data tree is normalized. It never throws.

## Plugin composition step 3: presets

### Mixed name and provider props

`createStructureRepresentationParams`, `createVolumeRepresentationParams`, and the structure representation builder
(`addRepresentation`, `buildRepresentation`) accept a name or a provider in each of `type`, `color`, and `size` and
resolve each field on its own, so `{ type: CartoonRepresentationProvider, color: 'element-symbol' }` works. The "any
string field selects the by-name path" dispatch is gone. The params types follow the field: a provider gives its params,
a built-in name gives `BuiltInParams`, and any other string gives `{}`. The prop types are now
`StructureRepresentationProps<R, C, S>` and `VolumeRepresentationProps<R, C, S>` with `R`, `C`, and `S` constrained to a
provider or a string (the new `*RepresentationRef`, `*ColorThemeRef`, and `*SizeThemeRef` types);
`*RepresentationBuiltInProps` is kept as the name-only alias. The builder methods are generic in `R`, `C`, and `S`
instead of in the whole props object, so the params of each field are checked against that field; a call that previously
compiled only because the props were checked against the union of every built-in's params can now fail to compile (for
example `colorParams: { palette }` on `illustrative`). Names that are not registered still warn and use the registry
default, and an empty registry still throws.

### Presets import what they run and export entries

Each representation preset passes provider objects (and the imported theme providers for its fixed color themes) to the
builder instead of names, and its module exports an entry next to the preset: `EmptyPresetEntry`, `AutoPresetEntry`,
`AtomicDetailPresetEntry`, `PolymerCartoonPresetEntry`, `PolymerAndLigandPresetEntry`, `ProteinAndNucleicPresetEntry`,
`CoarseSurfacePresetEntry`, `IllustrativePresetEntry`, `MolecularSurfacePresetEntry`, `AutoLodPresetEntry`, and
`MesoscalePresetEntry`. An entry lists its own preset first, followed by the representation entries (and theme
providers) it builds; `AutoPresetEntry` also includes the entries of the presets `auto` chooses between. The hierarchy
presets export `DefaultHierarchyPresetEntry`, `AllModelsHierarchyPresetEntry`, `UnitcellHierarchyPresetEntry`,
`SupercellHierarchyPresetEntry`, and `CrystalContactsHierarchyPresetEntry` with the color themes they apply (and, for
all-models, the default preset it delegates to). Hierarchy preset entries do not include a representation preset; that
is the configured one. `@molstar/plugin/registry/merge` exports `mergeRegistryEntries(...entries)` for building such
entries. `DefaultPresets` is built from these entries and lists the same 5 hierarchy and 11 representation presets in
the same order. The preset providers keep their names.

### Hierarchy preset `representationPreset`

The `representationPreset` param of every hierarchy preset is typed with the built-in id and alias unions (`import type`
from the catalog) and no longer defaults to `'auto'`: when it is absent or empty the preset applies
`PluginConfig.Structure.DefaultRepresentationPreset`, and the hierarchy presets no longer import the `auto` preset. The
configured preset is resolved through the registry, so a hierarchy preset applied to a plugin that has not registered it
fails with `Preset '<id>' is not registered in this plugin`. Calls that relied on the old `'auto'` param default and a
config that names another preset now get the configured preset.

### `presetSelectionComponent`

`presetSelectionComponent(plugin, structure, query, tag, params?)` takes a `StructureSelectionQuery` object and a tag
instead of a `StructureSelectionQueries` key; the component key stays `selection-<tag>`. Import the query from its
module (for example `protein` and `nucleic` from `@molstar/plugin/state/queries/structure/type`). The module no longer
imports the query catalog. Replace `presetSelectionComponent(plugin, s, 'protein')` with
`presetSelectionComponent(plugin, s, protein, 'protein')`.

### Delegating presets

The Viewer's `ViewerAutoPreset` looks up the model-archive quality-assessment and SB-NCBR partial-charges presets by id
through `plugin.builders.structure.representation.resolveProvider` and skips them when they are not registered, instead
of importing them; it still falls back to the imported `AutoPreset`.

## Plugin composition step 3: formats and post-load presets

### Format entries

Every built-in data format module exports a registry entry next to its provider: the provider name without `Provider`
(`SdfProvider` and `Sdf`, `MmcifProvider` and `Mmcif`, `Ccp4Provider` and `Ccp4`, `RelionStarParticlesProvider` and
`RelionStarParticles`). An entry holds `formats: [provider]` and the actions `DefaultActions` lists for the format
(CCP4: `ParseCcp4`, `VolumeFromCcp4`; mmCIF: `ParseCif`, `TrajectoryFromMmCif`; PDB, PDBQT, and PQR:
`TrajectoryFromPDB`; SDF: none), so a slim plugin registers `registry: [Sdf, Ccp4, ...]`. Volume and particle format
entries also include the representations and themes their `visuals` step builds: `Isosurface` for CCP4, DSN6, DX, Cube,
density-server CIF, structure-factor CIF, and MTZ, `Segment` for segmentation CIF, and the particle spacefill
representation or the spacefill, fibers, and target representations with their color themes for the particle formats.
Shape, topology, coordinates, and trajectory formats carry no representations. The `visuals` steps pass provider objects
(`IsosurfaceRepresentationProvider`, `UniformColorThemeProvider`) instead of names, and
`VolumeRepresentation3DHelpers.getDefaultParams` and `getDefaultParamsStatic` accept a provider object or a name for the
representation and the color and size themes. `DefaultFormats` and `DefaultActions` are unchanged, and the built-in
catalogs (`BuiltIn*Formats`) still list the providers.

### Post-load presets come from config

Trajectory format `visuals` and the `DownloadStructure` action apply `PluginConfig.Structure.DefaultHierarchyPreset`
(default `preset-trajectory-default`) instead of `'default'`. `LoadTrajectory`, `AddTrajectory`,
`StructureHierarchyManager.updateStructure`, and the Cube format apply
`PluginConfig.Structure.DefaultRepresentationPreset` (default `preset-structure-representation-auto`) instead of
`'auto'`, and `DownloadStructure` takes both the default of its representation select and the empty-preset check from
config. A config value that is not registered fails the load with `Preset '<id>' is not registered in this plugin` (the
actions log it and revert, as for any other error in their transaction); there is no fallback to the default preset. A
plugin that does not register the default presets sets the two config items.

### `DownloadStructure` depends on the registered formats

`DownloadStructure` offers only the sources whose format is registered (`plugin.dataFormats.has`): PDB, PDB-IHM,
AlphaFold DB, and Model Archive need `mmcif`, SWISS-MODEL needs `pdb`, PubChem needs `mol`, and URL needs at least one
trajectory format. The URL format list is the registered formats of the trajectory category in registration order
(formats that extensions register in that category, such as `g3d` in the Viewer, now appear in it;
`BuiltInTrajectoryFormat` no longer types the `format` param, which is a `string`), and its default is `mmcif` when
registered. The "Load all entries into a single trajectory" option is hidden unless `mmcif` is registered. With no
usable source the action's source select has a single empty option.

### `builders.structure.parseTrajectory(blob)`

The blob overload no longer imports the mmCIF transformers. It calls `parseBlob` of the registered `mmcif` format and
fails with `parseTrajectory(blob) requires the 'mmcif' data format to be registered in this plugin.` when mmCIF is not
registered. `parseBlob` is a new optional method of `TrajectoryFormatProvider` implemented by `MmcifProvider`, and
`parseMmcifBlob(plugin, blob, params)` in `@molstar/plugin/state/formats/trajectory/mmcif` is the same function.

## Declaration contracts

`ExternalModules['jpeg-js']` exposes the injected codec's `encode` contract instead of the entire codec module type.
`JpegBufferRet` is owned by headless and describes the returned `width`, `height`, and `Buffer` data.

Standalone MVS `validateTree` accepts the logger contract `{ log: { error(message: string): void } }` rather than
requiring `PluginContext`. Existing plugin arguments still satisfy that contract.

`Clip.Rotation` is a new exported alias for `{ axis: Vec3, angle: number }`. Explicit parameter-schema annotations
retain the existing rotation value shape. The background image asset declaration changes from TypeScript `export =` to a
default export to match ESM asset imports.

Standalone MVS also exposes package-local `object`, `json`, and `color-names` helpers. The builder's `ColorNames` values
are membership flags (`true`), not the numeric colors in core's color-name table; use core when a color value is needed.
The MVS JSON expression shape remains structurally compatible with MolScript expressions, but these builder helpers do
not provide compiler semantics.

## Distribution and import paths

Compiled packages use ESM and new package subpaths. Old monolithic CommonJS and source paths require migration; use
`migration-map.json` for ownership mappings and `migration-symbols.json` for symbols whose module was split or removed.
Classic Viewer/MVS Stories asset paths and globals remain available, subject to the shape API changes above. Package
exports exclude tests and build caches.

`molstar` now packages only browser distributions. Library/server/CLI consumers must install the owning `@molstar/*`
packages instead of importing `molstar/lib/...` or `molstar/lib/commonjs/...`. There is no CommonJS build or `require`
export condition. Published library modules and the browser ESM distribution target ES2022. All 44 public packages share
one release version.

Native `gl`/`canvas` remain injected by headless callers and are optional peers. `@molstar/plugin-headless/native` adds
lazy `loadNativeModule` resolution for the MVS rendering CLI; importing it or requesting `--help` no longer eagerly
loads native modules. A rendering command still needs the caller to install them explicitly.

## Version metadata

`PLUGIN_VERSION` is generated from `version.json` rather than falling back to `'(development)'` outside a bundle. The
plugin startup log displays only this version. `PLUGIN_VERSION_DATE` and the build timestamp override have been removed.
Additional build identification metadata is deferred in the v6 checklist. The old `__MOLSTAR_PLUGIN_VERSION__`
definition does not control the generated module's value.

## CLI generated imports and entrypoints

`cifschema --moldataImportPath` now defaults to `@molstar/core/data`, replacing `molstar/lib/mol-data` in generated
TypeScript. Pass the option explicitly to generate another import path. Its field-name CSV presets read canonical root
data in the workspace and staged package assets when installed.

CLI/server entrypoints live in their owning packages and use checked-in `bin/*.mjs` launchers. Their existing command
argument parsing is retained; repository maintenance generators such as `syminfo` target the new package/test paths.
Model-server TAR helpers now expose `decodeLongPath`, `encodePax`, and `decodePax` as named ESM/TypeScript exports
instead of CommonJS `exports` assignments.

## MVSX ZIP options

`createMVSX(data, assets, options?: { zip?: ZipOptions })` accepts fflate ZIP options under `options.zip`, including
`mtime` and compression level. The default uses the current time rather than v5's fixed ZIP timestamps. Callers needing
reproducible archives should provide `zip.mtime`; use a local-calendar date when identical ZIP date fields across time
zones matter.

The return value remains an asynchronous `Uint8Array` archive with the same document/assets layout. Compression
implementation and timestamp defaults changed, so matching archive contents does not imply byte identity with v5.
`MVSData.toMVSX` currently does not forward ZIP options; the direct `createMVSX` helper exposes this control.

## Follow-up: automatic MVS asset loading

`MVSData.toMVSX` now uses platform `fetch` instead of core's `ajaxGet` task. Explicit `options.assets` still bypasses
automatic fetching.

The options argument is optional and now accepts `fetch?: typeof globalThis.fetch`. The callback receives the resolved
URI and returns a platform-compatible `Response`. Omitting it uses platform `fetch`. Cached assets are consulted before
fetching, including empty content, so a shared cache avoids repeated requests and permits exports when the source is
unavailable. Failed responses are not cached.

Under Node, platform `fetch` does not load `file://` assets. Callers can provide a file-aware `options.fetch` adapter
(for example, returning a `Response` containing bytes read with Node's `fs` APIs), or supply `options.assets` directly.
The old `setFSModule(fs)` core IO hook no longer controls standalone builder fetching. This explicit adapter contract
keeps the builder independent of plugin/rendering code. Tests cover default and custom fetching, resolved file URIs,
cache reuse, explicit assets, skipped external URIs, and HTTP errors.

## Standalone query language and MVS MolQL validation

Language modules (`language/*`, `script/*`, `transpile`, and `transpilers/*`) move from `@molstar/model/script` to the
corresponding `@molstar/query-language` subpaths. Molecular `Script` operations and `runtime/query/*` remain in model.
The text-only `@molstar/query-language/script` API provides `toExpression` and explicit language registration without
molecular imports. Existing serialized symbol IDs are unchanged.

The duplicated `@molstar/mvs-builder/expression` module is removed. Import the canonical `Expression` from
`@molstar/query-language/language/expression`. MVS authoring gains `@molstar/mvs-builder/molql`, exporting
`MolScriptBuilder` and `compileScript(language, source)` for MolScript, PyMOL, VMD, and Jmol. This facade requires no
registration imports and does not modify the plugin's language registry.

The MolScript parser now rejects repeated argument keys instead of keeping the last value. This includes collisions
between positional arguments and numeric named arguments such as `:0`.

`MVSData.validationIssues` and `isValid`, and the `mvs-validate` CLI, now reject malformed expressions, unknown callable
symbols, and undeclared argument names. The schema codecs remain vocabulary-independent. Standalone callers validating
custom symbols must provide `getMolQLSymbol(name)`, including fallback to the standard `SymbolMap` where needed. Runtime
validation resolves registered custom-property symbols and retains executable compiler checks. Missing required
arguments, argument/value types, and selection-result types remain outside standalone validation.

`MolScriptBuilder.re(pattern)` omits the absent flags argument instead of including `undefined`, so the result remains a
valid expression after JSON serialization. Passing flags explicitly preserves the previous array shape.

## Planned: with the classic UI's retirement (Prism)

Not made in v6: both remove public API, Prism ([its spec](../../packages/prism/viewer/_spec/README.md)) needs neither, and the
UI-only maps are still filled and read by the classic UI and its extensions. They land when the classic UI is retired,
each with a migration note.

- **Dead interactivity commands.** `PluginCommands.Interactivity.Structure.Highlight` and `.Select` are removed (nothing
  subscribes to them; `Interactivity.Object.Highlight` and `ClearHighlights` stay), and so is `customControl` on
  `PluginSpec.Action` (`packages/plugin/core/src/spec.ts`), a React component slot that nothing reads.
- **UI-only maps leave core `PluginContext`.** `customStructureControls`, `customImportControls` and
  `genericRepresentationControls` (React component classes and UI callbacks in a non-UI package) move to the classic
  UI's context, or go with it; Prism's contributions replace them. Extensions that fill them register with the UI they
  target instead.
- **No `molstar.Viewer` in the Prism viewer.** The Prism single file (`prism-viewer.js`) is an ES module exporting
  `createPrismViewer`, `FullPrismViewerSpec()` and the `<prism-viewer>` element; it has no classic-script global, no
  `Viewer.create(id, options)` and no `PluginExtensions` (a `Viewer` shim is not trivial: v5's options and `load*`
  signatures do not map one to one). Scripts move to `viewer.api` (the operations) or keep the classic viewer, which
  ships beside it for a major version (parity.md §9 maps each member).

Compatible core follow-ups made alongside Prism (a registry-changed event, richer `cell.stateUpdated`, settings events,
a running-task behavior, string ids for `PluginCommands`) are additive and are not listed here.

## MVS selectors through model indices: the MolQL generator `element-schema`

Tables of `StructureElement.SchemaItem` rows (`label_asym_id`, `auth_seq_id`, `beg_label_seq_id`, ...) are resolved
through the indices of the model instead of a MolScript query on every atom for every row. The public changes:

- **MolQL.** New symbols `structure.generator.element-schema` (`sel.atom.element-schema`; one set of the elements which
  satisfy at least one row), `structure.type.schema-item` (`schema-item`; a row, or a table whose fields are lists, the
  columns, with the other fields the same in every row) and `core.type.empty` (`empty-value`). `Expression.Literal`
  includes `null` (a missing value) and lists of `string | number | boolean | null`, as an argument of an application;
  the compiler passes a list through (it is the way a column of any size costs nothing), and
  `expressionValidationIssues` accepts them there (not as an expression of their own) and only for an argument which can
  be a list (declared as a list, or of any type): `atomGroups({ 'residue-test': [false] })` is still rejected.
  `formatMolScript` writes them as `[a b c]` and `empty-value`. The arguments of `element-schema` have to be
  `schema-item`.
- **`StructureElement.Schema.toExpression` returns the generator** (a table of 10 or more rows by columns) instead of
  the per-atom `atom-groups` expression, so `toLoci`, `toBundle`, `Loci.fromSchema` and `Bundle.fromSchema` use the
  indices (about 17 ms for 1000 rows on 720k atoms, several seconds before) and take no `queryContext` (no code in the
  repository passed one; `Loci.fromExpression` and `fromQuery` keep theirs, so a time limit is
  `Loci.fromExpression(structure, Schema.toExpression(schema), queryContext)`, and the generator checks it). A table
  given by columns (`items` as an object of arrays) has columns of the same length (they are rejected otherwise;
  `toLoci` and `toStructure` are the same here; they used to take the length of the first), and no columns or no values
  are no rows. `SchemaItem` has `beg_auth_ins_code` and `end_auth_ins_code`. New `Schema.toStructure` (the substructure,
  the structure itself when every element is selected).
- **Semantics differ from the previous expression of the same rows.** A residue without `label_seq_id` (ligands, waters
  in many files) never satisfies a `label_seq_id`, `beg_label_seq_id` or `end_label_seq_id` criterion (MolScript read a
  missing value as 0, so `{ end_label_seq_id: 50 }` selected them); criteria which only exist for atoms (`auth_*`,
  `label_comp_id`, atom ids, ...) select nothing on coarse elements with a warning, as annotation rows always did (the
  per-atom expression threw); a range with `beg > end` selects nothing, also on coarse elements (the expression selected
  the elements which span the two bounds); the result is one set of the elements of all rows, and the structure itself
  if all are selected. The values are matched as they are, with the type of their field: a number given as text
  (`"10"`), an identifier given as a number, or a fractional residue number on a coarse element select nothing,
  whichever criterion finds the candidates (JSON annotations convert their values when they are parsed); `null`,
  `undefined` and `NaN` are missing values (not criteria) at every entry, and `null` for a bound is an open end. An
  empty `pdbx_PDB_ins_code` means a residue without an insertion code, as in MolScript. `operator_name` (the unit
  operator's name) and `instance_id` are honored. A range with `beg_auth_ins_code` or `end_auth_ins_code` is a range of
  residues in the order of the input file, in each auth chain (the chains which share an `auth_asym_id`, e.g. a polymer
  and its ligands): from the first residue with `beg_auth_seq_id` (and that code, `''` for none) to the last with
  `end_auth_seq_id` (and its code); a bound without a code is the whole number, and it is looked for in the chains which
  are not waters first (waters are numbered again from the start in many files), and in the waters only if no other
  chain has it; a missing bound or a range which ends before it begins selects nothing from the auth chain, and so does
  a code without its number. The chains which the row selects get the part of the range which is theirs (so the waters,
  which come after the polymer, get none of a range between residues of the polymer).
- **MVS.** `SelectorParams` (`@molstar/mvs/components/selector`) is unchanged, `componentPropsFromSelector`
  (`@molstar/mvs/load-helpers`) returns a `Selector` instead of the `type` of `StructureComponentParams`, and rows and
  arrays of rows are `{ name: 'expression', params: Schema.toExpression(...) }`, so a `component` node of MVS is the
  plugin's own `StructureComponent` (the classic UI can modify it). `rowToExpression` and `rowsToExpression`
  (`@molstar/mvs/helpers/selections`) are removed: use `StructureElement.Schema.toExpression`. `prettyNameFromSelector`
  returns compact names (`A`, `A432`, `A:10-50`, `120 selections`) instead of `Custom Selection: {...}`. The dispose of
  `StructureComponent` no longer disposes the properties of a structure which is not its own (a component which is the
  source structure: the component keeps whether it is the source, as `MVSAnnotationStructureComponent` does).
- **Selection helpers move to the model package.** The resolution of rows through model indices does not depend on MVS
  and lives next to `StructureElement.Schema` in `@molstar/model/model/structure/structure/element/schema-selection/*`
  (`element-ranges`, `indexing`, `selections`, `row-structure`, `fields`); they were `@molstar/mvs/helpers/*` modules,
  which are not kept as re-exports. `IntBuffer`, `MultiMap`, `NumberMap` and `isDefined` are in `@molstar/core/util`
  (`int-buffer`, `multi-map`, `number-map`, `object`). Rows are `StructureElement.SchemaItem` there.
  `getAtomRangesForRow(s)`, `getSphereRangesForRow(s)` and `getGaussianRangesForRow(s)` take the unit's operator
  (`{ instanceId, name }`) where they took an `instanceId` string; `getAtomRangesForRowInto` writes into reusable
  buffers; `ElementRanges.union` is implemented by sorting the starts and the ends; `IndicesAndSortings.get` is shared
  by the models of a trajectory (`HierarchyScope`); `Sorting.findIndexRange` replaces a linear scan;
  `substructureFromRows` and `ElementRangesCache` are new. `MVSAnnotationRow` (annotation schema `all_atomic`) has
  `operator_name`. Annotation values given as text for numbers (and numbers for identifiers) in JSON annotations are
  converted to the type of the field.
- The annotation color theme and the annotation tooltips already used these semantics for missing values before.

## Prism UI: duplicates removed, strict command ids

- `@molstar/prism-ui`: `StatusModel` is gone (use `NotificationsManager`, same class); `shell/status.ts` is removed.
- `ShellModel.layout` and the `ShellLayout` type are gone; read `viewState.layout` (a `ViewLayout`) or save it with
  `viewState.save()`/`restore()`, and `navigation.save()` for what is open. `ShellOptions.mode` and `setMode` take
  `ViewLayout['mode']`.
- `CommandManager.register` and `registerGroups` throw on an id that is already registered (or repeated in the call),
  naming it, where the later command used to replace the earlier one silently. Override before registering.
- The viewer's Structure tab has no inline preset or style lists: `StructureTab` loses `presets`, `quickStyles` and
  `onAction`; Preset and Styles open the command menu through `openCommands`.

## Prism UI: the shell's layout

- `ViewLayout` has `hiddenTabs` (the tabs `ShellModel.hideTab` took out of view); a `ViewLayout` built by hand needs
  it (`defaultViewLayout()` has it). A kept layout without it restores with none hidden.
- The `embedded` mode is a layout of its own (shell.md §3): the compact rail, a floating right panel, nothing pinned,
  and compact below `compactWidth` (it was the desktop layout under another name).

## Seq\*: the system registry in core

- `@molstar/seq-viewer` no longer exports `SystemRegistry`: the registry is core's (`@molstar/seq-core/translate`,
  `SystemRegistry`, with no `version`); the viewer's `ctx.systems` is a `ViewerSystems` (core's registry plus the
  reactive `version`, exported by the viewer). Code that named the viewer's class names `ViewerSystems`.
- The `Coordinates` type is core's (`@molstar/seq-core/translate`); the viewer re-exports it and `isCoordinates` as
  before.
- Prism's saved sequence state (`prism-viewer.sequence`) is version 3 (adds `edits`); versions 1 and 2 still read.

## MVS loading: `UpdateTarget` builds a serialized tree (breaking)

`@molstar/mvs/load-generic`: `UpdateTarget` no longer wraps a `StateBuilder`. Its `update` field and
`UpdateTarget.commit` are removed; it has `builder: VirtualStateBuilder`, which builds the serialized state tree directly
(the same output as `StateTree.toJSON(update.getTree({ useHashVersion: true }))`, `packages/mvs/runtime/_spec/load.md`).
A loading extension that used `updateTarget.update.to(updateTarget.selector).updateState(state)` calls
`UpdateTarget.updateState(updateTarget, state)`; `UpdateTarget.apply`, `tag`, `dependsOn`, `setMvsDependencies` and
`deleteChildren` are unchanged. `UpdateTarget.create` takes an optional ref-hash cache, `loadTreeVirtual` options gain
`refHashCache` and `reuse` (`LoadTreeVirtualOptions`). `UpdateTarget.selector` is a selector on the plugin's data state
whose cell does not exist (as before: nothing is committed).

`UUID.create22` (`@molstar/core/util/uuid`) returns a per-process random prefix and a counter: still 22 URL-safe base64
characters and unique, but consecutive ids share their first 14 characters and are not random per call.

`MolstarLoadingContext.nearestReprMap` (`@molstar/mvs/load`) is `{ get(node) }`, computed on its first `get`, instead
of a `Map`; `MolstarLoadingContext.units` and `LoadingExtension.local` are new (`packages/mvs/runtime/_spec/load.md`
§4–5).

MVS-loaded transform versions (`VirtualStateBuilder`) are `jsonHash` of the params: `hashMurmur128o` for small params
(up to 256 JSON entries), a hash of the parts' hashes for larger ones; equal params still have equal versions.
`ObjectMutator.create(plan, { equals })` takes the slot equality (default `deepEqual`).

`LoadTreeVirtualOptions.resetContext` (`@molstar/mvs/load-generic`) clears what loading wrote into the context before a
load is redone without `subtreeMemo`; `subtreeMemo` is used only with it (`MolstarLoadingContext.reset` for the MVS
loader).

## Mesh primitives are typed arrays

`Primitive` (`@molstar/graphics/geo/primitive/primitive`) has `vertices: Float32Array`, `normals: Float32Array` and
`indices: Uint32Array` instead of `ArrayLike<number>`, so `MeshBuilder.addPrimitive` reads one array kind per vertex.
Every built-in producer (`PrimitiveBuilder`, `Polyhedron`, the shapes in `geo/primitive/`) already returned them; a
primitive made of plain arrays wraps them (`new Float32Array(vertices)`, `new Uint32Array(indices)`).

## Cage vertices and mesh fan/strip vertices are Float32Array

`Cage.vertices` (`@molstar/graphics/geo/primitive/cage`) is a `Float32Array` instead of `ArrayLike<number>`, and the
`vertices` of `MeshBuilder.addTriangleStrip`, `addTriangleFan` and `addTriangleFanWithNormal` are `Float32Array`, so
their loops read one array kind through `Vec3.fromFloat32Array`. `createCage` still takes any array and copies a
non-Float32Array one; a cage object literal goes through `createCage`, and fan or strip vertices are wrapped
(`new Float32Array(vertices)`).

## Unit props are read-only; units know their parent

`unit.props` (`@molstar/model/model/structure`, `Unit.Atomic`, `Unit.Spheres`, `Unit.Gaussians`) is typed `Readonly`:
the static props are shared by symmetry copies and filled only by the units' own getters
(`packages/model/_spec/unit.md` §1). Code that wrote a prop (`unit.props.bonds = ...`) passes it to `Unit.create`
instead. `Unit.Base.parent` is new: the root unit a child was made from (§2). `remapModel` drops `rings` and
`resonance` when it drops `bonds` (§3).

## Empty sequence alignment module removed

`@molstar/model/model/sequence/alignment` was an empty module (as `mol-model/sequence/alignment.ts` was in v5) and is
gone; `align` and `AlignSequences` are in `@molstar/model/model/sequence/alignment/alignment` and `.../sequence`.

## Alpha factor removed

`Representation.State.alphaFactor`, `Visual.setAlphaFactor` (and the visuals' `setAlphaFactor`) and
`RenderableState.alphaFactor` are removed: nothing set them. `uAlpha` is the base alpha (`alpha`) clamped, set with it by
the geometry's value update, no longer in `Renderable.render`; code changing `values.alpha` directly sets `uAlpha` too.

## AssetManager is an interface

`AssetManager` (`@molstar/core/util/assets`) is an interface, implemented by `DefaultAssetManager` (the former class,
same behavior), so a plugin can take another one from `PluginSpec.assetManager`. Code that constructs one writes
`new DefaultAssetManager()`; code that only takes or calls one is unchanged. `get` and `assets` return
`AssetManager.Entry` (`asset`, `file`, `isStatic`, `tag`), without the reference count; `urlCacheSize` is on
`DefaultAssetManager` only.

## Particles are opt-in: `ParticlesRegistry`

`DefaultRegistry` (`@molstar/plugin/default-registry`) no longer holds particles: no particle formats, transformers,
representations, `particles` theme scope or particle trajectory animation, and its structure and volume scopes no
longer list the `particle-*` color and size themes. `DefaultParticleRepresentations` is gone (now
`ParticleRepresentations`). `@molstar/plugin/state/transforms/catalog` no longer imports the particle transformers.
An application that reads particles composes `[...DefaultRegistry, ...ParticlesRegistry]`
(`@molstar/plugin/particles-registry`), whose import also registers the particle transformers a session needs. The
classic viewer (`apps/viewer`) composes it, so it reads particles as 5.12 did. Prism loads it with its lazy particles feature when a particle file,
an mmCIF holding particles or a session with particles needs it. `looksLikeMmcifParticles` moved to
`@molstar/model/formats/particles/mmcif-sniff` (still re-exported by `.../mmcif`).

## The component manager's `clipObjects` no longer replace a representation's clip

`StructureComponentManager` options' `clipObjects` are the clip objects the manager puts on every structure
representation: an options change replaces the objects it put there before with the new ones and keeps the others (a
representation's own, from a clipping editor or code; `RepresentationClip.replaceOwned`,
`@molstar/plugin/state/helpers/representation-clip`). It used to overwrite `type.params.clip` whole, `variant`
included; a representation with objects of its own now keeps its `variant`.

## The snapshot manager's time has its own event

`PluginStateSnapshotManager.events.changed` no longer emits when only `state.currentAnimationTimeMs` moves (every tick
of a playing transition, a scrub): `events.animationTimeChanged` does. A view that shows the time subscribes to both
(`changed` still covers a new current snapshot, which resets the time). Every reader of the entries used to wake at the
frame rate.

## `Column.__array` and `CifField.__array` are `storage()`

`Column` (`@molstar/core/data/db`) and `CifField` (`@molstar/io/reader/cif`) no longer have `__array`. Their backing
storage is the optional method `storage(): Column.Storage | undefined`: `{ kind: 'array', array }` where `__array` was
an array, `{ kind: 'dictionary', indices, values }` (`value(i)` before the schema's transform is
`values[indices[i]]`) for a BinaryCIF string column and its views (which had a string array), absent or `undefined` for
a computed column. A column literal drops its `__array: void 0`; a reader of
`__array` calls `storage?.()` and handles the kinds it knows, falling back to `value` for the rest. The storage is not
always the values: an int column over a string field keeps the strings' storage, and a BinaryCIF array storage holds
the data under masked rows, so a reader first checks `typeof value(0)` against it, as `Column.view` does. Why: the storage can
now take another form (a BinaryCIF string column as indices into its distinct strings), which a reader of an array of
values would have misread; a method also says that asking may decode (a BinaryCIF field decodes on first access).
Readers in the repository were core's `Column` helpers and `sortAtomSite`.

## Prism: each structure's sequence, track providers

- `ViewerExtension.sequence` is `{ tracks?: SequenceTrackProvider[]; behaviors?: SequenceBehavior[] }`, where it was a
  list of behaviors: an extension that brought a behavior writes `sequence: { behaviors: [b] }`. A track (data to show
  along the sequence) is better a provider (sequence.md §3): the viewer then puts it in and takes it out per structure.
- `SequenceShown.after` is optional (a document whose tracks were all removed has none; an insert without `after` goes
  last).
- Prism's saved sequence state (`prism-viewer.sequence`) is version 4: `edits` by structure ref
  (`{ [ref]: { tracks?, chains? } }`); version 3's edits read as the followed structure's, versions 1 and 2 still read.
- `SequenceViewer`'s track edits take the structure they are for (`setTrackHidden(track, hidden, of?)`), and
  `SequenceTrack` has `removable`; `sequenceMenu`'s fourth argument is whether Edit in structure tree is offered (it was
  the API: the menu has no track rows).

## MVS and Seq\* value palettes in core

- MolViewSpec's value palettes moved to `@molstar/core/util/color/value-palette` (spec → props:
  `valuePalettePropsFromSpec`, evaluation: `makeValuePaletteFunction`, named lists: `ValuePaletteLists`),
  `value-palette-params` (the param definitions, `SplitColorProp`, `MaybeFloatParamDefinition`) and `split-color`
  (`SplitColor`). Seq\*'s copy is gone; its exports (`makePaletteFunction`, `categoricalList`, `paletteHash`,
  `PaletteDefaults`, `ColorLists`, `ColorNames`, `Color`) read core's.
- MVS's runtime keeps its names as aliases (`MVSCategoricalPaletteParams`, `SplitColorProp`, `SplitColor`,
  `MaybeFloatParamDefinition`, `MvsNamedColorLists`), and the `mvs-annotation` theme's params keep their shape (saved
  sessions load). `makeContinuousPaletteCheckpoints` (`@molstar/mvs/components/annotation-color-theme`) is removed:
  `makeValuePaletteFunction` evaluates a continuous palette, without sorting the props' colors in place as it did.
- `decodeColor` reads only the color table's own names: `constructor` and other `Object.prototype` names are no colors.

## `Download`, `ReadFile` and `DownloadBlob` give Response Data

`Download`, `ReadFile` and `DownloadBlob` (`@molstar/plugin/state/transforms/data/fetch`) give
`PluginStateObject.Data.Response` (`DownloadBlob`: `Data.ResponseBlob`) unless their `asResponse` param is `false`: the
bytes, held once, and the reader chooses the form (`data.bytes`, `data.text(ctx)`, `data.json(ctx)`,
`data.read(ctx, type)`). `asResponse` unset is `true`, so a snapshot saved before, and code that does not set it, now
gets Response Data; `isBinary` matters only with `asResponse: false`. Why: text CIF parses from the bytes, never decoded
to a string and encoded back (flat-192-cif: 388 → 270 MB retained).

- Every built-in parser takes Response Data: a text parser `from: [String, Binary, Response]` (Binary Data is decoded),
  reading `SO.Data.readText(ctx, a)`; a binary one `from: [Binary, Response]`, reading `SO.Data.readBytes(a)`
  (`{ copy: true }` when it writes into its input or keeps views of it). TypeScript: a binary parser can no longer be
  applied to a ref typed `String | Binary` (`RawData`'s): a `from` that neither contains nor is contained in its type
  does not compile; such a ref is typed by what it holds, or goes through a format provider's `parse`.
- **A custom parser on top of `Download` or `ReadFile` has to adjust**: add `SO.Data.Response` to its `from` and read
  through `readText` / `readBytes`; or set `asResponse: false` on the download; or put `ReadResponse { isBinary }`
  between them. Code that reads a download's `cell.obj.data` as a string or bytes reads `data.text(ctx)` or
  `data.bytes`.
- A response body holds one form at a time (`ResponseBody`): its bytes, or its text once a reader asks for it (then
  `bytes` encodes it again). Its `header` (64 KiB) and its sniffs (`sniffBytes`, `sniffText`) work on either form.
- A `DataFormatProvider.isApplicable` gets a response as itself (a `Sniffable`, core) when the provider sets
  `sniffsResponse`, as bytes when it reads the extension as binary, else as its text, which turns the body to text
  (`DataFormatProvider.applicable`).
- Types: `builders.data.download`, `readFile` and `downloadBlob` return a `Data.Response` (`ResponseBlob`) selector,
  or `String | Binary` (`Blob`) with `asResponse: false`. The transforms' output is `String | Binary | Response`: code
  that writes it out as `String | Binary` (`StateBuilder.To<PSO.Data.Binary | PSO.Data.String>`) adds `Response`.
- Memory: a text format other than CIF holds its text, and its bytes as the file it was read from or a blob (a
  download), as before; CIF holds its bytes only.

## `Download` and `ReadFile` inflate gzip by default

`Download`, `ReadFile` and `DownloadBlob` have a `decompress` param, `'auto'` when unset: gzip data, told by its bytes and not the name,
is inflated, for text, binary and Response Data alike. A server that sends `.gz` files without `Content-Encoding` (the
EMDB maps) used to give the gzip bytes (garbage as a string); it gives the inflated data now. Code that downloads such a
file as binary and inflates it itself gets inflated bytes: it drops its own `ungzip`, or sets `decompress: 'none'` to
keep the bytes as received. `DeflateData` passes bytes that are not gzip through, so a `Download` followed by it still
works. A file named `.gz` was already inflated by its name and still is.

## CSV reader follows RFC 4180

`parseCsv` (`@molstar/io/reader/csv/parser`) reads CSV as RFC 4180 has it, from bytes or a string. Behavior changes for
input that read differently before:

- Spaces and tabs in a field are kept (they were trimmed, and inner ones dropped from the value's end); numbers still
  parse past them.
- A doubled quote in a quoted field reads as one quote (it was kept doubled); only spaces and tabs may follow a closing
  quote, other text is an error with its line, as is an unterminated quoted field.
- An empty field is empty, also the last of a record (a trailing empty field, or an empty field after a quoted one,
  moved the following values into the wrong columns). A record shorter than the header reads its missing fields as
  empty, a longer one is cut to the header's count.
- A column is a CSV column, not a CIF field: `.` and `?` are values; only an empty value is not present.
- `comment` is off by default (it was `#`), and when set, it starts a comment only at the start of a line (a value
  starting with it dropped the rest of its line).
- A UTF-8 byte order mark is skipped; the last record is counted without a final line break.

## Text parsers read bytes

The text parsers on `TextByteTokenizer` (XYZ, GRO, MOL, SDF, MOL2, PSF, TOP, PRMTOP, LAMMPS, cube, DX, OBJ, PLY, Dynamo
tbl, CSV) take `StringLike | Uint8Array`; passing a string still works (encoded to UTF-8 once). Differences:

- Fixed-width columns (GRO, MOL/SDF, PRMTOP) count bytes, not chars: a non-ASCII char in such a line shifts the columns
  after it (the formats are ASCII).
- XYZ coordinates parse with Mol*'s number parser, not `+value`: some differ in the last bit (as every other text
  format's numbers); a field that is no number reads its leading number or 0, not NaN.
- `mol/parser`'s exported `handleAtoms`, `handleBonds` and `handlePropertiesBlock` take a `TextByteTokenizer`.
- A PLY whose format line is not ascii or binary (little or big endian) fails as a string did; the plugin's bytes read
  it as ASCII.

## Direct volume `jumpLength`

`DirectVolume.Params` (`@molstar/graphics/geo/geometry/direct-volume/direct-volume`) and so the volume representation's
`DirectVolumeParams` no longer have `jumpLength` (default 0, off): volumes skip bricks of 8³ cells the transfer function
makes transparent instead, exactly and faster. The Gaussian volume keeps it (`GaussianDensityVolumeParams`,
`UnitsGaussianDensityVolumeParams`, default 4 as before): its density is made on the GPU, and its own bricks measured
slower. States and code passing `jumpLength` to a volume's direct-volume representation have it ignored.
