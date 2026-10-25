# Plan: Viewer app on WebGPU

## Status

- **Current**: done, milestones 1–4. Every feature WebGL renders is ported; the WebGPU visual profile passes all cases
  without skips, the Viewer smoke test passes on Chrome and Edge, and the Viewer defaults to `graphics-backend=auto`
  ([backend-selection.md](backend-selection.md)).
- The features the WebGPU programs reported unsupported are ported where WebGL renders them: palettes (`dUsePalette`,
  e.g. `volume-value` colors; `vPaletteV` per vertex, `tPalette` in the material and emissive colors) and wiggle data
  (`tWiggle`, per instance or group instance, `assign_wiggle_wgsl`) in all geometry programs, merged render items
  ([webgpu-draw-overhead.md](webgpu-draw-overhead.md)). Visual cases `geometry/palette/*` (an isosurface colored by a
  gradient volume; one-pixel band edges differ, a known difference) and `renderer/wiggle/*` (rendered at time 0 with
  animation on; they fail with the wiggle data ignored). What `getUnsupported` still reports is not produced or not
  rendered on WebGL either: grid data outside meshes (color smoothing is for meshes), solid interiors of lines, points
  and text, 2d grid textures of direct volumes (WebGPU has 3d textures).
- Milestone 1 steps 1–4 done; step 5 done: `pnpm test:viewer-smoke` (`scripts/viewer-smoke.mjs`) builds the Viewer and
  checks per backend that a local structure renders, that there are no console errors (WebGPU errors are printed in
  debug mode), and that the hover labels of a 5 × 5 grid agree between WebGL and WebGPU; WebGPU is skipped where
  unavailable (`--require-webgpu` to fail instead). Milestone 2 done (checked by hand in the Viewer: hover labels match
  WebGL at 81 of 81 grid points, highlight looks the same; visual tests `geometry/cylinders/star` and
  `renderer/marking/highlight-select` pass on WebGPU). Milestone 3 done: multi-sample compose, SMAA, wboit, dpoit, lines
  and text on WebGPU, each matching WebGL with 0 differing pixels in its visual cases; the WebGPU visual profile has 48
  passing cases (was 33); in the default Viewer only bloom is still turned off. Milestone 4 in progress: bloom done
  (`emissive` render variant of all ported geometry, luminosity, blur and composite programs, bloom in postprocessing;
  both bloom visual cases match WebGL with 0 differing pixels, 50 passing WebGPU cases), so the default Viewer runs on
  WebGPU with no feature turned off. Clipping done: pixel, instance and primitive clipping for all ported geometry, and
  the clip caps of sphere and cylinder solid interiors (`clipExit`); new visual cases
  `renderer/clipping/impostor-interior` and `renderer/clipping/points-lines`, all clipping cases match WebGL with 0
  differing pixels (54 passing WebGPU cases); mesh solid interiors (stencil caps, `depth` variant) remain. Outlines and
  shadows done: `outlines` and `shadows` programs, both composited in postprocessing; `postprocessing/shadow/default`,
  `postprocessing/outline/default` and `postprocessing/combined/occlusion-outline-smaa` match WebGL with 0 differing
  pixels (57 passing WebGPU cases). Transparent outlines and transparent occlusion need the `depth` render variant
  (transparent depth pass; until then the transparent depth stays empty, with a one-time warning), which is next,
  together with mesh solid interiors. Depth variant and mesh solid interiors done: `depth` render variant for all ported
  geometry (opaque, transparent and all-mask depth; the back-depth pass of illumination is not ported with it), mesh
  solid interior caps (stencil mark and fill passes, near-plane and clip object caps, cap picking);
  `postprocessing/outline/transparent`, the four `transparency/*/solid-interior` cases and
  `renderer/clipping/solid-interior` match WebGL with 0 differing pixels (63 passing WebGPU cases). DOF, FXAA, CAS and
  backgrounds done: `dof`, `fxaa`, `cas` and `background` programs; WebGPU cube textures (`resources.cubeTexture`, faces
  decoded without color space conversion, tracked as pending work) and mipmaps (`Texture.mipmap`, a linear downsampling
  blit per level and layer, `webgpu/mipmap.ts`); new visual case `postprocessing/background/image-blur` exercises the
  mip levels; all DOF, anti-aliasing, sharpening and background cases match WebGL with 0 differing pixels (72 passing
  WebGPU cases, the remaining 13 skips are image and direct volumes, GPU compute and illumination). Image and direct
  volumes done: `image` program (all variants but tracing; cubic interpolation, palettes, trim shapes, iso-level
  cut-off) and `direct-volume` program (color variants, as on WebGL; 3d grid textures only, the 2d grid fallback is
  reported unsupported); WebGPU textures honor `TextureImage.flipY` (rows reversed on upload, as `UNPACK_FLIP_Y_WEBGL`);
  the visual scenes get the neutral `graphics` context, so the volume and CPU surface cases run on WebGPU:
  `geometry/image/slice`, `geometry/mesh/isosurface-cpu` and `geometry/mesh/gaussian-surface-cpu` match WebGL with 0
  differing pixels, `geometry/direct-volume/blobs` with 4 (76 passing WebGPU cases; the remaining 9 skips are GPU
  compute and illumination). Headless WebGPU and device loss done: headless rendering takes a WebGPU `GPU` from the
  caller (e.g. Dawn's Node.js bindings, the `webgpu` package; `HeadlessPluginContext.create` with `backend: 'webgpu'`),
  images read with `readRenderTarget` on both backends; a lost device is replaced in place (`GraphicsContext.restore`:
  new device from the same adapter preference, per-device state recreated, textures re-uploaded from their retained
  data, render items recreating buffers, programs and bind groups on use), `Canvas3D` restores on `contextLost` with its
  pass resets. Checked in Node (image byte-identical after a restore, `plugin/headless/src/_test/webgpu.test.ts`, runs
  with the optional `webgpu` package and `--experimental-vm-modules`) and in the Viewer smoke test (device destroyed,
  new device, same coverage). Next: GPU compute (Gaussian surface, marching cubes, grid3d) as compute shaders, planned
  in [webgpu-compute-shaders.md](webgpu-compute-shaders.md) (infrastructure) and [webgpu-compute.md](webgpu-compute.md)
  (algorithm ports). Priority over the remaining WGSL coverage order of [wgsl-strategy.md](wgsl-strategy.md) step 7.
  Builds on [backend-split.md](backend-split.md) (second WebGPU milestone: `Canvas3D` and `ImagePass` on WebGPU) and
  [backend-selection.md](backend-selection.md) (descriptors, async creation). Illumination done: `tracing` render
  variant of all ported geometry (shaded color, normal with emissive, material color with density), `trace`,
  `illumination-compose` and `illumination-multi-sample-compose` programs, traced marking shading in
  `overlay`, and the sphere impostor back depth of the thickness pass; fragment-position reads of the shared depth
  textures clamp their coordinates, so the 1×1 placeholders (passes without a depth texture, e.g. the back-depth pass)
  read their texel as on WebGL instead of an out-of-bounds value; `ImagePass` renders again when a render started work
  (not only when work is still pending after it). Eight illumination visual cases (shadows, primitives, transparency,
  marking, multi-sampling, orthographic with transparent background, occlusion and outlines) match WebGL within 17
  levels (89 passing WebGPU cases).

## Goal

The Viewer (`apps/viewer`) runs on WebGPU when asked for (URL option and plugin config), loads and shows structures with
the default representations, and is interactive (rotate, hover, select). Features that are not ported yet fall back
gracefully instead of throwing; the WebGL Viewer is unchanged.

## Gap

- **Creation**: the plugin creates `Canvas3DContext.fromCanvas` (sync, WebGL) in `_initViewer`; no backend option.
- **WebGL accessors**: 29 uses of `canvas3d.webgl` / `canvas3dContext.webgl` in plugin, Viewer and extensions, mostly
  `graphics: plugin.canvas3d?.webgl` for representation contexts (a `GraphicsContext` suffices), plus screenshots
  (limits, drawing buffer size, capabilities), the animation loop timer, headless screenshots (`readPixels`), geo-export
  and extensions (sb-ncbr tunnels, alpha-orbitals) with WebGL compute.
- **Frame submission and presentation**: WebGL flushes implicitly; the WebGPU context only submits on `flush`/readback
  and draws into an offscreen drawing buffer that is never presented, and the drawing buffer size doesn't follow canvas
  resizes.
- **Viewer defaults**: SSAO (ported), SMAA anti-aliasing, temporal multi-sampling (compose), wboit transparency, marking
  pass (highlight and selection outlines), picking (hover, click). Not ported except SSAO.
- **Default representations**: cartoon and molecular surface (mesh, ported), ball-and-stick (spheres ported, cylinders
  not), Gaussian surface (GPU path is WebGL compute; the CPU fallback is used on other contexts), labels and
  measurements (text), lines/points representations.

## Milestones

### 1. Boots and renders

1. Backend option (done): `PluginConfig.General.GraphicsBackends` (descriptors tried in order, WebGL by default), the
   plugin creates the context with `Canvas3DContext.create` (async init path); Viewer option `graphicsBackend` / URL
   `graphics-backend` (`webgl` or `webgpu`, the latter falls back to WebGL).
2. `canvas3d.webgl` uses in plugin, Viewer and extensions switch to `graphics`; WebGL-only paths narrow with
   `isWebGLContext` and degrade (no GPU compute, CPU fallbacks).
3. Frame submission: `Canvas3D` ends a frame with `graphics.flush()`; the WebGPU drawing buffer renders into the canvas
   texture (`getCurrentTexture` per frame) and follows the canvas size (`getDrawingBufferSize` from the canvas).
4. Graceful fallbacks for unported features on WebGPU (done): `GraphicsContext.isProgramSupported`; wboit/dpoit and
   illumination report unsupported (transparency falls back to `blended`); `limitToSupported` turns off marking
   (material color markers instead), multi-sampling, anti-aliasing, sharpening, outlines, shadows, bloom, DOF and
   backgrounds in the `Canvas3D` and `ImagePass` props, with a one-time console note; render items of unported programs
   draw nothing with a one-time warning; picking finds nothing. The visual harness uses
   `createWebGPUBackend({ strict: true })` so unported cases stay skipped.
5. Check: the Viewer loads a structure (local file) with default representations, cartoon and surfaces render, no WebGPU
   errors; a Playwright smoke test.

### 2. Interactive

1. Picking (done): `pick` variants of mesh, spheres, points and cylinders (multiple render targets); neutral
   `GraphicsContext.readTextureAsync` and async `PickBuffers` for backends without synchronous reads (`identify` then
   finds nothing until a read has finished, clicks wait for the async pick data).
2. Cylinders (done, direct port): the WebGPU program interleaves the per-vertex attributes into one buffer when there
   are more attributes than `maxVertexBuffers`; solid interiors without clipping.
   [primitive-expansion.md](primitive-expansion.md) step 3 would remove the interleaving. Known difference: the star
   test where all cylinders overlap at the origin (depth rounding).
3. Marking pass (done): `marking` variants, `edge`, `dim` and `overlay` programs (SSAO marking shading; traced shading
   unsupported).

### 3. Viewer defaults on

1. Multi-sample compose (temporal anti-aliasing) and SMAA (edges, weights, blend; area/search textures from images)
   (done). Program names made unique for the registry (`multi-sample-compose`, `illumination-compose`,
   `illumination-multi-sample-compose`). The SMAA WGSL reproduces a quirk of the GLSL port: its
   `SMAASampleLevelZeroOffset` converts the `ivec2` offset with `float()`, taking only the first component. The SMAA
   area and search textures load asynchronously; until they are loaded, SMAA has no effect on either backend (isolated
   visual test runs of the SMAA case can render before they are loaded).
2. wboit: `colorWboit` variants and `evaluate-wboit`; then dpoit (done). `floatRenderTarget` on WebGPU follows
   `float32-blendable` (WebGL enables `EXT_float_blend` with float render targets). Textures that are not set or have no
   storage yet bind a placeholder that reads zero, as on WebGL.
3. Text (labels, measurements), lines (done); points styles already ported.

### 4. Rest

In order of what Viewer users notice: bloom (done), clipping (done; mesh solid interiors remain), outlines and shadows
(done; transparent outlines need the `depth` variant), the `depth` variant and mesh solid interiors (done), DOF,
CAS/FXAA, background (done), image and direct volume (done), illumination (done); GPU compute paths (Gaussian surface,
marching cubes) as compute shaders ([wgsl-strategy.md](wgsl-strategy.md), Compute); headless WebGPU (done); device loss
handling (done).

## Validation

- Visual test profiles unchanged on WebGL; the `webgpu` profile gains cases with each port.
- A Viewer smoke test (Playwright): load the Viewer with `graphics-backend=webgpu` and a local structure, wait for the
  first render, check for WebGPU errors and a non-uniform canvas; later hover/pick and selection checks.

## Open questions

- Default backend of the Viewer once milestone 3 is reached: `auto` (WebGPU where available) or keep WebGL? See
  [backend-selection.md](backend-selection.md), Default backend: staged, once its blockers are done.
