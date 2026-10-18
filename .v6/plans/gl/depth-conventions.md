# Plan: Depth range, reversed-Z and the rect origin

## Status

- **Current**: done. Reversed-Z on both backends, default on without a prop (steps 2, 3 and 5, in a different shape than
  planned, see "Reversed-Z as implemented"). Checked by hand on the large scenes with SSAO banding (step 1): some SSAO
  artefacts are gone, some banding remains, from the SSAO blur's depth packed into two 8-bit channels (16-bit fixed
  point, where reversed-Z gains nothing), reduced to banding far from the camera by linear depth, the rest dropped
  ([ssao-blur-depth.md](ssao-blur-depth.md)). `EXT_clip_control` on WebGL done (step 4); WebGL without it
  (`webgl2-minimum`, browsers without the extension) renders reversed without the precision gain. The top-left origin is
  done (step 6, see "Top-left origin as implemented"). Details [neutral-conventions.md](neutral-conventions.md) steps 6
  (rect origin) and 7 (depth range) and [wgsl-strategy.md](wgsl-strategy.md) step 8 (reversed-Z).

## Where things stand

- **Depth range**: the renderer builds GL projections (clip depth in `[-1, 1]`); the WebGPU global bindings convert the
  four projection uniforms to `[0, 1]` (`ProjectionConversions` in `gl/webgpu/global-bindings.ts`). Window depth is
  `[0, 1]` on both backends and holds the same values, so depth textures, the shared depth formulas
  (`perspectiveDepthToViewZ` and the like, 17 WGSL files), depth packing (33 files) and picking agree without further
  conversion. The depth range on its own buys nothing; it matters as the prerequisite of reversed-Z.
- **Depth precision**: depth increases from near (0) to far (1), and the perspective mapping crowds far distances close
  to 1, where floating-point depth is coarsest (the scene's opaque depth is `depth32float-stencil8` where available,
  `DEPTH32F_STENCIL8` on WebGL2), so precision falls off with distance; reversed-Z puts the far distances close to 0,
  where floats are finest. Near and far follow the visible scene's bounding sphere, so typical scenes are fine; large
  scenes viewed from inside (assemblies, cell packings, cryo-ET maps) are where z-fighting would show.
- **Rect origin**: the neutral API uses GL's bottom-left origin. WebGPU flips in three places (viewport and scissor,
  `readTextureAsync`/`readRenderTarget`, the texture reader); WGSL converts fragment positions with `glCoords` (24
  uses). Bottom-left is also assumed by `Camera.viewport`, picking, the trackball, interaction, XR and Hi-Z on the CPU,
  and by 77 uses of `gl_FragCoord`/`uViewport` in GLSL.

## Decisions

1. **Rect origin: top-left in the neutral API, as in the browser** (CSS, DOM, pointer events, `ImageData`, WebGPU).
   Window positions and rectangles have their origin at the top-left with y down everywhere outside the WebGL backend:
   - **Context API**: `setViewport`, `setScissor`, the readback rects (`readRenderTarget`, `readTextureAsync`, texture
     readers; rows top to bottom as today's results); WebGPU passes them through, WebGL flips with the target height
     (viewport, scissor, `readPixels`, `copyTexSubImage2D`), as WebGPU flips today.
   - **Canvas3D**: `Camera.viewport` and the viewports of the passes and helpers (101 uses), `cameraProject`,
     `cameraUnproject` and `getRay` (screen y down), picking, the trackball, interaction, XR and Hi-Z; their conversions
     from pointer coordinates (`input.height - y`) disappear, as pointer events are top-left.
   - **WGSL** works top-left without conversions: uv from the fragment position directly, `read_<name>` without its flip
     (`glCoords`, 24 uses, and the hand-written flips in the background, Hi-Z and illumination programs go);
     screen-space uv from projected positions with y down (`uv = ndc * vec2(0.5, -0.5) + 0.5`).
   - **GLSL** keeps GL's bottom-left internally (`gl_FragCoord`, 61 uses; rendered textures stored bottom-up): uniforms
     in window pixels (the background's `uViewport`, illumination's `uBounds`, the Hi-Z offset) are converted where
     WebGL binds them (`y' = height - y - h`, as WebGPU converts projections today), and screen-space uv from
     projections keeps y up there.
   - Textures loaded from images or data keep their row order (upload order) on both backends.
   - Breaking, for the API ledger: `Camera.viewport` and `Viewport` are top-left, the context's rects and readbacks too;
     consumers that converted (e.g. `height - y`) drop the conversion.
2. **Depth: reversed-Z with a `[0, 1]` clip range, on both backends.** The neutral convention becomes WebGPU's (the
   renderer builds `[0, 1]` reversed projections, near → 1, far → 0); WebGL adapts: `EXT_clip_control` with
   `ZERO_TO_ONE` where available (WebGL2), otherwise the inverse of today's conversion on the projection uniforms
   (reversed `[-1, 1]`, correct but without the precision gain, as WebGL1 and WebGL2 without the extension today).

## Reversed-Z as implemented

- **Projections**: the camera keeps building GL projections (`Mat4` in double precision); uniform specs holding a
  projection or its inverse are tagged (`UniformSpec.depthRange`: `'projection'`, `'inverse'`; the global projections
  and those of the SSAO, outline, DOF, shadow and tracing passes) and each backend converts them in double precision
  where it binds them (`gl/depth-range.ts`: `splitValues` and the global bindings): WebGPU to `[0, 1]` reversed
  (`z' = (w - z) / 2`, replacing `ProjectionConversions`), WebGL to `[-1, 1]` reversed (`z' = -z`). The CPU (frustum,
  `cameraUnproject`, Hi-Z sphere bounds) keeps GL conventions; the picked depth is converted back
  (`PickBuffers.getDepth` returns `1 - d`).
- **Window depth** is `1 - d` of the GL depth `d` on both backends: the background and the cleared depth 0, packed far
  depth zeros (clears of packed depth targets, the renderer's placeholders, the marking depth).
- **Render states**: `greater` for the scene, `less` for the back-face (thickness) pass, `greater-equal` for occluded
  emitters, solid interior caps and fill depth, depth clears 0 (the thickness pass 1, with packed ones).
- **Shaders, both languages**: the depth formulas for reversed depth, background tests on 0, flipped compares and
  min/max (outlines, postprocessing, compose, bloom, DOF, marking, emissive, depth variants, Hi-Z reduction), wboit
  weights from the reversed depth, dpoit peeling on `1 - depth` (its blending assumes increasing depth), impostor
  near-plane caps and clip markers, vertex near/far-plane tricks (direct-volume, mesh solid interiors); WGSL passes do
  their screen-space math in `[0, 1]` depth (`screenSpaceToViewSpace`, ray offsets).
- **WebGL clip control**: with `EXT_clip_control` (Chrome exposes it on WebGL1 too) the context sets `ZERO_TO_ONE` once
  and again after a restore, binds the WebGPU conversion (`WebGLContext.depthRangeConversion`) and compiles GLSL with
  `clipZeroToOne`: `clipZToDepth`/`depthToClipZ` (`calcDepth`, `screenSpaceToViewSpace`, the tracing ray offset), the
  far-plane vertex tricks (direct-volume, mesh solid interiors), the line near-plane estimate and the clip z of the GPU
  Gaussian density and color smoothing (were `[-1, 1]`, clipped half the atoms). WebGL images change by single pixels;
  the drawing buffer's 24-bit depth gains nothing, the scene's float depth targets do.
- **Cylinder impostors**: box corners in front of the near plane are clamped onto it (were moved back by twice the
  length, which could cross the far plane and drop a cylinder whose tip touches the near plane;
  `geometry/cylinders/star` lost one cylinder on both backends, its known difference is gone).
- **Checks**: the visual suite passes on all profiles; WebGL images against the previous ones differ in the cylinder
  star (the fix above) and in a few pixels of SSAO, outlines and illumination (depth quantization in packed channels);
  the Viewer smoke test checks that picked positions lie at the picked element; the render benchmark with `--gpu-timing`
  shows no runtime cost.

## Top-left origin as implemented

- **Context API**: rectangles top-left on both backends; WebGPU passes them through (its flips went), WebGL keeps the
  viewport and scissor as given and applies them flipped with the height of the render pass target (framebuffers learn
  it from their attachments), flips readback rectangles and the texture reader's rows, and the source and destination of
  `copyFromFramebuffer`. Its own compute on data textures converts rectangles in texel rows from the bottom internally
  (`setTexelViewport`, `setTexelScissor`; the pass target height is internal too); the Gaussian density and color
  smoothing atlases place tiles by texel rows from the bottom explicitly (unchanged layouts on both backends).
- **CPU**: `Camera.viewport`, project/unproject/rays, picking (no pointer flip; pick buffers top to bottom), trackball
  and interaction checks, XR (the viewports of its WebGL layers converted from GL's origin), the image pass crop, Hi-Z
  (levels anchored at the top-left of the viewport, indexed top to bottom), the viewport params.
- **Uniforms**: `UniformSpec.windowRect` tags (`uBounds` of SSAO, shadows, DOF and tracing; `uViewport` of the
  background and SMAA; the Hi-Z `uRect`; `uOcclusionOffset`), flipped by WebGL where it binds them.
- **WGSL** top-left: `fragmentUv`, reads without flips, projected uv with y down, `screenSpaceToViewSpace` with y down,
  the SSAO normal's "up" neighbor, the DOF box offsets, the background gradient, Hi-Z; SMAA keeps GL's window space
  internally (oriented searches and area lookups), and the SSAO noise and tracing seeds use bottom-up coordinates so
  that they match WebGL's.
- **Checks**: WebGL byte-identical on all profiles, WebGPU passes all cases, the smoke test's picking (positions at the
  picked element) on both. Partial viewports covered since by cases captured from the canvas (`capture: 'canvas'`, a
  Canvas3D frame read in the task of the frame, as the image pass draws the full size): `renderer/viewport/canvas`,
  `relative-frame` and `static-frame` (frames placed from the top-left), `renderer/camera/stereo` (side by side) and
  `renderer/hi-z/viewport` (culling in a frame; the same image as without Hi-Z on WebGL2 and WebGPU, checked once); XR
  is not covered.

## What reversed-Z touches

- **Projections and the CPU**: `Mat4.perspective`/`ortho` variants for the reversed `[0, 1]` range (the camera builds
  them), `Frustum3D.fromProjectionMatrix`, `cameraUnproject` and the picked position (depth read back from the pick
  target), Hi-Z culling (the depth bounds of a sphere).
- **Render states**: depth compares flip (`less` → `greater`, `less-equal` → `greater-equal`, the back-face pass
  `greater` → `less`; 5 states in the renderer), depth clears 1 → 0 (16 places).
- **Shaders, both languages**: depth linearization and reconstruction (`perspectiveDepthToViewZ`, `depthToViewZ` and
  their inverses), fragment depth written by impostors (spheres, cylinders: 13 files), background tests on depth 1
  (`isBackground` in 18 files, about 22 comparisons with 1.0 or 0.99), depth compares in SSAO, outlines, DOF, fog,
  marking, wboit and dpoit (depth peeling blends min and max of depths), the Hi-Z reduction (`max` → `min`), depth
  packing where the background value matters.
- **Formats**: precision gains need float depth; the WebGPU drawing buffer's `depth24plus-stencil8` and 16-bit depth
  textures (`ushort`, the renderer's placeholder) gain nothing but must still work.

## Steps

1. **Measure**: by hand, on the large scenes with SSAO banding (before and after step 3); a large-scene visual case
   (e.g. 60 instanced copies of 7qpd on a ring, camera inside) to keep it checked.
2. **Projection helpers**: reversed `[0, 1]` perspective and orthographic matrices with unit tests (the near plane at
   depth 1, far at 0, the inverse), the frustum extraction and unprojection for them; the depth formulas as one shared
   neutral definition (GLSL and WGSL chunks generated or kept in sync by a test).
3. **WebGPU first** (behind a renderer prop, default off): reversed projections without `ProjectionConversions`, the
   flipped compares and clears, the shader changes above; the visual suite on WebGPU with the prop on against WebGL2
   (identical within today's tolerances), the step 1 case shows the gain. The render benchmark with `--gpu-timing` with
   the prop off and on: no runtime cost expected (only the compare direction and the clear value change).
4. **WebGL**: `EXT_clip_control` where available, the projection conversion otherwise, the same shaders; the visual
   suite on all WebGL profiles with the prop on.
5. **Default on**, the prop removed; picking, Hi-Z, OIT, illumination and the Viewer smoke test checked; the conversion
   in the WebGPU global bindings gone.
6. **Top-left origin** (independent of 2–5, after them, or before if it eases the shader changes): the context API first
   (WebGL flips, WebGPU's flips go), then `Camera.viewport` and the CPU code, then WGSL without `glCoords` and the GLSL
   uniform conversions; the visual suite on all profiles unchanged (WebGL byte-identical), picking (the smoke test's
   hover and identify grid), screenshots and the readback checks (`gpu-paths/*`) unchanged; the breaking changes
   recorded.

## Open questions

- Whether the drawing buffer on WebGPU should use `depth32float-stencil8` where available (today `depth24plus-stencil8`)
  so that scenes drawn directly to the canvas gain too.
- Custom shaders of extensions or users that compare depth with 1.0: a breaking change to record, with the depth helpers
  as the migration.
