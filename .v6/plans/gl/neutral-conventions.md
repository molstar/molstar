# Plan: Neutral conventions follow WebGPU, WebGL adapts

## Status

- **Current**: steps 1–7 done. 7 (depth range) with reversed-Z: window depth is reversed on both backends, the
  camera's projections stay GL projections on the CPU and each backend converts the tagged projection uniforms to its
  clip range. 6 (rect origin): top-left in the neutral API, WebGL flips at its boundary. Both in [depth-conventions.md](depth-conventions.md) (decided: the neutral
  origin becomes top-left as in the browser, WebGL adapts). The readback API landed with bottom-left rects
  ([backend-narrowing.md](backend-narrowing.md) step 4), so step 6 is no longer tied to it. Viewports partly outside the
  target are cropped on both backends (current WebGPU spec), scissors intersected with the target
  (`intersectTargetRect`).
- Steps 1-7 done. Related: [wgsl-strategy.md](wgsl-strategy.md) (decision 3 on depth range, step
  6 prototype), [backend-split.md](backend-split.md) (step 5 readback API), [texture-bindings.md](texture-bindings.md)
  (decision 8, alpha textures), [minimum-capabilities.md](minimum-capabilities.md).

## Principle

WebGL peculiarities and legacy are worked around in the WebGL backend, not in the WebGPU backend or the neutral API.
Where WebGL can adapt (formats, defines, coordinate flips, clears), the neutral API uses WebGPU semantics and the WebGL
backend converts. WebGL-only fallbacks (packed depth, emulated vertex ids, WebGL1 paths) stay behind capabilities or
inside `webgl/`.

## Inventory

Workarounds of WebGL conventions found in the WebGPU prototype and the neutral API:

| Item                 | Today                                                                                                                                                                                                                | Proposal                                                                    |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `rgb` formats        | 8-bit `rgb` for `tColor`, `tPalette`, `tSize`, color smoothing, SMAA area; `rgb` float in texture-mesh specs. WebGPU has no 3-channel formats and expands to `rgba` on every upload (`expandRgb`)                    | Drop `rgb` from `TextureFormat`, producers write 4 channels                 |
| `alpha` formats      | WebGL1 `ALPHA`, WebGL2 `ALPHA` (ubyte) or `RED` (float/int); GLSL reads `.a` (18 chunk sites) or switches on `dValueChannel` (marching cubes); `resources.getValueChannel`; WebGPU stores `r8unorm`, WGSL reads `.r` | One-channel textures in red everywhere, format renamed `r`                  |
| Vertex id emulation  | `aVertex` has `usedIf: ['!noNonInstancedActiveAttribs']`, a WebGL context define; the WebGPU program sets that define to drop the attribute                                                                          | A positive WebGL backend define for the emulation                           |
| `tDepth` placeholder | Renderer binds a 1x1 `rgba` texture to a slot whose spec is `depth` (packed depth fallback of WebGL1 without depth textures); WebGPU can't bind it to a `texture_depth_2d`                                           | Placeholder of the spec format, packed only where `depthTexture` is missing |
| Clears               | `clear` respects the scissor and can come at any point (GL semantics); WebGPU emulates with load ops, partial clears unsupported                                                                                     | Clears are pass load ops; partial clears audited                            |
| Rect origin          | `setViewport`, `setScissor`, `readPixels`, `copyFromFramebuffer` use bottom-left; WebGPU flips (`toTargetRect`), image readback flips (`PixelData.flipY`)                                                            | Consider top-left (step 6)                                                  |
| Depth range          | Camera and renderer build `[-1, 1]` projections; WebGPU converts four matrices in its global bindings                                                                                                                | Consider `[0, 1]` neutral (step 7)                                          |

Not WebGL cruft, kept as is:

- std140 padding of uniform blocks: WGSL uniform layout rules require the same 16-byte array strides.
- Float `aGroup`/`aInstance` attributes: a transport choice (WebGL1 has no integer attributes); WebGPU may later use
  `instance_index` and integer data where it pays off.
- WebGL-only fallbacks behind capabilities (packed depth, WebGL1 paths): WGSL ignores them (wgsl-strategy decision 1).

## `rgb` (step 1, done)

- Done: color, palette and size data with 4 channels (colors with alpha 255), merged items, mesh and texture-mesh color
  smoothing, SMAA area texture, alpha-orbitals coefficients (only 2 of 3 channels were used), texture-mesh specs,
  geo-export strides; `rgb` removed from `TextureFormat`, both backends and WebGPU's expansion. Visual tests
  byte-identical on all WebGL profiles.

- Hardware rarely has 3-channel 8-bit formats: D3D11 (ANGLE on Windows), Metal and most Vulkan drivers have no sampled
  `RGB8`, so ANGLE allocates `RGBA8` and converts every upload with a CPU copy. GPU memory is the same with `rgba`, and
  uploads skip that conversion; only the CPU-side arrays grow by a third.
- `RGB32F`/`RGB16F` are not color-renderable in WebGL2 either; the texture-mesh specs (`tPosition`, `tNormal`, kind
  `texture`) only describe textures that compute allocates elsewhere, so their spec format becomes `rgba`.
- Producers: `color-data.ts` (`createTextureImage(..., 3, ...)`, `Color.toArray(..., i * 3)`), `size-data.ts` (24-bit
  size packed into rgb, unchanged packing, alpha unused), palette arrays, `merged.ts` (`itemSize: 3` for
  `tColor`/`tSize`), mesh and texture-mesh color smoothing, SMAA area texture.
- Consumers reading the arrays: `geo-export/mesh-exporter.ts` (`tSize`, `tColor`, palette strides), `size-data.ts` max
  size, merged copies. Shaders read `.rgb` and are unchanged.
- Alpha written as 255 so the data is a valid opaque color.
- Then `rgb` leaves `TextureFormat`, `getFormat`/`getInternalFormat`, byte counts, and WebGPU's `expandRgb`.

## Single-channel textures in red (step 2, done)

- Done: format `r` (`R8`/`R16F`/`R32F`/`R32I` on WebGL2, `LUMINANCE` on WebGL1), GLSL reads `.r` (marker, transparency,
  emissive, clipping, wiggle and grids, direct-volume transfer, image, text font, marching-cubes triangle tables,
  alpha-orbitals without its WebGL1 switch). `getValueChannel` is gone; `dValueChannel` stays because the channel is a
  property of the data, not of the backend: Gaussian densities keep the value in `.a` with the group packed into `rgb`,
  so callers of `extractIsosurface` pass the channel. `readAlphaTexture` became `readSingleChannelTexture` (values in
  `.r`). Fixed on the way: `image` read its float `alpha` texture from `.r`, which was 0 on WebGL1; `tGroup` of texture
  meshes was specified `alpha` but holds a packed `rgba` group. Visual tests byte-identical.

- WebGL2: `R8`/`RED` for ubyte (now legacy `ALPHA`), as already for float and int.
- WebGL1: `LUMINANCE` instead of `ALPHA` (ubyte, and float/half-float with the texture float extensions); sampling
  returns `(L, L, L, 1)`, so `.r` reads the value as on WebGL2 and WebGPU.
- GLSL reads `.r`: marker, transparency, emissive, clipping, wiggle and the grid variants (18 chunk sites), marching
  cubes without `dValueChannel`; `getValueChannel` and `dValueChannel` go away.
- Render targets with format `alpha` already fall back to `rgba` on WebGL1; audit that writers put the value in red
  (SSAO, hi-z, others creating `alpha` targets) and that readback (`pixelPack('alpha', ...)`) reads red.
- Rename the format `alpha` → `r` (WebGPU channel naming: `r8unorm`, `r32float`), `rg`/`rgba` unchanged.
- Validation: visual tests byte-identical on all profiles, including `webgl1-minimal`.

## Vertex id emulation (step 3, done)

- Done: `emulateVertexId` in the GLSL ES 1.0 vertex prefix and in the 3.00 prefix where `gl_VertexID` is unusable,
  `aVertex` used if it is defined; the WebGPU program has no backend define set. Visual tests byte-identical.

- `aVertex` gets `usedIf: ['emulateVertexId']`, defined by the WebGL context define lookup when `gl_VertexID` is not
  usable (WebGL1, the Safari active-attribute bug); GLSL picks `VertexID` from the same define.
- The WebGPU program loses its `BackendDefines` set.

## Depth placeholder (step 4, done)

- Done: `tDepth` held two kinds of data, the opaque scene depth (read with `getDepth`) and packed depth with alpha
  (front transparent depth for emissive, the marking depth, read with `unpackRGBAToDepthWithAlpha`); WebGPU can't bind
  both to one `texture_depth_2d` slot. The packed data moved to a new global `tDepthAlpha` (`rgba`, used by the emissive
  and marking variants), `tDepth` is only scene depth. The renderer's `tDepth` placeholder is a 1x1 depth texture
  cleared to 1 through a framebuffer pass where `capabilities.depthTexture`, the packed `rgba` placeholder otherwise;
  `tDepthAlpha` keeps the packed one.
- WebGPU gained framebuffers (attachments for passes; `attachFramebuffer`/`detachFramebuffer`, depth-only targets,
  sampling views with the depth aspect of depth-stencil formats) and deferred destruction: resources destroyed while
  commands are recorded are released after the next submit (WebGPU invalidates destroyed resources for all later
  submits, GL drivers defer deletes).
- Visual tests byte-identical on all WebGL profiles, `webgpu` cases unchanged (the placeholder clear runs in every
  case).

- The renderer's empty `tDepth` matches its spec: a 1x1 depth texture cleared to 1 where `capabilities.depthTexture`,
  the packed `rgba` texture only without it (the packed path is already a capability decision of the passes).
- Prerequisite for the WebGPU depth, marking and emissive variants, which declare `tDepth` as `texture_depth_2d`.

## Clears (step 5, done)

- Audit (45 sites): compute passes clear right after `beginRenderPass` with a scissor over their target; postprocessing
  passes (CAS, FXAA, SMAA, DOF, compose, marking) and the renderer clear within the camera viewport, which is smaller
  than the target for stereo eyes and viewport exports, so partial clears are a real need, not WebGL legacy; two clears
  happen between draws (solid-interior stencil, OIT depth resolve).
- Done: the neutral semantics stay "clear within the scissor rectangle, ignoring write masks" (documented on
  `GraphicsContext.clear`; `RenderPassOptions.clear` is `clear` right after `beginRenderPass`). WebGPU uses load
  operations for whole-target clears before the first draw and draws a full-screen triangle with a clear pipeline
  (color, depth via position z, stencil via `replace`) for all other clears, so no clear is unsupported anymore.
  Verified by forcing every clear through the draw path: `webgpu` cases unchanged (the stencil part is not exercised
  yet).
- Not done: moving clears into `beginRenderPass` options. Most sites set viewport and scissor between `beginRenderPass`
  and `clear`, so the move would apply the previous scissor; WebGPU gets load operations without it.
- Original notes:

- Audit the 19 `clear` call sites (renderer, passes) for clears after draws and clears with a scissor smaller than the
  target (viewport rendering, stereo eyes, tiled export).
- Full-target clears move to `beginRenderPass(target, { clear })` where they open a pass; `clear` stays for clears
  inside a pass and starts a new pass on WebGPU.
- Partial clears, if any remain: a neutral clear-quad helper drawing through the normal pipeline on both backends, not a
  scissored `gl.clear` that WebGPU can't express.

## Rect origin (step 6, consider)

- WebGPU, the DOM, `ImageData` and pointer events use a top-left origin; only GL uses bottom-left. With top-left neutral
  rects, WebGL flips in `setViewport`, `setScissor`, `readPixels` and `copyFromFramebuffer` (it knows the target
  height), and the flips in image readback and in the WebGPU backend disappear.
- Scope question: `Camera.viewport`, picking and trackball (`input.height - y`), hi-z projection and the `uViewport`
  uniform used with `gl_FragCoord` in GLSL also assume bottom-left. Converting only the context API keeps those as they
  are, with one conversion where canvas3d calls the context; converting the camera viewport is a larger, separate
  change.
- Best done together with the async readback API of backend-split step 5, which changes the `readPixels` signature
  anyway.

## Depth range (step 7, consider)

- wgsl-strategy decision 3 keeps `[-1, 1]` on WebGL and converts on WebGPU. The inverse is possible: the renderer builds
  `[0, 1]` projections (WebGPU, D3D, Metal, Vulkan convention) and WebGL converts, or uses `EXT_clip_control` with
  `ZERO_TO_ONE` where available (`capabilities.clipControl` exists).
- `EXT_clip_control` would also enable reversed-Z on WebGL2 (wgsl-strategy step 8) with better depth precision; the
  CPU-side frustum extraction (`Frustum3D.fromProjectionMatrix`) and depth reconstruction in GLSL would follow the
  chosen convention.
- Not started; revisit when reversed-Z is planned.

## Validation

- Each step: visual tests byte-identical on all WebGL profiles (or documented differences), the `webgpu` cases
  unchanged, graphics jest and the full type check.
- Steps 1 and 2 need no guard test: `TextureFormat` and `GraphicsResources` no longer have `rgb`, `alpha` or
  `getValueChannel`, so the type check rejects them.

## Open questions

- Step 6 scope: context API only, or also the camera viewport and shader uniforms?
