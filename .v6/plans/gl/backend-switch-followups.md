# Plan: Backend switching follow-ups

## Status

- **Current**: done (steps 1–4).
  - **Step 1**: `GeometryUtils.isCompatible`/`makeCompatible` (mesh: varying groups need `flatInterpolation`,
    `uniformTriangleGroup` otherwise). Visuals create incompatible geometry again in `prepareUpdate`, and their
    producers make it compatible with the new context (they keep their own capability check: isosurface and segment
    meshes don't split triangles). After a switch, Canvas3D updates every representation, not only those with dropped
    render objects; visuals with nothing to do change nothing. Making geometry compatible in place in
    `Visual.withPreparedGeometry` was tried and removed: for kept geometry it changed the mesh without the visual
    updating its values.
  - **Found on the way**: `Scene.setContext` queued the scene's render objects for adding again, and scenes that are not
    committed per frame (the camera, handle and pointer helpers) stayed empty after a switch (the axes missing). The
    render objects that were drawn are now created again right away; queued adds stay queued.
  - **Checks**: `pnpm test:viewer-smoke --switch` adds a molecular surface colored by chain after the switches and
    switches to WebGL1: no render object may keep varying groups (2 meshes made uniform on the `auto` and `webgl` pages;
    with the visual check disabled it fails with "2 meshes keep varying groups"). A WebGL2 Viewer switched to WebGL1
    renders that surface identically to a fresh WebGL1 Viewer (0 pixels differ, checked once).
  - **Step 2**: the tunnels have no context in their params. Their mesh is created with varying groups, and the shape
    representation makes it compatible: `ShapeBuilder.graphics` (the plugin's shape transformer passes the current
    context), `setGraphics`, and `makeCompatible` before creating or updating values. Unit tests
    (`repr/shape/_test/compatible-geometry.test.ts`). Not run against the ChannelsDB server.
  - **Step 3**: Canvas3DContext keeps the backend list it was created from (`create`; `fromGraphics` takes it as an
    optional argument, the WebGL convenience passes none and doesn't fall back) and falls back
    (`Attribs.contextFallback`, default on) when a context is lost for good: a restore that resolved false, no restore
    within `ContextRestoreTimeoutMs` (3 s) of visible page (paused while hidden), or more than `MaxContextLosses` (3)
    losses within `ContextLossWindowMs` (60 s), restored or not. No loops: only to backends further down the list, each
    backend at most once (lost ones and those that couldn't be created are not tried again), none while a change is in
    progress, a restore resolving after the fallback ignored, and when nothing is left one error is logged and nothing
    else happens. `contextChanged` carries `reason: 'switch' | 'fallback'` (Canvas3DContext and Canvas3D); the plugin
    logs fallbacks as warnings and emits `backendChanged`, the viewport settings show the backend in use.
    `pnpm test:viewer-smoke --fallback`: on a WebGPU page the device destroyed with device requests failing (falls back
    to WebGL2), WebGL2 lost without a restore (WebGL after the timeout), WebGL lost (nothing left: logged, exactly 2
    changes in all); on a WebGL page losses restored each time until the fourth within the window (WebGL, exactly 1
    change); each fallback renders with the coverage of the first render. The pause while hidden is not covered by a
    test.
  - **Step 4**: a restore replaces the lost context with a new one. `GraphicsContext.recreate()` returns a context of the
    same backend with the same options on the same canvas once it can be created (WebGL: after `webglcontextrestored`,
    the restored `WebGLRenderingContext` wrapped anew; WebGPU: `requestDevice` from the options, then a new context with
    that device), or undefined. This replaced the planned restore through `changeBackend`: the context knows how it was
    made, so `fromGraphics` contexts without a backend list restore too. Canvas3DContext puts the new context in place
    (`replaceContext`, shared with switch and fallback: `contextChanged` with `reason: 'restore'`, then
    `contextRestored`) and drops the lost one with `destroy({ lost: true })` (no GPU calls). Losses keep counting across
    restores, so the fallback after repeated losses still applies. The headless plugin's Canvas3DContext rebuilds the
    same way. The WebGL resources of the lost context skip their GL deletions when destroyed (`createResources`'
    `isLost`): the restored context shares the `gl` object, and deleting the lost generation's objects through it is an
    `INVALID_OPERATION`.
  - **Release**: `GraphicsContext.release()` loses the context on purpose to free its GPU resources (WebGL
    `WEBGL_lose_context`, WebGPU the device destroyed, its loss not logged as an error), `recreate` brings it back (WebGL
    asks with `restoreContext`). `Canvas3DContext.release()` / `restore()`: a released context is reported as lost but
    not restored, counted as a loss nor fallen back from until `restore`; a restore that fails leaves it released. Used
    by Prism's canvas presence (viewers out of view release their contexts within a page budget). Tested by the visual
    case `renderer/context/release-restore` (every profile with `releaseContext`, the scene rebuilt and compared) and
    on Dawn (`webgpu.test.ts`).
  - **Measurements** (4.1–4.3; Chrome, ANGLE D3D11; median of 5, loss to the first frame after the restore; the smoke
    scene and a heavy one with 300k spheres and a GPU Gaussian surface):

    | Scene         | in-place | rebuild | rebuild / in-place |
    | ------------- | -------- | ------- | ------------------ |
    | WebGL2 small  | 58 ms    | 82 ms   | 1.42               |
    | WebGL2 heavy  | 310 ms   | 360 ms  | 1.16               |
    | WebGPU small  | 58 ms    | 65 ms   | 1.12               |
    | WebGPU heavy  | 132 ms   | 158 ms  | 1.20               |

    Rebuilt images were identical to those before the loss on both backends and both scenes. The in-place restore
    rendered the WebGPU heavy scene wrong and its resource counts grew with each restore. WebGL2 small misses the 1.2
    criterion (24 ms more); rebuild was chosen anyway for correctness and the removed code. Resource counts after a
    rebuild stay flat over repeated restores (lazily created resources such as pick buffers come back on use).
  - **Removed** (4.4): `contextRestored`, `handleContextRestored` and `restore` on the contexts; the `reset` members of
    the WebGL resources (extensions, state, buffers, textures, programs, shaders, framebuffers, renderbuffers, vertex
    arrays, uniform buffers), of the WebGPU restorables and `replaceDevice` (the device is now fixed per context, the
    uniform ring has no generation); `reset` on the neutral `Texture`, `Framebuffer`, `RenderTarget`, `Renderbuffer` and
    `UniformBuffer`; the passes' and helpers' `reset`; `values.meta.reset` and the five producers' closures. Texture
    readers and `PickBuffers` keep `reset` (they cancel pending reads).
  - **Found on the way**: a texture mesh created on the new context could keep the isosurface draw arguments of the
    lost one (blank WebGPU image). `TextureMesh.DoubleBuffer.setDrawArgs` records the context of the draw arguments, and
    `get` drops arguments of another context like foreign textures.
  - **Checks**: `pnpm test:viewer-smoke --switch --fallback` on `auto`, `webgl` and `webgpu` (the fallback check counts
    backend changes, restores on the same backend are expected before it); visual tests unchanged against the
    previous commit on all profiles; the headless WebGPU test expects a new context and device after a device loss.

## Where things stand

- **Switching** keeps the render objects. Textures loaded from data are transferred. Visuals whose render object holds
  textures or storage buffers (GPU surfaces in the vertex arena) of another context create their geometry again
  (`hasForeignResources`), and Canvas3D updates their
  representations. Everything else is drawn again as it is.
- **Geometry that depends on capabilities is not caught.** CPU surface visuals (molecular, CPU Gaussian, blob,
  isosurface and segment meshes, the particle blob surface) set `varyingGroup` when the context has `flatInterpolation`
  (`dVaryingGroup`: group ids vary across a triangle). Otherwise they make the groups uniform per triangle
  (`Mesh.uniformTriangleGroup`). WebGL1 lacks flat interpolation. After a switch from WebGL2 or WebGPU to WebGL1 these
  meshes keep `varyingGroup`, so group colors bleed across triangles at group borders until the geometry is created
  again. The other direction is harmless (uniform groups draw correctly everywhere).
- **The sb-ncbr tunnels** keep the graphics context in their transform params (`graphics: PD.Value<GraphicsContext>`,
  set from `plugin.canvas3dContext?.graphics` by the actions, the behavior and the examples). They read only
  `capabilities.flatInterpolation` from it, for the same `varyingGroup` choice. After a switch they keep the choice of
  the previous context, and a context object in state params doesn't belong in snapshots either.
Before this plan:

- **A context lost for good stays lost.** Canvas3DContext restores on loss (`graphics.restore`). The WebGPU restore asks
  for a new device and resolves false if none comes (adapter gone, repeated losses). The WebGL restore waits for the
  browser's `webglcontextrestored` event, which may never come, with no timeout. Either way the canvas stays empty;
  nothing tries the other backends. Canvas3DContext also doesn't keep the backend list it was created from (`create`
  takes it, `fromGraphics` has none).
- **Two recovery mechanisms.** Restore works in place: resources reset and upload their retained data again (`reset`),
  GPU-derived values compute again (`values.meta.reset`), helpers reset. Switching rebuilds on a new context and
  transfers or recreates. Both must cover every resource kind.

## Design

### Render objects incompatible with the new context

- `GeometryUtils` gets an optional `isCompatible?(geometry, graphics)`. The mesh returns false for `varyingGroup`
  without `flatInterpolation`. It sits next to `prepare` ([geometry-preparation.md](geometry-preparation.md)): the
  geometry kind knows what it relies on.
- The visuals' `prepareUpdate` creates the geometry again when the render object holds foreign textures _or_ its
  geometry isn't compatible with the update's context. Canvas3D updates the representations of such render objects after
  a switch, as for dropped ones. Incompatible objects stay drawn meanwhile: they are only less exact, not invalid.
- The render-object check in Canvas3D works on values (`dVaryingGroup` and the context's capabilities), so it needs no
  geometry access.

### Tunnels without a context in params

- Remove `graphics` from the tunnel params. The tunnel mesh is created with `varyingGroup` (as when the capability is
  present) and goes through the same compatibility check. This needs the shape representation to take part (it has no
  context today): `Representation.setGraphics` on the shape representation, which records the context, and the
  `isCompatible` check there. Alternatively the tunnel mesh always uses uniform groups (correct everywhere, at the cost
  of the per-vertex group interpolation on capable backends). Recommended: the shape representation path, as other
  shape-based extensions can rely on it too.

### Fallback when a context is lost for good

- Canvas3DContext keeps its backend list (`create`; `fromGraphics` with an optional list) and gets a fallback policy:
  `fallback: 'next' | 'none'` (default `next`).
- WebGPU: `restore` resolving false triggers `setBackend` with the backends after the current one in the list. WebGL: a
  restore that doesn't happen within a timeout (e.g. 3 s) counts the same; a later restore event is ignored once
  switched. The same backend is tried once more first where the list has nothing after it.
- Events: `contextChanged` as for any switch, plus a `fallback` reason on it
  (`{ previous, reason: 'switch' | 'fallback' }`), so the plugin can log or toast "Graphics context lost, continuing
  with WebGL 2".
- XR sessions end on loss; nothing extra.

### Restore through the switch path (step 4)

Implemented with `GraphicsContext.recreate` instead of `changeBackend` (see Status); the rest as designed below.

**Today (in-place restore)**: the graphics context object survives the loss and every resource repairs itself.

- WebGL: on `webglcontextrestored` the context re-queries extensions and capabilities, resets its state tracking, and
  calls `reset` on every resource it created: buffers, textures (uploading the data of their last load), programs and
  shaders (compiled and linked again), framebuffers, renderbuffers, vertex arrays, uniform buffers, texture readers.
- WebGPU: `replaceDevice` swaps in a new device on the same canvas and calls `reset` on every restorable (textures,
  framebuffers, renderbuffers, render items, compute resources, uniform rings), recreating samplers, timers, mipmap and
  clear helpers; pipelines are created again on use.
- Canvas3D: `extraResets` reset the draw, pick and illumination passes (and through them SSAO, wboit, dpoit,
  postprocessing, tracing), the pick and ray helpers and Hi-Z, then `values.meta.reset` recomputes GPU-derived geometry
  (five producers: Gaussian surface and density, volume isosurface).
- Size: `reset` members in about 37 files (WebGL: 13, WebGPU: 7 plus 18 uses of the restorable set, Canvas3D passes and
  helpers: 12, neutral interfaces: `Texture`, `Framebuffer`, `RenderTarget`, `Renderbuffer`, uniforms), plus the five
  `meta.reset` producers. Every new resource kind has to implement it, and the smoke test's loss is the only check that
  it does.

**Proposed (restore as rebuild)**: a restore creates a fresh graphics context object and runs the switch path
(`changeBackend` with `reason: 'restore'`):

- WebGL: wait for `webglcontextrestored` as today (a new context on another canvas may fail while the GPU resets), then
  wrap the restored `WebGLRenderingContext` in a new `WebGLContext` (fresh resources, state, extensions) on the same
  canvas; no canvas replacement.
- WebGPU: request a device as today, then create a new WebGPU context object on the same canvas with it.
- Then as for a switch: Canvas3D rebuilds renderer, passes and helpers; render objects stay; textures loaded from data
  are transferred from their retained data; visuals with rendered or computed values (GPU surfaces) or incompatible
  geometry recreate them through their representations; image passes follow; `contextChanged` with `reason: 'restore'`.
- The previous context object is dropped without destroying its objects: they belong to a lost context (deleting them is
  pointless and on WebGL reports errors).
- The fallback policy (step 3) stays as it is: it decides between restore and fallback, the rebuild only replaces how a
  restore is carried out.

**What could then be removed** (step by step, each step checked by the smoke test's loss and the fallback cases):
Canvas3D's `extraResets` and the passes' and helpers' `reset`; `values.meta.reset` and the five producers' closures (the
representations recreate GPU-derived values, as after a switch); the WebGL resources' and WebGPU restorables' `reset`
and `replaceDevice`; the `reset` members of the neutral interfaces (a breaking change for implementers).

**Costs and risks**:

- Rebuild time: renderables are created again and representations with GPU-derived values update; programs (WebGL) and
  pipelines (WebGPU) are compiled again in both approaches. Expected similar, to measure.
- Identity: `canvas3dContext.graphics` changes on restore (today it stays). Code holding the context across a loss must
  follow `contextChanged`, as it already must for switches and fallbacks.
- Render objects without a representation that hold GPU-derived values are dropped instead of recomputed (`meta.reset`
  covers them today). Mol\* itself has none outside representations (helpers hold CPU data).
- Restore during an XR session: sessions end on loss in both approaches.

### Step 4 in detail

4.1. **Measure today's restore**: time from loss to the first restored frame, and resource counts before and after, on
both backends, for the smoke scene and a heavy one (300k spheres and a GPU Gaussian surface). Instrumented through
`contextLost`, `contextRestored` and `didDraw`; a restored image compared with the image before the loss (canvas
capture, a new smoke check that applies to both approaches).

4.2. **Prototype the rebuild** behind an attribute (`restore: 'in-place' | 'rebuild'`, default `in-place`): the restored
WebGL context wrapped anew, the WebGPU context created anew with the new device, the switch path with
`reason: 'restore'` and without canvas replacement. The smoke test's loss, the fallback cases and the image comparison
pass with both settings.

4.3. **Compare and decide**: rebuild at most about 1.2 times the in-place time and identical images, then `rebuild`
becomes the default; otherwise the prototype is dropped and in-place restore stays (recorded with the numbers).

4.4. **Remove the in-place code** (only after 4.3 chose the rebuild), in reviewable steps: Canvas3D's extra resets and
the passes' and helpers' `reset`; `meta.reset`; the backends' resource resets and `replaceDevice`; the neutral
interfaces' `reset` members (breaking changes recorded). The `restore` attribute goes away with the last step.

## Steps

1. Done: `GeometryUtils.isCompatible` (mesh: `varyingGroup` needs `flatInterpolation`), used by the visuals' update
   check, Canvas3D updating every representation after a switch; the switching check from WebGL2 to WebGL1 with a
   molecular surface colored by chain.
2. Done: tunnels without `graphics` in the params; the shape representation with `setGraphics` and the compatibility
   check; snapshots without a context object.
3. Done: fallback with the backend list on Canvas3DContext, the policy, the restore timeout, `contextChanged` with a
   reason, the plugin logging it; smoke test cases for a failing WebGPU restore, WebGL lost without restore, nothing
   left, and repeated losses.
4. Done: restore as rebuild (`GraphicsContext.recreate`), measured against the in-place restore and made the only
   restore; the in-place reset code removed.

## Decisions

- **Fallback timeout**: 3 s of visible page (`ContextRestoreTimeoutMs`); the timer pauses while the page is hidden, so a
  background tab whose restore the browser delays doesn't fall back.
- **Explicit restores wait longer**: `restore` of a released context (`release`) has 30 s of visible page
  (`ReleasedRestoreTimeoutMs`): no fallback waits for it, and under load a new WebGPU device took over 6 s (the visual
  suite's release cases on 8 parallel shards); not restored, the context stays released.
- **The plugin's Graphics Backend setting** follows a fallback: it shows the backend in use (`backendChanged`).
- **A lost backend is tried again only on the next page load**: within the page a fallback goes only down the list and
  never back; nothing about the loss is persisted, so a new page starts from the full list. Switching back by hand
  (`setGraphicsBackends`, the setting) stays possible.
