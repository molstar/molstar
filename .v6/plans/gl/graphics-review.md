# Plan: follow-ups from the graphics layer review

## Status

- **Current**: batches A and B done, B verified in CI (the Dawn job and the visual smoke run passing on the runners);
  batch C (the multi-sample savings, measured first) done: P1, P2 and P3, P5 kept as is; batch D done: S5-S8,
  W1-W3, P6-P9 done or decided, S2 and S3 done as far as they cut. The review is closed. Follow-ups
  outside the review: [geo-split.md](geo-split.md) and [core-gpu-compute.md](core-gpu-compute.md), both deferred.
  - **S2 done**, `Renderer.create` split by its seams into `gl/renderer/`: the pass tables (`pass-draws.ts`: per
    pass its checks, draws and render states, solid interior states; module-level, no closure state) and the draw lists
    (`draw-lists.ts`: snapshots, chunks, record and replay, behind `DrawListHost`, what they read of the renderer, its
    frame state by getters) and the solid interior caps (`solid-interior-caps.ts`: the cuts of the near plane and clip
    objects, per segment, marked and filled, behind `SolidInteriorCapsHost`) and GPU culling with the sphere raster
    (`gpu-culling.ts`, one module: each reads the other, the raster's split culled with and its second phase run with
    the culling's; the opaque passes drive it by named steps, `resolveRaster`, `beginOpaque`, `finishRaster`).
    `renderer.ts` 2,924 to 1,668 lines, the allocation guards unchanged. Left in the closure: visibility and the ghost
    cull (they set the bounding sphere and facts every draw reads), the uniforms' update, the passes' entry points.
  - **S3 done**, `Canvas3D.create` split the same way: its modules take one `Canvas3DParts` (`canvas3d/parts.ts`:
    the canvas's parts, getters for those a context change replaces), not a host each. First the preparation
    (`canvas3d/prepare.ts`: `prepare`, `precompile`, `prepareObjects`, the warm-up, the pass targets' release).
    Then the representations' membership (`canvas3d/representations.ts`: each representation's render objects in the
    scene of its layer, kept in step, moved with the layer, dropped by a context change; the mark index) and the debug
    stats (`console-stats.ts`), the marking (`marks.ts`). Then the camera reset (`camera-reset.ts`: requested by the
    API or a commit that moved objects out of view, resolved after the whole scene committed) and what a frame draws
    (`frame-render.ts`: the marking-only redraw, illumination or the draw pass, the Hi-Z levels, more marking
    samples; when to draw stays, its flags shared with the scheduling), its marking redraw and occlusion callback
    made once (hover 10.4 kB a step on WebGL, 15.8 on WebGPU, from 10.5 and 15.9). Then the parts bound to the
    context (`context-parts.ts`: Hi-Z, renderer, shader manager, pick and ray helpers, multi-sample helper), made by one
    function at the start and on a context change (`setContext` made them a second time, by hand). `canvas3d.ts` 3,077
    to 2,200 lines. S3 done as far as it cuts; left in the closure on purpose: picking (thin over its helpers), frame
    scheduling (`tick`, `draw`, `animate`) and the rest of `setContext` (both reassign the closure's state: a module
    would turn it into a host of setters).
  - **P6 done**: a WebGPU texture reader (picking) copies the textures of a read into one staging buffer, each at its
    offset, in one command buffer and submit (four of each before), maps it once, and copies its rows straight into
    the caller's arrays (an intermediate array per texture before: about 2 MB a dirty pick at 480×270); its staging
    buffers kept for the next reads, two at most while a pick each frame waits for its map. Test:
    `gl/_test/webgpu-texture-reader.test.ts` (the stub device's calls: no buffer made again, one encoder a read).
  - **P8 done**: WebGL records what each attachment point of a framebuffer holds (by GL framebuffer: its GL texture or
    renderbuffer and layer, `webgl/framebuffer.ts`, `isAttached`/`setAttached`); attaching the same again returns
    before ending the open pass or any GL call, so the driver checks no framebuffer again (passes attach their targets
    every frame). Any change attaches: another texture, point or layer, a renderbuffer on a texture's point, a detach;
    the combined depth-stencil point clears the separate ones and back, as the storage tracking. Visual suite on WebGL2
    (291) and WebGL1 (290) unchanged; `gl/_test/webgl-attachments.test.ts` counts the GL calls (medium: native `gl`,
    the shim's GL objects are numbers, which the tracking leaves out).
  - **W1 done**: framebuffers are internal to the backends. A pass names what it draws into as a whole, as WebGPU's
    render pass descriptor: `RenderPassAttachments` (`gl/render-pass-attachments.ts`, `color` textures or 3d texture
    layers in draw buffer order, `depth` a texture or renderbuffer), made by `resources.renderPassAttachments` and kept
    by the pass (not a literal per begin: no allocation a frame), its fields changed between passes (ping-pong) and
    destroyed with the pass. Each backend draws one with a target of its own (WebGL a framebuffer, WebGPU an attachment
    set, `webgpu/attachment-set.ts`), brought in line on begin (`createBackendRenderPassAttachments`; attaching what a
    point holds is skipped, P8). `RenderTarget` is a texture and a depth renderbuffer with the attachments naming them
    (no framebuffer, no `bind`). The draw pass names the depth at each begin (`pass(target, depth)`: its own, none,
    the opaque depth texture, or packed, the color target's depth-stencil), no record of what earlier passes left;
    five begins that only changed the depth of a pass nothing drew in are gone. `Texture.copyFromTexture(source, …)`
    replaces `copyFromFramebuffer` (WebGL reads through a framebuffer the source keeps). The WebGL compute keeps its framebuffers (`WebGLResources.framebuffer`,
    `WebGLContext.beginRenderPass` takes one). Kept from the design: WebGPU's attachment set target cache (a pass into
    unchanged attachments makes nothing). The pass targets record what passes draw into by object, not by key strings
    (`forEachRenderPassAttachment`), and WebGL's attachment refs are made once: the hover 10.5 kB a step on WebGL
    (10.6 before), 15.9 kB on WebGPU (17.2), the frame cases unchanged. Visual suite unchanged (291/290/296).
  - **S8 done**, steps 1 and 2: the context's named caches (`namedTextures`, `namedFramebuffers`,
    `namedComputeRenderables`, string-keyed dictionaries anyone wrote) are gone; a module keeps its scratch resources in
    a `createScratch<T>(label)` of its own (`gl/scratch.ts`: a `WeakMap` per context, by name within the module, so
    no two modules collide), 11 modules (the WebGL compute, the GPU Gaussian density, the texture mesh color smoothing,
    the Gaussian surface's density textures, the volume segment). Every one was an intermediate a computation fills and
    reads (none a result callers keep), so the lifetime is the context's, as before. Step 3, idle release (as the scratch
    buffers have), decided against for now: scratch resources stay for the context's lifetime. If revisited: per module
    (a cached compute renderable holds the textures it was made with), after measuring what they hold after a session
    (the Gaussian surface's density textures the largest).
  - **S7 done**: the folders of graphics in layers, lowest first `gl`, `theme`, `geo`, `formats` and `repr`,
    `canvas3d` and `props`, `util` (`src/_test/layers.test.ts`: a folder imports only from lower ones, relative imports
    included). `theme` sits below `geo` (geometries take a theme), so what the themes took of `geo` moved to them (the
    color and size types; `isPositionLocation` from model). The renderer reads a camera of `gl`'s own (`gl/camera.ts`,
    `RenderCamera`), `Viewport` from core; the marching cubes tables in core (`geo-split.md` step 1);
    `getVolumeTexture2dLayout` in `geo/util` and `ShapeGetter` in `geo/shape/provider`, not `repr`; `namedCatalog` in
    core, `util` the top layer (gizmos and their visual sets).
  - **S6 done**: the backend folders' exports source-only (`molstar-src`): the published package doesn't export
    `gl/webgl/*` and `gl/webgpu/*`, the repository (its harness, tests and builds, all on `molstar-src`) still reaches
    them. Contexts through `gl/contexts/webgl` and `gl/contexts/webgpu` (factories, guards, types, backend descriptors),
    what an extension uses through `gl/extend/webgl` (`registerGlslProgram`) and `gl/extend/webgpu`
    (`registerWgslKernel`, and `registerWgslProgram`, new: an extension's program drew on WebGL only before). The GLSL
    sources, chunks and assembly stay internal (an extension brings its own GLSL). The registries in modules of their own
    (`webgl/program-registry.ts`, `webgpu/shader/program-registry.ts`, `kernel-registry.ts`), the built-in tables
    registering into them: `gl/extend/*` reaches only the registries (`extend-imports.test.ts`), so registering an
    extension's WGSL bundles none of WebGPU into a WebGL-only build, nor its GLSL any of WebGL into a WebGPU-only one. A Biome rule keeps in-repository
    code outside graphics, tests and `tests/` out of the backend folders; `backend-neutral.test.ts` needs no allowlist.
    `breaking-v6-changes.md` (Backend folders internal).
  - **P9 done**: a WebGPU item's draw block made by a function of its own, once (`createDrawBlock`), not through two
    closures a draw: a direct frame of 8 objects (draw lists off) 5,273 to 4,249 bytes; the mark buffer emptied in place
    each render, not replaced. The classic uniforms' `.slice()` went with S5. Guarded: `frame: opaque pass drawn directly`
    (`frame.quick.ts`, draw lists off), 128 bytes on WebGL, 4,249 on WebGPU.
  - **P7 decided, kept**: a hover pick is made only when the pointer moved (`handleMove`), so a spin or rock with the
    pointer still draws none; with the pointer moving, a pick a tick (30 Hz) for a pose one or two frames gone, nearly
    right. Skipping them while `controls.isAnimating` would stop hover labels and highlights during a slow rock, for
    little saved.
  - **The representations' move broke a Prism test's timing** in the whole quick lane (5 of 7 runs, none before; no
    behavior change found): Prism's axe test of the trajectory bar now waits for its condition, the cause open in
    `packages/prism/TODO.md` (Decisions open).
  - **Flaky Dawn tests**, older than this batch: `webgpu-bundle-frag-depth.test.ts` ("draws writing frag_depth
    replay ...") fails about one run in four on Windows (replays 20 against 18): it compares replay counts across three
    plugin sessions, which depend on when pipelines finish compiling. `webgpu-points.test.ts` ("picks a ring point on its
    ring") and `webgpu-sphere-resolve.test.ts` ("shades alike from near and far") fail now and then in the whole
    suite, passing alone.
  - **W2 done**: the draw lists' context members one optional `GraphicsContext.drawLists` (`DrawLists`: `create`,
    `submitCount`, `resourcesVersion`, `capture`), there with all four or not at all; the renderer checks it once, not
    `captureDraws` apart. `RenderItem.renderFullscreen` and `getDataLayer` stay apart: features of their own (the
    sphere raster's resolve, kernels reading an item's data), not a set.
  - **W3 done** as agreed: `DrawMode` without `line-loop` and `triangle-fan` (no user; WebGPU has neither, so its
    mapping no longer throws), the context's doc comment above its own declaration. The rest of the review's W3 list
    left as is: `doNotForceWebGLContextLoss` (a WebGL option of the neutral destroy, named so), `readTextureData`'s
    row rule (documented; rendered content read by `readTextureAsync`), the XR projection in GL clip depth (as the
    camera builds every projection, see `depth-conventions.md`), the WebGL1 schema fields (with WebGL1, W4).
  - **S5 done**: the guards run on the shared harness (`test-harness/`, `@molstar/test-harness`, moved out of
    Seq\*), the graphics `test:perf` tier (`*.quick.ts` under `src/`, baseline `src/_test/perf-baseline.json`):
    `gl/_test/frame.quick.ts`, a frame at rest (the view unchanged) of 8 objects, its update and its opaque pass, on
    WebGL (the GL shim: WebGL1, every uniform classic) and on WebGPU (a stub device, `gl/_test/stub-gpu.ts`: the draw
    list replayed), steady: update 96 and 112 bytes a frame, opaque pass 128 and 784 (3,704 before WebGPU's: a pass state
    per pass, now one rebound; the global bindings' state a key string per replay, now a reused object and a number per
    texture set; the dynamic state a key string, now compared; a pass descriptor per pass, now one per target). Into a
    render target (an offscreen pass) 128 and 784 (2,320 before: its pass target made per pass, now kept per
    framebuffer, made again where its attachments, their views or its size change). Of the 784, about 470 are Node's
    `queueMicrotask` (an `AsyncResource` a frame; the context's submit at the task's end). On the way: the classic
    uniforms' comparison boxed every float it read (megamorphic, 34 KB a frame; compared as `Float64Array`s now, P9's
    `.slice()` gone), the global bindings' missing textures allocated per program change, and the update and each pass
    boxed their matrices and vectors anew (`ValueCell.update`, 744 bytes a frame; set only where changed now, no
    version bump or uniform block write at rest). The hover (`gl/_test/hover.quick.ts`: a highlight moved between
    spheres and a Canvas3D tick, the marking redrawn): 10.8 KB on WebGL, 17.2 KB on WebGPU, spread over the marking
    path. Candidates, each 0.2 to 0.9 KB a hover: the marker cells boxed anew (`ValueCell.updateIfChanged` in the
    visual's `mark`), the sort of the texture changes' ranges (`mergeRanges`), the marking layer's values per sample, the
    shader manager's `getRequired` per mark and commit, the shape's group iteration; on WebGPU `getPipelineFor` and
    `bindProgram` per draw of the marking passes.
  - **P5 decided, kept**: antialiasing stays on every multi-sample (SMAA, FXAA). Measured against 32 samples without
    antialiasing (a temporal level 5 reference, mean absolute difference): SMAA on every sample at level 2 0.77, on the
    unjittered one only 0.46 (WebGL2, WebGL1); FXAA (Dawn) 2.17 against 1.25. But at the default level 2 the edges of
    the four samples' average alone look feathery: the slight blur of antialiasing each sample is preferred, which the
    mean difference to a 32 sample average does not capture. Its cost: SMAA 0.10 ms GPU a sample on WebGPU, 0.18 on
    WebGL2 at 2560×1440 (FXAA 0.06–0.08). `renderer/multisample/smaa` and `renderer/marking/multi-sample-smaa` added
    (the harness runs without antialiasing). SMAA's image in headless Dawn is wrong (mean difference 230, with or
    without multi-sampling; a guess: its lookup textures not loaded there): not looked into.
  - **P2, lazy part done**: a moving frame of temporal multi-sampling to the drawing buffer presents the draw pass's
    color target; the compose target is written from it when sample 1 follows (nothing else draws into it between).
    Measured (`passes/temporal-moving`, 2560×1440, spacefill, RTX 4070 Ti): 0.55 → 0.49 ms GPU a moving frame on
    WebGL2 and WebGPU, as without multi-sampling. Not measured on an integrated GPU: Chrome on Windows kept the
    discrete one (`--force_low_power_gpu`, `--use-adapter-luid` ignored). Images: the moving frame as mode off's
    (`webgpu-temporal-multi-sample.test.ts`), settled as before on every WebGL profile (images saved and compared),
    `renderer/multisample/canvas-*` new (the canvas captured; the harness now flushes before reading it: WebGPU submits
    a task's commands in a microtask, a multi-sampled frame read blank).
  - **P2, running average done**: temporal multi-sampling's hold replaced by a running average in the compose
    target (constant blend `1/(k+1)`, as tracing and the marking layer). Settled: as hold's (at most 1 level apart;
    against mode `on` on Dawn, level 2: max 4 either; level 5: max 27 against 26, mean 0.34 against 0.36, `on`'s
    uneven weights the distance there). The 8 bit compose (WebGL profiles without float targets) better: against `on`,
    level 5 mean 0.63 against 1.99, no banding. Intermediate: closer to the settled image at level 4 (mean 0.43 against
    0.72 after the first frame), farther at level 2 (0.56 against 0.42: hold weighted the unjittered frame more).
    Speed: 0.03 to 0.04 ms GPU a frame at 2560×1440. Memory: the hold target, 4 B/px (14.7 MB at 2560×1440, 33 MB at
    3840×2160), in the canvas's pass and the image pass's each. Its pipeline planned by `prepare` (a test: no draw
    skipped settling). Illumination's multi-sample compose the same (its iterations the samples; settled images unchanged on WebGL2, WebGL1 and WebGPU; 4 B/px less with multi-sampling on).
  - **P1 done**: the first frame after a change (unjittered) is sample 0, every level's first offset being none: its
    scene not drawn again (it added nothing to the float compose, the hold's own share to the 8 bit one, seeded with
    it now). A settle draws the scene `n` times, not `n + 1`: at the default level 4, not 5, a frame sooner. The
    jittered samples reuse that frame's occlusion (`reuseOcclusion`), as they did sample 0's. Test
    (`webgpu-temporal-multi-sample.test.ts`): the frames past those of mode `on`, 3 at level 2 (4 before), and the
    settled image as `on`'s (as before: 26 of 9216 values more than 2 apart, by 4 at most, the accumulations' rounding).
  - **P3 done**: measured, a moving view with Hi-Z and bloom made 36 bind groups a frame (about 0.35 ms of CPU on
    native Dawn, more in a browser: each a call to the GPU process): Hi-Z's levels (16, a frame of every scene with
    objects culled on the GPU, `hiZ.gpu`) and bloom's blur (20, scenes with emissive objects), their inputs swapped
    every frame. A render item keeps up to 8 groups a layout for the textures and buffers they bind
    (`alternateGroups`), found again on a swap; its programs invalidated only where the linear-filtered textures
    changed. None a frame then (`webgpu-pass-bind-groups.test.ts`; 34 a frame without).
  - **V1 done**, two unrelated causes:
    - The WebGL2 stall was the harness's load, not the case: debug mode's synchronous GL queries (each waits for the
      GPU process, the commands of every page) made cases with many programs and draws wait over a minute. The compile
      and link checks now run once linked (`finalize`, after the parallel compile's completion), the GL error check
      per pass, not per draw. The default run: no stalls (2 to 9 a run before), 74 s, not 133 s. Pages in a GPU
      surface case still answer the harness's poll only after up to 29 s (readbacks behind the other pages' work, not
      hung): its hang limit 90 s, not 30 (two more cases moved one past 30).
    - The WebGPU failure: a frame drawn straight to the drawing buffer leaves the opaque depth texture as an earlier
      frame drew it (an earlier case's, offscreen for its overlay), and Hi-Z built its levels from it; an image pass
      culled with them (no second phase) and lost instances (`merged-clip-pixel`, `merged-clip-instance`,
      `merged-grid-lod` after the same cases). The levels now only of a frame that drew the texture, cleared
      otherwise (`DrawPass.opaqueDepthDrawn`; `webgpu-hi-z.test.ts`).
  - **S4 done**: `gl/_test/shader-parity.test.ts`, per graphics program the schema names its GLSL and its WGSL read
    (comments not; the WGSL with every feature of its variants on), the differences listed and why (WebGL mechanics,
    grids only on meshes, the sphere raster's uniforms); a new one-sided name fails, and so does a listed one gone. A
    name check, not a translator: a feature on one backend only shows, a use dropped from one file of a program does
    not.
  - **Batch B, CI** (S1): a `WebGPU tests (Dawn)` job in `node.yml` (Mesa's lavapipe Vulkan, the headless WebGPU
    files, `MOLSTAR_REQUIRE_WEBGPU`: no Dawn or no adapter fails instead of skipping) and `visual-tests.yml` (the
    smoke set nightly, by hand and on a PR labelled `visual tests`, `--require-backends`: a profile without its backend
    fails). Measured here (Windows, Chromium): `--angle swiftshader` leaves Chromium without WebGPU (every WebGPU case
    skipped, the run passing); `--use-webgpu-adapter=swiftshader` gives a SwiftShader WebGPU adapter beside WebGL on the
    GPU. Whether the runner gives both on Mesa (ANGLE on Vulkan, Dawn on Vulkan) is for the first run to show; the
    flags may need a turn.
  - **Found on SwiftShader WebGPU, fixed**: every GPU marching cubes surface shredded (11 cases, 2.8k-15.8k px; the
    smoke set's `isosurface-gpu` 4973 px). A probe reading back the extraction's buffers: counts, pages, vertices and
    the workgroups' triangle bases as on a GPU, but 3386 of 4200 triangles never written. The kernels added to
    workgroup counters (`wgVertices`, `wgTriangles`, `wgOffset`) relying on WGSL's zero-initialized workgroup memory,
    which SwiftShader's adapter skips (the previous dispatch's values): triangles past their limit, dropped. Reset by
    the kernels themselves now (as `cell-cull` and `sphere-raster` already did); a static test that every atomic
    workgroup counter is stored first (`gl/_test/wgsl-workgroup-init.test.ts`). The 16 GPU surface cases pass on
    SwiftShader and on the GPU.
  - **The first CI runs**: the Dawn job ran (one test past the quick lane's 30 s on the CPU: the job's timeout raised);
    the visual run had both backends on Mesa's Vulkan, 17 of 19 WebGPU cases matching WebGL2: the isosurface above,
    and `geometry/text/labels` (141 px of anti-aliasing and a box edge between the two software rasterizers: a
    tolerance of its own). Its third label was mojibake since the harness began (`Ã…Ã„Ã–` for `ÅÄÖ`), fixed.
  - **The Dawn job canceled at about 2 min** (out of memory: 15.9 of 16 GB, rising with no jump): every destroyed
    WebGPU context stayed reachable, its device's `uncapturederror` listener (whose scope is the whole context) kept by
    Dawn's Node.js bindings, which hold the device by a handle of their own also when destroyed; a run's contexts and
    their resources (in RAM on lavapipe) held to its end. The context removes its listener on destroy; a heap snapshot
    gave the retainer path, `webgpu-device-listener.test.ts` the regression test.
  - **Batch A done**: C1 (WebGL uploads return once destroyed, the cube's faces counted and not uploaded), C2 (a
    skybox face whose source fails reports an error, `onload` once, its object URLs revoked), C3 (a background load
    marks its own entry, also one reported before the entry exists), C4 (the buffer pages reclaimed after each
    submission and once idle after a free, `BufferPages.reclaim`), C5 (WebGL material buffers released with their last
    user), C6 (lines past their index limit not drawn, as boxes), C7 (no WebGL1 mipmaps of NPOT cube faces), C8 (the
    grown index buffer made before the old one is released), C9 (`Passes.dispose`: a canvas releases its old passes'
    background assets once the new ones drew). Tests: `gl/_test/webgl-late-loads.test.ts`, the buffer pages' in
    `gl/_test/webgpu-buffer-pages.test.ts`; C3, C5-C9 without one of their own.

The review (`molstar-graphics-review.md`, outside the repository) is a static reading of `packages/graphics` at a commit
not in this repository, most likely `master` before the gizmo, clip ghost and marking work. Its line numbers drift;
each item below was checked against `gl-gizmos-clip-ghosts-mesoscale`. It did not read the merged objects, the culler,
the sphere raster or `scene.ts`.

## Triage

| Id | Verdict on this branch | Batch |
| --- | --- | --- |
| C1 | Confirmed, fixed: WebGL `load`, `loadData`, `loadRows` and the cube `load` have no `destroyed` check (`gl/webgl/texture.ts`); WebGPU's has it | A |
| C2 | Confirmed, fixed: a skybox face's `source.then` has no rejection handler, so `onload` never runs and `pending.whenReady()` hangs; its blob URLs are never revoked | A |
| C3 | Confirmed, fixed: a superseded background load writes `this.skybox.loaded` (the current entry), not its own (`canvas3d/passes/background.ts`) | A |
| C4 | Confirmed, fixed: WebGPU `BufferPages.free` only queues; pages are reclaimed by a later `allocate` on the same pages | A |
| C5 | Confirmed, fixed: the WebGL material uniform cache evicted only at context destroy | A |
| C6 | Confirmed, fixed: lines past the shared index limit (only boxes guarded) | A |
| C7 | Confirmed, fixed: WebGL1 NPOT skybox faces with mipmaps (an incomplete cube map) | A |
| C8 | Confirmed, fixed: `QuadIndices.get` releases the old buffer before the new one is made | A |
| C9 | Confirmed, fixed: `BackgroundPass.dispose` without callers, `Passes` without `dispose` | A |
| S1 | Confirmed, fixed: CI installs the native backend without `--webgpu`, so the Dawn suite skips; visual tests by hand only | B |
| S4 | Done: `gl/_test/shader-parity.test.ts` | B |
| P1 | Done: the first frame after a change is sample 0, not drawn again | C |
| P3 | Done: bind groups kept per bound textures, programs only on a filter change | C |
| V1 | Found here, older: two visual cases failing only after others in a page (`merged-clip-pixel` on WebGPU, `isosurface-gpu-negative` hanging on WebGL2) | done |
| P2 | Done: moving frames compose nothing, hold replaced by a running average (see Status) | C |
| P4 | Done on this branch: marking reads the scene's depth, per jittered sample too (marking-performance.md step 3) | - |
| P5 | Decided: kept on every sample (feathery edges at level 2 without; see Status) | C |
| P6 | Done: a pick read in one staging buffer and command buffer, the buffer kept, rows copied once | D |
| P7 | Decided, kept: hover picks while the camera animates (see Status) | - |
| P8 | Done: WebGL skips attaching what an attachment point holds already (see Status) | D |
| P9 | Done: no closures per direct WebGPU draw, the mark buffer emptied in place (see Status) | D |
| S2, S3 | Done as far as they cut: `Renderer.create` and `Canvas3D.create` split along their seams (see Status) | D |
| S5 | Done: a frame at rest and the hover guarded on WebGL and WebGPU (see Status) | D |
| S6 | Done: the backend folders internal, contexts through `gl/contexts/*`, extensions through `gl/extend/*` | D |
| S7 | Done: graphics' folders import one way, held by `src/_test/layers.test.ts` (see Status) | D |
| S8 | Done: the named caches a scratch per module (`gl/scratch.ts`, see Status); kept for the context's lifetime, no idle release | D |
| W1 | Done: render pass attachments, framebuffers internal to the backends (see Status) | D |
| W2, W3 | Done: the draw-list members grouped (`GraphicsContext.drawLists`), the unused GL draw modes removed | D |
| W4 | Decided: WebGL1 stays for now (as `minimum-capabilities.md` keeps it) | - |

## Steps

1. **Batch A, correctness** (small, one change each, tests where a stub reaches): C1-C3 as one pattern, a callback
   outliving what it was started for: a generation or `isCurrent` check, with one test destroying a texture before its
   decode resolves and one failing a skybox face; then C4 (reclaim after a submit, an idle canvas's follow-up), C8, and
   C5-C7, C9 once checked.
2. **Batch B, verification**: a Dawn job in CI's medium lane (`native:install --webgpu`, a software adapter; whether
   Dawn's prebuilt binaries run on the runner is the first question) and a `--tags smoke` visual run nightly; the GLSL
   and WGSL schema-name cross-check as a quick test.
3. **Batch C, performance, measured first** (`--gpu-timing`, the benchmark cases): P1, P3, P2; P5 by image comparison.
   Each kept only where it measures. P2's larger variant (hold replaced by a running average, blend `1/(k+1)`, as
   tracing and the marking layer do) measured for all three: speed (the frame time per sample: hold's second quad),
   memory (the hold target's 4 B/px and the targets left, at the canvas sizes measured) and quality (the settled image
   against hold's, and the intermediate ones, by the visual tests' comparison).
4. **Batch D, structure**: S5 first (guards make the splits safe), then S2 and S3 by their seams, W2 and W3, S6 before
   v6 freezes; W1.
5. **The backend folders internal (S6)**: contexts made through helpers outside them, `gl/contexts/webgl.ts` (a context
   from a GL context, a GL context from a canvas: `createContext`, `getGLContext`) and `gl/contexts/webgpu.ts`
   (`requestWebGPUDevice`, a context from a device, `WebGPUContextOptions`), each re-exporting its backend descriptors;
   `canvas3d/webgl.ts` and `plugin/headless/src/screenshot.ts` import them. Not in `gl/backends/*`: those are the
   backend sets a single-backend build aliases to the kept one (`scripts/esbuild/backend-sets.mjs`), so a WebGPU helper
   there would be the WebGL module in a WebGL build. Then nothing outside `graphics` imports `gl/webgl/**` or
   `gl/webgpu/**` (inside, `backend-neutral.test.ts` already holds it). The repository's `tests/` harness reads backend
   internals (compute kernels, programs, constants): exempt as the in-repository harness, or what it needs behind
   `gl/contexts/*` too, then the exports narrowed. To choose when this is taken up.

## Decisions

- WebGL1 stays for now (W4).
- The backend folders are internal, contexts made through `gl/contexts/*` (S6, step 5).
- P2's larger variant taken: hold replaced by a running average (measured, see Status; `pass-targets.md` updated).
