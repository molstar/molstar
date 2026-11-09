# Plan: follow-ups from the graphics layer review

## Status

- **Current**: batches A and B done, B verified in CI (the Dawn job and the visual smoke run passing on the runners);
  next: batch C (the multi-sample savings, measured first): P1, P2 and P3 done, next P5.
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
    skipped settling). Illumination's multi-sample compose keeps its hold (the same change possible there).
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
| C1 | Confirmed: WebGL `load`, `loadData`, `loadRows` and the cube `load` have no `destroyed` check (`gl/webgl/texture.ts`); WebGPU's has it | A |
| C2 | Confirmed: a skybox face's `source.then` has no rejection handler, so `onload` never runs and `pending.whenReady()` hangs; its blob URLs are never revoked | A |
| C3 | Confirmed: a superseded background load writes `this.skybox.loaded` (the current entry), not its own (`canvas3d/passes/background.ts`) | A |
| C4 | Confirmed: WebGPU `BufferPages.free` only queues; pages are reclaimed by a later `allocate` on the same pages | A |
| C5 | Confirmed, fixed: the WebGL material uniform cache evicted only at context destroy | A |
| C6 | Confirmed, fixed: lines past the shared index limit (only boxes guarded) | A |
| C7 | Confirmed, fixed: WebGL1 NPOT skybox faces with mipmaps (an incomplete cube map) | A |
| C8 | Confirmed: `QuadIndices.get` releases the old buffer before the new one is made | A |
| C9 | Confirmed, fixed: `BackgroundPass.dispose` without callers, `Passes` without `dispose` | A |
| S1 | Confirmed: CI installs the native backend without `--webgpu`, so the Dawn suite skips; visual tests by hand only | B |
| S4 | Done: `gl/_test/shader-parity.test.ts` | B |
| P1 | Done: the first frame after a change is sample 0, not drawn again | C |
| P3 | Done: bind groups kept per bound textures, programs only on a filter change | C |
| V1 | Found here, older: two visual cases failing only after others in a page (`merged-clip-pixel` on WebGPU, `isosurface-gpu-negative` hanging on WebGL2) | done |
| P2 | Done: moving frames compose nothing, hold replaced by a running average (see Status) | C |
| P4 | Done on this branch: marking reads the scene's depth, per jittered sample too (marking-performance.md step 3) | - |
| P5 | To decide by images: SMAA per jittered sample | C |
| P6-P9 | Small; with nearby work, P9 behind allocation guards | - |
| S2, S3 | Agreed, long-term: `Renderer.create` and `Canvas3D.create` split along their seams | D |
| S5 | Agreed: allocation guards on render items, the pass walk, draw-list replay | D |
| S6 | Decided: the backend folders (`gl/webgl`, `gl/webgpu`) are internal; contexts made through `gl/contexts/*` | D |
| S7, S8 | Agreed: layering inversions, owned kernel caches | D |
| W1 | Agreed, long-term: target descriptors instead of mutable attachments, pass by pass (also P8) | D |
| W2, W3 | Agreed: WebGPU-only draw-list members grouped, unused GL draw modes removed | D |
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
   v6 freezes; W1 pass by pass when a pass is touched anyway.
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
