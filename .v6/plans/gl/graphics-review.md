# Plan: follow-ups from the graphics layer review

## Status

- **Current**: batch A done; next: batch B (WebGPU in CI, the GLSL/WGSL cross-check).
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
| S4 | Agreed: no check that a program's GLSL and WGSL read the same schema names | B |
| P1 | Confirmed that every jitter level starts at `[0, 0]`; whether the first sample after a stop is wasted to be measured | C |
| P3 | To measure: WebGPU quad renderables invalidating programs and bind groups on texture swaps (Hi-Z, bloom) | C |
| P2 | To measure: hold and compose written on every moving frame; its larger variant for quality, speed and memory | C |
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
- Open, after measuring quality, speed and memory (step 3): P2's larger variant, hold replaced by a running average
  (changes intermediate images; `pass-targets.md` kept hold without a recorded cost).
