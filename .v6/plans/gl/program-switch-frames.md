# Plan: Frames while programs change (first use of a feature)

## Status

- **Current**: A done (step 2); B done for props changes through
  [pipeline-readiness.md](pipeline-readiness.md) (`holdFramesUntilPrepared`: `Canvas3D.setProps` prepares the props
  from params, and while that work is pending frames are skipped, the last frame stays shown; the Viewer smoke test
  reads back every frame after turning on illumination: none of only the background, 13 of 50 before), for at most
  1 s (`MaxFrameHoldMs`), so work that never settles doesn't freeze the canvas. Step 4 decided: new representations and
  define changes keep presenting frames without their draws until the pipelines compiled (they appear when ready,
  which feels right in use).
- A: `uSteps`, `uRefineSteps` and `uRendersPerFrame` are loop bound uniforms
  (defines on WebGL1, where the tracing pass updates the renderable when they change). Turning on illumination now
  compiles 11 pipelines and 9 modules on WebGPU (1 `trace`; was 51 and 49, 40 `trace`), the last done after 315 ms (was
  1.4 s); WebGL links 6 programs, the last after 360 ms (was 12, 1.9 s). Images unchanged (WebGL2 and WebGL1
  byte-identical, WebGPU against WebGL2 as before). The trace pass on WebGL2 takes longer with uniform bounds: 0.16 →
  0.24 ms at the first iteration's quality (`render/7qpd-illumination-rotating`; half of it the refine loop, which ANGLE
  no longer unrolls), WebGPU unchanged (0.13 ms); kept, as adapting the quality no longer stalls WebGL with synchronous
  links. Still open: the almost empty WebGPU frames for about 150–250 ms after the switch (B), and WebGL's one
  synchronous link of the illumination programs (about 900 ms without a new frame). Follows
  [resource-readiness.md](resource-readiness.md) (pending work, background pipeline compiles, collect mode
  `preparePrograms`).

## Problem

Turning on a feature for the first time (e.g. illumination) on WebGPU shows an almost empty canvas for a few frames: the
draws whose pipelines are still compiling are skipped, but the frame is presented. WebGL keeps the previous image
instead (`Canvas3D.draw` skips the whole frame while the required programs aren't linked), at the cost of a frozen UI
while it links synchronously.

Measured in the Viewer (1cbs, 800 × 600, Chrome, NVIDIA), turning on illumination
(`canvas3d.setProps({ illumination: { enabled: true } })`), canvas screenshots in a loop:

|                               | WebGPU                                                                                                                                | WebGL                                                                           |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| frames until a complete image | about 150–250 ms almost empty (0.2–0.7% coverage of 7.6%), then nearly final                                                          | the previous image stays, no new frame for about 800 ms                         |
| compiles                      | 51 pipelines and 49 shader modules in 1.4 s, the slowest 300–400 ms (geometry `tracing` variants into the three illumination targets) | 12 program links in 1.9 s, each synchronous; frames 200–800 ms apart until then |
| `trace` program               | 40 modules and pipelines                                                                                                              | part of the 12 links                                                            |

Two causes:

1. **Skipped draws are presented (WebGPU).** The first frame with illumination needs about 11 new pipelines: the
   `tracing` variants of the geometry programs, `trace`, `illumination-compose`, `copy`. Their draws
   are skipped while they compile in the background, and the frame is shown without them.
2. **Adaptive tracing recompiles (both backends).** `TracingPass.getAdjustedProps` adapts `steps`, `refineSteps` and
   `rendersPerFrame` per frame to reach `targetFps`, and all three are defines (`dSteps`, `dRefineSteps`,
   `dRendersPerFrame`; WGSL override constants on WebGPU). Every new combination is a new program: 40 `trace` programs
   in the first 1.5 s on WebGPU, synchronous links on WebGL that stall frames while it converges.

## Goal

The first use of a feature shows the previous image until the new one is complete, without freezing the UI: on WebGPU
the pipelines compile in the background while the previous image stays; adapting the tracing quality compiles nothing.

## Design

### A. Adaptive tracing without recompiles

`steps`, `refineSteps` and `rendersPerFrame` become uniforms marked `loopBound` (the mechanism of
[backend-narrowing.md](backend-narrowing.md) step 1): WebGL2 and WebGPU loop over uniform bounds, WebGL1 turns them into
defines as today (its GLSL 100 loops need constant bounds). One `trace` program per configuration of the remaining
defines; the adaptive loop changes uniform values only. `dBounces`, changed by the user rather than per frame, stays a
define.

- Risk: loops with constant bounds unroll and may be faster; measure the trace pass GPU time (timing mode, `render/*`
  benchmark with illumination) with constants and with uniforms. If a bound matters, keep a define for the maximum and
  loop to the uniform inside it.

### B. Hold the frame until its pipelines exist (WebGPU)

Like WebGL's `finalizeRequired`, but asynchronous: a frame that would skip draws isn't presented; the canvas keeps the
last presented image (a WebGPU canvas presents only when the frame takes its current texture), and the frame is redrawn
when the pipelines compiled (`pending.settled` requests it).

- **Detect before drawing**: when the frame may need new pipelines (programs created since the last complete frame, e.g.
  a define change or new representation; render props changed, e.g. illumination, transparency mode, postprocessing),
  `Canvas3D.draw` first runs the frame in collect mode. A synchronous variant of `preparePrograms` returns whether the
  run requested anything missing (`resources.preparePrograms` resolves later; the new call answers now, e.g.
  `resources.hasPrograms(run): boolean`). If something is missing, the frame is skipped and `forceNextRender` set, as on
  WebGL.
- **Cost**: the collect run is a frame's CPU work without encoding, only on frames after such changes. A context counter
  of created programs and requested pipeline configurations tells `Canvas3D` when a frame needs the check, so steady
  frames don't pay for it.
- **Fallback**: draws skipped in a presented frame anyway (a pipeline the collect run didn't predict) are counted; the
  frame requests a redraw as today, and the count is reported in debug mode so missing predictions surface.
- **Picking and marking** keep their own paths (picking renders offscreen, the marking warm-up prepares it).
- **`ImagePass`** keeps rendering again while a render started work (it renders asynchronously, see resource-readiness).

## Steps

1. **Measure**: the script used above (a Viewer page, a prop switch, screenshots per frame against the final image,
   pipeline and module counts) as a check in the Viewer smoke test or a separate script: frames after turning on
   illumination never fall below the previous image's coverage, and the number of `trace` programs.
2. **A**: tracing loop bounds as `loopBound` uniforms; trace GPU time before and after; visual suite (WebGL2 and WebGPU
   images unchanged within the illumination tolerance, WebGL1 byte-identical); the switch measurement shows one `trace`
   program per configuration and no WebGL links after the first frames.
3. **B**: the synchronous readiness check in the WebGPU context, the program and configuration counter, `Canvas3D.draw`
   holding frames; check that a WebGPU canvas keeps its image across skipped frames (Chrome, Firefox); the switch
   measurement shows the previous image until the complete one; visual suite unchanged.
4. **Other first uses**: the same measurement for transparency modes (wboit, dpoit), postprocessing (outlines, SSAO,
   bloom, DOF), a new representation and a theme change (define changes), to find what B doesn't cover.

## Open questions

- Whether holding the frame should have a limit: yes, 1 s (`MaxFrameHoldMs`).
- Whether the collect run should also prepare the variants a frame will need next (e.g. the illumination accumulation's
  later passes) or rely on holding each frame until it is complete.
