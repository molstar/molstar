# Plan: illumination performance on mesoscale models

## Status

- **Current**: done. Option 2 done for the tracing input (WebGPU), option 3 measured (a fixed thickness of 256 Å looks as
  auto on the mesoscale models; auto stays the default, too thick on other models); option 1 shelved, options 4 and 5 dropped. Benchmark cases `illumination-mesoscale/{postsyn,mg}-{whole,zoom0.3}`
  with `-large` (2560x1440), `-still` and thickness variants (`-fixed` 4 Å, `-fixed16`, `-fixed64`, `-fixed256`;
  PetWorld and MG_1184 as the viewer draws them, Hi-Z on, the mesoscale explorer's `rayDistance: 1024`; on demand, the
  models in `.benchmark/`, `--gpu-timing`).
  - **Option 2 done, the tracing input only**: the tracing pass culls on the GPU before its input
    (`Renderer.prepareCull`, the previous frame's pyramid), and after the input's first part runs the second phase
    against the pyramid of its depth so far (`HiZPass.planFrame`, `prepareCullSecond`, `renderTracingSecond`), as the
    draw pass's opaque pass. `Canvas3D` builds the levels after the inputs, not read back. Only with `hiZ.gpu` and
    `twoPhase`: the image accumulates the inputs, so what an earlier frame's pyramid hides and the view reveals must be
    drawn in the same frame; no CPU test (it lags behind, and WebGL, with no GPU culling, stays without occlusion).
  - **The back depth is never culled by the front's occlusion** (its pyramid or lists): it is the farthest back faces
    of everything in view, also of what the front hides whole; culled with the front's lists it lost those (the
    thickness of the visible layer only, another shading). It is drawn under a frame of its own, so the merged objects
    draw the CPU culling's lists (frustum and levels of detail, no occlusion test in illumination): its draw calls and
    CPU submit as before. A frustum-only GPU cull of its own would need a second set of the culler's lists.
  - Tried and dropped: a CPU occlusion test for the inputs with a redraw without it at rest (zoomed into MG 67.7 ->
    40.0 ms of inputs on WebGL2, but the whole model's CPU submit 9.7 -> 21.7 ms: a test per cell, split runs).
  - Test: `webgpu-illumination-occlusion.test.ts` (a wall hiding a row of copies, a jump to the side: the frame, the
    back depth and the iterations after as without occlusion; the back depth of the front view holding the hidden
    copies).
  - **Measured**, WebGPU, rotating, GPU ms per frame (CPU submit): before, after option 2 (auto thickness, the back
    depth unculled), and with `thicknessMode: 'fixed'` (no back depth):

    | case                  | before      | after       | fixed       |
    | --------------------- | ----------- | ----------- | ----------- |
    | postsyn whole         | 2.60 (0.77) | 1.81 (1.38) | 0.76 (0.97) |
    | postsyn whole 2560    | 5.11 (0.72) | 3.64 (1.15) | 1.89 (0.75) |
    | postsyn zoom 0.3      | 8.83 (5.83) | 5.70 (6.18) | 2.40 (0.91) |
    | postsyn zoom 0.3 2560 | 34.4 (6.49) | 22.8 (6.45) | 10.3 (0.92) |
    | mg whole              | 5.61 (10.3) | 4.07 (9.91) | 1.63 (1.20) |
    | mg whole 2560         | 14.8 (9.83) | 11.8 (11.0) | 6.21 (1.18) |
    | mg zoom 0.3           | 17.1 (8.84) | 12.9 (8.10) | 6.27 (1.16) |
    | mg zoom 0.3 2560      | 72.7 (10.2) | 51.0 (7.35) | 23.8 (1.18) |

    The tracing input 2-4 times cheaper (MG zoomed at 2560x1440 42.0 -> 17.0 ms); the culling kernels and the pyramid
    add 0.1-0.4 ms. Part of the gain is GPU culling itself (the input drew the CPU culling's lists before: no levels of
    detail by the GPU, a draw call per segment). Images against before: scattered pixels (PetWorld zoomed 2.9% of
    pixels off by more than 8 levels with the back depth culled, the rest under 0.4%), nothing missing.
  - **Option 3 measured**: still frames (about 120 iterations), against auto, mean difference in levels: fixed at 4 Å
    21-36 zoomed (flat, light, creases gone, light through the gaps), 16 Å 15-22, 64 Å 8-10, 256 Å 0.14 (MG zoomed),
    1.8 (PetWorld zoomed: alike by eye), 0.14-0.19 whole. So on these models auto amounts to an effectively unbounded
    thickness: the back depth is the back of the cell. With a fixed thickness the trace is also 15-20% cheaper (no
    back depth reads: MG zoomed 14.6 -> 12.2 ms, PetWorld zoomed 11.8 -> 9.6 ms a still frame).

## Where things stand

- **Passes**: illumination draws two inputs, then the trace (`TracingPass`):
  - `tracing input`: the opaque scene into color, normal and shaded targets with the opaque depth.
  - `thickness` (`thicknessMode: 'auto'`, the default): the farthest back faces of the opaque scene
    (`renderDepthOpaqueBack`, depth reversed) into a packed depth target at full resolution. The trace reads it only
    where a ray step lands behind the front depth: the thickness accepted for a hit is
    `max(minThickness, (front z - back z) * thicknessFactor * density)` (`trace.frag.ts`, `illumination.wgsl.ts`).
  - `trace`: per pixel `rendersPerFrame` rays of up to `bounces` bounces, each a march of `steps` depth reads and
    `refineSteps` refinements. The quality adapts per frame to `targetFps`: the first iteration (`iteration === 0`,
    every frame while the camera moves) at its minimum (half the rays, half the steps, one refinement), later
    iterations climbing to the props' maximum.
  - The inputs are drawn on the first iteration only (or a new multi-sample offset); a still camera iterates the trace
    alone.
- **Occlusion culling**: on WebGPU the tracing input is culled on the GPU in two phases (option 2), the back depth not
  by occlusion; on WebGL neither (no GPU culling, no CPU test: it lags behind and the image accumulates what it
  misses).
- **Measured before option 2** (RTX 4070 Ti, WebGL2 through ANGLE's Vulkan backend; GPU ms per frame; rotating: the
  first iteration with its inputs, still: the iterations at the adapted maximum, the inputs drawn once):

  | case (rotating unless still)  | frame | tracing input | thickness | trace | compose | draw calls |
  | ----------------------------- | ----- | ------------- | --------- | ----- | ------- | ---------- |
  | postsyn whole                 | 2.60  | 1.22          | 1.06      | 0.12  | 0.14    | 308        |
  | postsyn whole, still          | 2.28  | -             | -         | 1.96  | 0.14    |            |
  | postsyn whole 2560            | 5.07  | 2.48          | 1.67      | 0.43  | 0.35    | 308        |
  | postsyn whole 2560, still     | 6.10  | -             | -         | 5.53  | 0.31    |            |
  | postsyn zoom 0.3              | 8.91  | 4.82          | 3.10      | 0.37  | 0.55    | 14k        |
  | postsyn zoom 0.3, still       | 10.5  | -             | -         | 9.68  | 0.50    |            |
  | postsyn zoom 0.3 2560         | 33.8  | 18.1          | 11.8      | 1.50  | 2.13    | 22k        |
  | postsyn zoom 0.3 2560, still  | 41.3  | -             | -         | 37.8  | 2.03    |            |
  | mg whole                      | 6.09  | 2.95          | 2.55      | 0.21  | 0.31    | 21k        |
  | mg whole, still               | 5.73  | -             | -         | 5.08  | 0.32    |            |
  | mg whole 2560                 | 15.1  | 7.45          | 5.62      | 0.81  | 1.05    | 26k        |
  | mg whole 2560, still          | 21.7  | -             | -         | 20.0  | 1.11    |            |
  | mg zoom 0.3                   | 17.1  | 9.17          | 6.78      | 0.43  | 0.64    | 30k        |
  | mg zoom 0.3, still            | 13.4  | -             | -         | 12.1  | 0.66    |            |
  | mg zoom 0.3 2560              | 71.8  | 41.1          | 26.6      | 1.39  | 2.51    | 26k        |
  | mg zoom 0.3 2560, still       | 49.9  | -             | -         | 45.0  | 2.73    |            |

  Without illumination, the same views still (`sphere-raster-zoom/postsyn-viewer-zoom0.3-off`,
  `sphere-raster-mg/zoom0.3{,-large}-off`, the impostors on WebGL2, Hi-Z on), the opaque pass: postsyn zoomed 4.38 ms
  (11.5k draw calls), MG zoomed 7.47 ms (11.8k), MG zoomed at 2560x1440 21.0 ms (11.6k).

  - While the camera moves, the inputs are 85-95% of the illumination's GPU time and the back depth 39-47% of the
    inputs: a second full geometry pass, sphere impostors computing the far intersection and writing depth (no early
    depth test), every overlapping sphere drawn, nothing occluded culled. The trace at its minimum quality is small.
  - The tracing input alone costs about as much as the opaque pass with occlusion culling at 1280x720 (4.82 / 4.38 ms,
    9.17 / 7.47 ms) and twice as much zoomed into MG at 2560x1440 (41.1 / 21.0 ms); with the back depth, a moving frame
    with illumination costs 1.8-3.2 times the opaque pass. The draw calls double (the back depth draws every segment
    again): on WebGL the CPU submits 9-58 ms a moving frame zoomed (11-21 ms without illumination, CPU culling per
    segment), as much as the GPU's work or more.
  - Still, the trace is almost all of it, scaling with the pixels the scene covers: 2 ms (PetWorld whole, 1280x720) to
    45 ms a frame (MG zoomed, 2560x1440) at 16 rays, 32 steps, 4 refinements. The benchmark's frames are submitted
    without waiting for the GPU, so the adaptation saw short CPU intervals and climbed to the maximum; what a browser's
    frame pacing makes of it is to check.
  - Over a remote desktop session Chrome's WebGPU adapter is WARP (software); WebGL2 reaches the GPU there with
    `--browser-arg --use-angle=vulkan`. The WebGPU numbers before option 2 (its table) match WebGL2's within about 10%.
    The sphere raster does not draw the inputs (the tracing and depth variants are impostors).

## Options

1. **The back depth at half resolution** (shelved): the thickness target at half the size per axis, read by the trace at the same
   normalized coordinates. Fragment work of the thickness pass divided by about 4; its vertex work and draw calls
   unchanged (the gain where it is fragment-bound: zoomed, large drawing buffers). The texels are coarse where objects
   meet: a back depth from a neighbour paired with a front depth gives too thin or too thick a thickness. A depth-aware
   pick in the trace (the four texels around, the nearest back depth behind the pixel's front depth; none behind:
   `minThickness`), as the SSAO upsampling (`getSsaoUpsampled`) weighs texels by depth. Expected: up to 3/4 of the
   thickness pass's GPU time where it is fragment-bound (zoomed at 2560x1440: 11.8 and 26.6 ms), less whole; the CPU
   side (its draw calls) stays. Images against full resolution to judge.
2. **Occlusion culling for the inputs** (done for the tracing input, WebGPU, see Status): the pyramid built from the
   tracing input's depth, two-phase GPU culling. Not for the back depth (it is the back of what the front hides).
3. **Is the auto thickness worth its pass on mesoscale models?** The farthest back face along a pixel of a dense model is
   far behind the front surface, so `front z - back z` is large almost everywhere and the hit test accepts nearly every
   step behind the depth. Compare `thicknessMode: 'fixed'` (no thickness pass) and auto by image on PetWorld and
   MG_1184; if alike, the mesoscale explorer (or the mesoscale preset) sets `fixed` and the pass is gone there, its
   GPU time (39-47% of the inputs) and its draw calls (half of the moving frame's).
4. **The back face of the front primitive** (dropped: of limited use in practice): the tracing input writes the far intersection of the primitive it draws
   (spheres and cylinders know it analytically) into a spare channel, a per-primitive thickness at no extra pass.
   Thinner than the object's (an atom, not a protein): different shading; for scenes that are mostly spheres.
5. **The trace at lower resolution for still frames** (dropped: the still frames' cost is acceptable as it is): the trace at half resolution with a depth-aware upsample, or
   fewer rays per frame where the drawing buffer is large; quality against the denoiser to judge. The larger cost when
   still (38-45 ms zoomed at 2560x1440), but it accumulates and stops after `2^maxIterations` frames.

**Decided**: auto stays the default, the mesoscale explorer's too: on other models 256 Å is too thick, and no fixed
value fits all. The numbers stand for a user setting `fixed` on dense models. Option 1 shelved, 4 and 5 dropped: nothing left to do.

## Steps

1. Option 2: done for the tracing input (WebGPU); the back depth unculled.
2. Option 3: measured (Status); auto stays the default everywhere.
3. Option 1: shelved (if taken up: the thickness target at half resolution in `TracingPass`, the depth-aware pick in
   the trace in GLSL and WGSL, measured with the cases above and compared by image with full resolution).

## Open questions

- WebGL: occlusion for the inputs would need a test without lag (none: no GPU culling there); left off.
