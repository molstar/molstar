# Plan: SSAO performance

## Status

- **Current**: options 1 and 2 tried. Benchmark cases `render/assembly-60-spacefill-ssao` and
  `render/assembly-60-spacefill-ssao-multiscale` (60 copies of 7qpd as spacefill, 2560x1440 for stable GPU times).
  - **Option 1 done**: below the full resolution the composite (postprocessing, illumination compose) upsamples the
    occlusion depth-aware (`getSsaoUpsampled`: the four nearest texels, bilinear weights over the difference of their
    linear depth to the pixel's). Against the full-resolution image the pixels off by more than 16 levels drop from 233
    to 29 (`multi-scale-half-resolution`, its WebGPU tolerance gone); multi-scale SSAO at 0.5 costs 0.19 ms instead of
    0.71 ms (2560x1440, WebGL2; WebGPU 0.18 / 0.66), the upsampling nothing measurable. `resolutionScale` stays 1 by
    default, then made the default (0.5, the scale at most that of CSS pixels, `min(resolutionScale, 1 /
    pixelRatio)`); the visual cases cover full resolution explicitly (`multi-scale-full-resolution`,
    `renderer/marking/dim-occlusion`). At the default, 2560x1440, 16 instead of 32 samples: SSAO 0.09 -> 0.08 ms
    (single scale), 0.22 -> 0.13 ms (multi-scale, WebGL2; WebGPU 0.20 -> 0.12). The marking overlay (dimming and fill shaded with the occlusion) upsamples the same way, with the scene's
    opaque depth (`MarkingShading` carries it): in `renderer/marking/dim-occlusion-half-resolution` (dimmed dense
    spheres, multi-scale) the blocky dark creases of the bilinearly filtered packed occlusion are gone; against full
    resolution 2937 pixels off by more than 8 instead of 3298, the rest the coarser occlusion of thin creases.
  - **Resolution-independent occlusion**: below the full resolution the SSAO pass read a downsampled depth copy
    (nearest texels), which widens silhouettes: samples next to objects hit them more often and the occlusion grew
    with lower resolutions (multi-scale, against full resolution: mean -0.21 levels at 0.5, -0.68 at 0.25, 623 pixels
    off by more than 16 at 0.25). The level selection thresholds and the blur were ruled out. The pass now reads the
    full-resolution depth at any SSAO resolution (the downsampled copies and their pass are gone; the half and quarter
    levels copy from the full-resolution depth): -0.09 at 0.5, -0.11 at 0.25 (198 pixels), the rest rings at
    silhouettes from the normals and the upsampling on the coarser grid. Cost at 0.5: about +0.03 ms (scattered reads
    from the larger texture), still about 3 times cheaper than full resolution.
  - **Option 2 measured**: against a 256-sample reference (RMS of the image difference / pixels off by more than 8):
    32 samples 0.83/0 (random rotation) and 0.80/0 (4x4 ordered rotation, now in use); 16 samples 1.86/802 with either
    rotation, so the rotation is not what limits 16 samples, the kernel is; 16 of a 64-sample kernel interleaved over
    2x2 pixels 1.39/6 (half the cost of 32, large errors gone, visibly noisier than 32 in RMS); 32 of 128 interleaved
    no better than 32. Interleaving needs dynamic uniform indexing (not on WebGL1). Decided: the 4x4 ordered rotation
    stays, 32 samples stay the default (16 saves little at half resolution, see below).

## Where things stand

- **Passes per frame**: the SSAO pass (per pixel: normal reconstruction from 13 depth reads, then `samples` samples per
  level, 32 by default, one level or four with multi-scale; each sample a projection, a depth read and an unprojection),
  the separable depth-aware blur (2 passes of `blurKernelSize` taps, 15 by default), with multi-scale the half and
  quarter depth levels (copies), with transparency all of it again for the transparent depth.
- **Measured** (RTX 4070 Ti, 1280x720, both backends alike):

  | case                             | frame   | SSAO pass | blur    |
  | -------------------------------- | ------- | --------- | ------- |
  | default (32 samples)             | ~1.0 ms | 0.03 ms   | 0.03 ms |
  | multi-scale (4 x 32 samples)     | ~1.1 ms | 0.13 ms   | 0.03 ms |
  | multi-scale with 16 samples      | ~1.1 ms | 0.07 ms   | 0.03 ms |

  The pass scales with the sample count (the fixed per-pixel part, normals included, is about 0.02 ms): the cost is the
  scattered depth reads, not arithmetic. Computing the sample view z with two multiply-adds and a division instead of a
  full inverse projection changed nothing measurable here. On integrated and mobile GPUs (roughly 10-20 times slower)
  multi-scale SSAO is in the range of 1-3 ms at this size, and more on large displays (`ssaoScale` already halves the
  resolution per axis at a pixel ratio of 2).

## Options

1. **Half resolution with depth-aware upsampling**: `resolutionScale` exists, but the occlusion is packed into RG and
   upsampled bilinearly in the composite (the `multi-scale-half-resolution` known difference, halos at depth edges).
   Unpacked occlusion (`r8unorm`, filterable everywhere) and a bilateral upsample (four texels weighted by their depth
   difference) would make 0.5 usable as a default for large drawing buffers: about 4 times fewer samples and blur taps.
   The biggest gain; quality to check on edges.
2. **Fewer samples with interleaved noise**: the per-pixel rotation (`getNoiseVec2`, a hash) is random; a 4x4
   interleaved pattern (or blue noise) spreads the samples evenly over neighbourhoods the blur covers, so that 16 samples
   look like 32. Halves the pass (measured above); a default change to judge by eye.
3. **Samples per level**: multi-scale takes `samples` per level; the large levels read the quarter-resolution depth
   and vary slowly, fewer samples there (e.g. 32, 16, 8, 8) cut about half of the multi-scale cost.
4. **A float depth copy at the SSAO resolution** (`r32float`), read by the pass and the blur: reads from the
   depth-stencil texture may be slower than from a color texture on some GPUs; it also gives the blur exact depth
   ([ssao-blur-depth.md](ssao-blur-depth.md) option 4) for one extra full-screen copy. To measure.
5. **Blur**: a quarter of the cost; fewer taps where the noise is interleaved (option 2), or on WebGPU a compute blur
   with the tile in workgroup memory (one read per pixel instead of 15).
6. **Arithmetic** (view z from depth, sparse projection, cheaper normals): no gain on this GPU; possibly on ALU-bound
   mobile GPUs, to measure there before doing it.
7. **Temporal reuse while rotating** (reprojection of the previous frame's occlusion): large gain, large change and
   ghosting risks; not now.

**Recommendation**: 1 and 2 first (the largest gains, both quality-checked by eye on the usual scenes and the large
ones), then 3; 4 and 5 after measurements on a weaker GPU.

## Steps

1. Measure on a weaker GPU (an integrated adapter, `powerPreference: 'low-power'`, or a laptop) with the two benchmark
   cases, to confirm the ranking.
2. Option 1: unpacked occlusion, bilateral upsample in postprocessing and compose (GLSL, WGSL), visual cases at 0.5.
3. Option 2: interleaved noise, 16 samples compared with 32 by eye; default change recorded.
4. Option 3: per-level sample counts.

## Open questions

- Whether half resolution should become automatic above a drawing buffer size (as `ssaoScale` does for pixel ratios).
