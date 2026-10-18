# Plan: SSAO blur without packed depth

## Status

- **Current**: dropped after option 1. The linear depth in BA (option 1, done) removed the banding at usual distances
  and is good enough (checked by hand); the 16-bit packing still bands beyond about 400 Å from the camera. Reading the
  scene depth in the blur (option 4) would add a full-resolution depth fetch to each of the 15 taps of both passes
  (scattered at the default half SSAO resolution), estimated +50–100% of the blur (about +0.02–0.05 ms at 2560x1440 on
  an RTX 4070 Ti, not measured), for the far banding only. To pick up again if large scenes show it. The linear depth
  also feeds the depth-aware upsampling of [ssao-performance.md](ssao-performance.md).

## Where things stand

- **The SSAO pass** (`ssao.frag`, `ssao.wgsl`) reads the scene depth itself (the depth texture: float on WebGPU and
  WebGL2, packed into RGBA8 as 32-bit fixed point without depth textures; the multi-scale levels float `r` or packed
  RGBA8) and writes RGBA8: the occlusion packed into RG and the pixel's window depth packed into BA
  (`packUnitIntervalToRG`, 16-bit fixed point, steps of 1/65280).
- **The blur** (`ssao-blur`, two separable passes, `blurKernelSize` taps, default 15) reads only that RGBA8 target: per
  tap the occlusion from RG and the depth from BA, converted to view z, and drops taps whose view z differs from the
  pixel's by `blurDepthBias` or more (default 0.5 Å). The background and the near plane are recognized by their packed
  bytes (`isBackgroundPacked`, `isNearClipPacked`). The BA depth is read nowhere else (postprocessing, compose and the
  marking overlay read only the occlusion).
- **Precision**: perspective window depth is about `n / |z|` (reversed: near at 1), so a depth step `Δd` is a view z
  step of `|z|² / n · Δd`. With 16 bits (`Δd ≈ 1.5e-5`) and the usual near plane inside large scenes (`minNear` 5 Å):

  | distance from the camera | view z step (16-bit packed) | 32-bit float depth |
  | ------------------------ | --------------------------- | ------------------ |
  | 100 Å                    | 0.03 Å                      | < 0.001 Å          |
  | 500 Å                    | 0.8 Å                       | < 0.001 Å          |
  | 2000 Å                   | 12 Å                        | < 0.001 Å          |
  | 10000 Å                  | 300 Å                       | 0.001 Å            |

  The float column holds with a `[0, 1]` clip depth (WebGPU, WebGL with `EXT_clip_control`); WebGL without it computes
  the window depth from a `[-1, 1]` clip depth and loses the float precision far away (about `|z|² / n · 6e-8`: 0.05 Å
  at 2000 Å, 1.2 Å at 10000 Å), still about 250 times finer than the 16-bit packing.

  Beyond about 400 Å the step exceeds the 0.5 Å bias: neighbours on the same surface fall into different steps, the blur
  stops at the step edges and the unblurred noise and the step contours show as bands. Fixed point is equally coarse for
  reversed and standard depth, so reversed-Z could not help here (it helps the float depth the SSAO pass itself reads).

- **The bias is absolute**: far from the camera a pixel covers several Å, so on sloped surfaces neighbours differ by
  more than 0.5 Å even with exact depth; the blur then keeps the noise (grain rather than bands). A second, smaller
  effect, to check once the steps are gone.

## Options

1. **Linear view depth in BA** (`(-z - n) / (f - n)`, 16 bits): uniform steps of `(f - n) / 65280`, 0.15 Å for a 10000 Å
   deep scene, 1.5 Å for 100000 Å (cellular scenes): better, still above the bias in the largest scenes.
2. **Log depth in BA** (`log(-z / n) / log(f / n)`): relative steps of `log(f / n) / 65280`, 1.5 Å at 10000 Å with
   `f / n = 2e4`: no better for large scenes.
3. **Float SSAO target** (occlusion and view z as `rg32float` or `rg16float`): half float is too coarse (relative 5e-4);
   `rg32float` is not filterable without `float32-filterable` / `OES_texture_float_linear`, but postprocessing upsamples
   the occlusion with linear filtering at half resolution; WebGL without float render targets keeps packing.
4. **The blur reads the scene depth** (the depth texture the SSAO pass reads, at the SSAO resolution: the depth texture,
   or the downsampled depth target when `ssaoScale < 1`): the precision of the depth buffer (table: float), with packed
   32-bit fixed point without depth textures (0.005 Å at 10000 Å). One more texture read per tap, the BA channels and
   the packed background/near-plane tests go. The transparent SSAO blurs with the transparent depth (packed 24-bit with
   alpha, as the SSAO pass reads it).

**Decision: option 4.** It is the only one with full precision on all backends and profiles, needs no new formats, and
removes code (the depth packing in the SSAO output and the packed-byte tests). The cost is a second read per tap; to be
measured.

## Steps

1. **Blur reads the scene depth**: `tDepth` (and the transparent depth) bound to the blur renderables at the SSAO
   resolution, `getDepth` as in the SSAO pass (packed or not, `dPackedDepth`), background `depth == 0`, near plane
   `depth == 1`; the SSAO output keeps the occlusion in RG (BA unused, cleared). GLSL and WGSL.
2. **Measure**: the render benchmark with `--gpu-timing` (`render/viewer-7qpd-rotating`,
   `render/7qpd-cartoon-postprocessing`): `SSAO.blurOpaque` before and after, both backends; the visual suite (WebGL
   images change only where the packed depth stepped, WebGPU against WebGL2 within tolerances).
3. **Check by hand** on the large scenes with banding; if grain remains far from the camera, step 4.
4. **Optional, scale-aware bias**: the bias relative to the pixel footprint (`max(blurDepthBias, k * pixelSize)`, the
   pixel size the blur already computes), so that sloped surfaces far away still blur; a behaviour change of the
   `blurDepthBias` param to record.
5. **A large-scene case** to keep it checked (the depth-conventions plan's step 1 suggestion: 60 copies of 7qpd on a
   ring, camera inside, occlusion on), compared with saved WebGL images.

## Open questions

- Whether `uBlurDepthBias` should become relative altogether (step 4) rather than a lower bound.
- Half-resolution SSAO also upsamples the RG-packed occlusion bilinearly (the `multi-scale-half-resolution` known
  difference); independent of this plan.
