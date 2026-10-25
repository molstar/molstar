# Plan: Half-precision shading on WebGPU (`shader-f16`)

## Status

- **Paused** after step 1 (measured, nothing kept): the impostor lighting in f16 is 6-11% faster on an integrated
  AMD GPU (RDNA2, packed math) and 4-12% slower on NVIDIA (RTX 4070 Ti), also with the light uniforms and the color
  varying in f16; images alike. Not neutral on NVIDIA, so it would need choosing per GPU (`adapter.info`, which
  browsers may leave empty, or timing both pipelines at runtime): not worth two numeric paths for about 10% on
  integrated GPUs now. Revisit with a GPU-choice mechanism of other uses, or where integrated GPUs are the target.
  Follows [sphere-impostors.md](sphere-impostors.md).

## Why

- WebGPU's optional `shader-f16` feature (`enable f16;` in WGSL) is widely available (D3D12 with native 16-bit
  shader ops, Metal, Vulkan with `shaderFloat16`). Where a GPU runs f16 at twice the f32 rate (packed math: AMD, Apple,
  Intel, mobile) arithmetic gets cheaper; everywhere f16 values take half the registers, which can raise occupancy of
  register-heavy fragment shaders. NVIDIA GeForce cards run plain f16 arithmetic at the f32 rate: there only the
  register savings count, so the gain differs by vendor and is measured per GPU.
- Lighting is most of the impostor's fragment work (PBR terms, light loop, fog, color mixing), and its inputs
  (normalized directions, colors in [0, 1], roughness, metalness) fit f16's range and an 8-bit output's precision.

## Design

- **Capability**: request `shader-f16` where the adapter has it (with the other optional features), a capability
  `shaderF16`; off where missing, and a renderer prop to turn it off (comparisons, problems on a GPU).
- **One source**: the WGSL preprocessor selects a type alias, `alias half = f16;` with `enable f16;` where available,
  else `alias half = f32;` (and `half3`, `half4`); shading chunks use `half` types and `half(...)` constructors instead
  of literal suffixes, so the f32 build is today's code. Pipelines of both kinds are variants (cached like others).
- **What moves to `half`**: the lighting chunk (normals, view and light directions, colors, PBR terms, the light loop),
  color, marker, overpaint, substance and emissive mixing; later full-screen postprocessing (SSAO combine, outlines,
  blur, SMAA, fog compositing) where measured worth it.
- **What stays f32**: positions and depth (the impostor ray-sphere intersection, `frag_depth`, clip space), distances
  in Å (fog, lod, clip objects), texture coordinates into large data textures, picking ids, accumulations (OIT
  weights, multi-sample sums), anything exponentiated beyond f16's range (65504; specular powers clamped first).

## Steps

1. **Upper bound** (done, not kept): `applyLight` as an f16 copy (inputs converted once, the BRDF, the light loop and
   the ambient terms in f16, `EPSILON` 1e-4), `enable f16;` prepended to the modules using it, `shader-f16` requested
   (both adapters have it in Chrome on D3D12). Opaque pass, f32 / f16:

   | case                                 | RTX 4070 Ti     | integrated Radeon (RDNA2, 2 CUs) |
   | ------------------------------------ | --------------- | -------------------------------- |
   | globule 300k, 2560x1440              | 13.9 / 15.2 ms  |                                  |
   | globule 300k, 1280x720               | 4.15 / 4.51 ms  | 299 / 268 ms                     |
   | 7qpd                                 | 0.08 / 0.08 ms  | 4.83 / 4.36 ms                   |
   | PetWorld zoomed, pixel LOD (indexed) | 5.71 / 5.89 ms  | 363 / 339 ms                     |

   Images: mean differences under 0.4 of 255, isolated edge pixels at sphere seams up to 93 (a neighbor winning the
   pixel); no banding seen. NVIDIA pays for the conversions without faster arithmetic; AMD's packed math gains.

   Where NVIDIA loses: the f32 lighting with the same values rounded to f16 and back (twice the conversions) costs
   4.71 / 15.7 / 6.09 ms on the RTX (the rows above: 4.15 / 13.9 / 5.71 f32): the conversions, about 20 per fragment,
   not the f16 arithmetic. Fewer of them: the color varying as `vec4<f16>` (allowed: f16 inter-stage IO and fragment
   outputs compile and draw on Dawn, D3D12, NVIDIA and AMD, perspective, linear and flat interpolation) and the light
   uniforms as f16 (stood in for by f16 constants, an optimistic bound); against the same constants in f32, opaque
   pass f32 / f16: RTX globule 300k 4.01 / 4.49 ms (720p), 13.4 / 15.0 (2560x1440), PetWorld zoomed 5.68 / 5.89;
   integrated Radeon globule 292 / 260 ms, 7qpd 4.74 / 4.25, PetWorld zoomed 357 / 334. NVIDIA still loses: what
   remains (the normal and view direction from the f32 ray-sphere intersection, roughness, metalness, the result,
   possibly f16 `exp2`, `sqrt` and `normalize`) can't be removed.
2. **Capability and alias**: the feature requested, `shaderF16`, the preprocessor alias, all shading still f32; a visual
   profile `webgpu-f16` (forced on).
3. **Lighting in `half`**: the lighting and mixing chunks converted; visual suite against `webgl2` within the usual
   tolerance, the measurements of step 1 again.
4. **Postprocessing** passes where measured worth it.

## Open questions

- Answered: it is slower on GeForce (4-12%), so not on everywhere; gating needs the vendor (`adapter.info`, unreliable
  on the web) or a runtime timing of both pipelines.
- Precision of PBR in f16 on smooth gradients (banding at 8-bit output); whether clamping specular terms suffices.
- Whether Tint's f16 code on D3D12 (DXC, 16-bit types) is as good as native HLSL `min16float`, and on Metal `half`.
