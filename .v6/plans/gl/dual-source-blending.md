# Plan: Dual-source blending where available

## Status

- Not started. A capability for features, not a performance plan: no pass is held back by its blending today.

## What it is

The fragment shader writes two colors to one target (`@blend_src(0)` and `@blend_src(1)` at location 0); the blend
factors can read the second (`src1`, `one-minus-src1`, `src1-alpha`, `one-minus-src1-alpha`): `dst * f(src1) + src0`
with a factor per channel. WebGPU: the optional feature `dual-source-blending` (broad on desktop D3D12, Metal,
Vulkan); WebGL: `WEBGL_blend_func_extended` (narrower). A pass using it needs a fallback.

The constraint: a pipeline blending with `src1` has one color target, the only output at location 0. Passes writing
several targets can't use it.

## The passes

- **WBOIT**: two targets (weighted color and revealage; the weight sum): five values and multiple render targets. No.
- **DPOIT**: max-blended depth and back color targets, the front color by hand: multiple render targets. No.
- **Blended transparency, background composite, direct volumes**: alpha or premultiplied blending; the alpha is the
  factor already.
- **Marking overlay, path tracing accumulation, multi-sample compose**: a running average with the frame's weight as
  the blend constant (`setBlendConstant`, between passes; it flushes the pending render bundles) or additive with
  weighted colors. With `src1` the weight comes from the shader: no blend constant, and a weight per pixel. Tracing
  writes the denoiser's weight into alpha, so the alpha can't carry the factor there.

## Candidates

1. **Colored transmission** in blended transparency: `dst * (1 - src1.rgb) + src0`, a tint per channel (colored glass,
   absorption through a translucent surface), one draw, no extra target. Blended mode only (WBOIT, DPOIT: multiple
   render targets).
2. **Per-pixel accumulation weights** in the marking samples, path tracing accumulation and multi-sample compose:
   adaptive accumulation (converged pixels weighted less) without an extra target. A quality feature.

## Steps

- None until a candidate is wanted; then the feature or extension requested where available
  (`capabilities.dualSourceBlending`), the blend factors added to `BlendFactor` (both backends), the pass with its
  fallback, checked by the visual suite with and without the capability.
