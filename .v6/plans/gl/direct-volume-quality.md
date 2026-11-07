# Plan: direct-volume rendering quality

## Status

- **Current**: step 4 (pre-integration) and option 6 (normals from the values) done; steps 1 (gradient-magnitude
  shading) and 2 (the occlusion volume) built and removed: neither worked well on our data. Shadows not for now. Follows
  [direct-volume-performance.md](direct-volume-performance.md): the volume pass costs 0.4–1.2 ms at 2560x1440 (RTX
  4070 Ti) on the benchmark's cases, so there is room for quality.
- **Removed**: `gradientShading` (lit by the gradient's magnitude) with the 2d grids' stored magnitudes, which are unit
  normals again; `occlusion` and `occlusionRadius` (`repr/volume/occlusion.ts`: transmittance over 14 directions per
  transfer function, into a 3d texture). What they cost, before going: occlusion 65-90 ms per transfer function change
  for a 256³ grid, about 300 ms first, the volume pass 0.34 → 0.35 ms.
- **What the shading is now** (`direct-volume.wgsl.ts`, the GLSL alike): per sample, the transfer function's color and
  opacity, and `applyLight` with the gradient as the normal (central differences of the values; WebGL1's 2d grids store
  unit normals); front to back, jittered starts (temporal multi-sampling averages the jitter). No
  shadows, no occlusion: a dense volume reads flat, its depth only from the lighting and fog. `renderMode: 'maximum'`
  instead keeps the brightest sample per ray, unlit, as bright as the transfer function makes its value (a maximum
  intensity projection); `attenuated` dims each sample by the density in front of it, `additive` sums them.
- **Constraints**:
  - the lights are in view space (`uLightDirection`, camera relative, `getTransformedLightDirection` for XR), up to
    `MaxLightCount`: lighting precomputed in the volume is stale when the camera turns; only what doesn't depend on the
    lights (occlusion) can be computed once per transfer function;
  - volumes draw in the transparent pass, after the opaque scene, stopping at its depth; SSAO and the screen-space
    shadows read the opaque depth only, so volumes neither cast nor receive them;
  - illumination (tracing) draws volumes as rasterized, outside the traced light.

## Options

Ordered by quality gained against cost and work; costs are estimates, to be measured.

1. **Gradient-magnitude shading** (built, removed: not good on our data): blend between shaded and unshaded by the gradient's magnitude (a prop, as
   other renderers do), so homogeneous regions and noise (cryo-EM maps, tomograms) don't take a random normal's
   lighting. One `length` per sample. Also fixes the zero-gradient case without the camera-facing fallback.
2. **Volumetric ambient occlusion** (built, removed: not good on our data): per voxel, the transmittance averaged over N directions
   within a radius (a few cells), through the transfer function's opacity; a compute pass into a texture beside the
   grid (half resolution is likely enough), redone when the transfer function, the grid or the radius change, not per
   frame; one more texture read per sample. The largest depth cue for dense volumes, and independent of the
   camera-relative lights. Needs compute (WebGPU), or a fragment pass per slice into a 3D texture (WebGL2).
3. **Shadows by an opacity map from the light** (per frame): the volume ray marched from the first light's direction
   into a few depth layers of accumulated opacity (deep or opacity shadow maps), at a fraction of the screen's
   resolution; each sample reads its transmittance from the layers. Its cost is about a second volume pass at the
   map's resolution, every frame the camera moves (the light moves with it). Directional, self-shadowing; the same map
   could shadow the opaque scene by the volume (not the reverse).
4. **Shadow rays per sample** (simplest, costliest): from each sample above an opacity threshold, march toward the
   first light with a coarser step and the bricks' skipping, until opaque or out. Cost about (shadow steps) times the
   primary samples that shade: 5–20x the pass for dense volumes, less for sparse maps. A quality reference for 2 and 3
   and an option for screenshots (`ImagePass`), not for interaction.
5. **Pre-integrated transfer function** (done): step 6 of [direct-volume-performance.md](direct-volume-performance.md),
   on by default (`preIntegration`): no slabs or rings at sharp transfer functions, no measurable cost.
6. **Better gradients** (done for 3d grids): normals from the filtered values (central differences, step 7 of the
   performance plan: value-only textures), not stored per cell; tricubic sampling as an option for screenshots not
   done.
7. **In illumination**: volumes in the traced light (shadowing traced surfaces, lit by the tracer); a separate plan
   if wanted.

## Recommended order

3 only if directional shadows are wanted (the lights following the camera make them less convincing); 4 as the
screenshot-quality reference.

## Open questions

- Whether shadows should also come from opaque objects onto the volume (the screen-space shadows only see opaque
  depth).
