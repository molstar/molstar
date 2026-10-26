# Plan: Sphere impostors rasterized by a compute kernel (WebGPU)

## Status

- **Current**: proposed, nothing done. Small opaque spheres drawn without quads: a kernel writes each covered pixel's
  nearest sphere into a visibility buffer with one 64-bit `atomicStoreMax` of depth and id (`atomic-vec2u-min-max`,
  Chrome 157), a resolve pass shades each pixel once and writes its depth. Large spheres, transparency and the other
  variants keep the impostors. Only where the device has the feature; elsewhere (WebGL, WebGPU without it) the
  impostors as today. Expected to be fast enough only with 64-bit atomics (the 32-bit alternative needs a second pass
  over every sphere): to measure. Step 1 (a prototype measured against the opaque pass) decides whether the rest is
  done.

## Why

- **Many small spheres are bound by geometry, not fragments** (`shader-analysis.md`, sphere impostors on NVIDIA: the
  PetWorld model, 3.68 M atoms, `merging/postsyn-gpu-hiz-zoom-*`):
  - 16× the pixels cost 2.3× the time;
  - no lighting changes nothing (13.0 ms), the approximate impostors are slower (15.4);
  - GPU Trace: the vertex stage stalls on allocating its output memory (ISBE) 47% of the time, the primitive engine's
    L2 traffic at 77%;
  - fewer vertex outputs gained nothing (the driver prunes them), one triangle per sphere little
    (`sphere-impostors.md`, option 8).
- **The per-sphere cost is fixed** (`screen-space-lod.md`, `merged-render-objects.md`): about 0.45 ns a sphere, about
  14.6 ms fixed plus 4 ms per Mpx zoomed, about 48 spheres a pixel; with trivial fragments 18.1 of 20.4 ms, with no
  fragments at all (quads collapsed) 13.4 ms. Levels of detail are the lever today because they draw fewer spheres;
  this plan makes each drawn sphere cheaper instead.
- **Quads overshade small spheres**: 2×2 quad shading, 21% of the quad outside the circle, and `frag_depth` disables
  early-Z (`sphere-impostors.md`, option 3), so every covered fragment is shaded, hidden or not. A resolve shades each
  pixel once.
- Compute rasterization of points is a known win where primitives are pixel-sized (point clouds drawn by atomics into a
  64-bit depth and id buffer, several times faster than the hardware's point primitives). A sphere a few pixels across
  is the same case with an exact depth per pixel.

## Today

- **Spheres** (`geo/geometry/spheres/spheres.ts`, `gl/renderable/spheres.ts`): one `tPositionGroup` element per sphere
  (center, group), pulled by the vertex stage (`assign_quad_position_group`, `vertexIndex / 6`), indexed quads from one
  shared index buffer (`quad-indices.ts`). The vertex stage sizes the sphere (`assign_size`, `uSizeFactor`, the LOD
  scale and fade), clips instances (`dClipPrimitive`), projects its bounding rectangle (`sphereProjection`, or
  `quadraticProjection` for asymmetric and XR projections) and reads every theme layer into varyings (color, marker,
  overpaint, transparency, emissive, substance, clipping). The fragment stage intersects the ray (`sphereImpostor`),
  writes `frag_depth` and shades, with per-fragment discards (miss, `fadeLod`, `clip_pixel`, `check_transparency`).
- **Drawn in every pass**: opaque color, the depth passes (`depthOpaque`, `depthOpaqueBack`, `depthTransparent`),
  pick, marking, emissive, tracing and the transparency paths (WBOIT, DPOIT, blended).
- **Compute** (`gl/compute.ts`, `gl/webgpu/compute.ts`): kernels from schemas with storage buffers (`atomic<u32>`),
  storage textures and indirect dispatch; the GPU culler (`cell-cull`, `instance-copy`, `cull-args`) compacts merged
  objects' instances per frame. Kernels cannot bind the render data pages (`tPositionGroup`, size and theme layers) or
  the object blocks; the culler packs its own instance copy (`packInstanceData`).
- **Limits of WebGPU**: no texture atomics, depth formats not storage-writable (`StorageFormats`,
  `gl/webgpu/texture.ts`). Gaussian density resolves "nearest wins" in two passes with 32-bit atomics
  (`gaussian-density.wgsl.ts`: `atomicMax` of the nearness, then `atomicMax` of the atom where it matches).
- **64-bit atomics, min and max only**: the optional feature `atomic-vec2u-min-max` (WGSL
  `enable atomic_vec2u_min_max;`, shipping in Chrome 157, positive signal from Safari; proposal
  `gpuweb/proposals/atomic-64-min-max.md`) adds `atomicStoreMax` and `atomicStoreMin` on `atomic<vec2u>` in storage
  buffers, compute and fragment stages, with no return value. They need the hardware's 64-bit atomics: Vulkan
  `shaderBufferInt64Atomics` (about 88% on Linux, 31% on Android, per the proposal), D3D12 shader model 6.6 with
  Int64 atomics, Metal on Apple 9. Not requested today (`OptionalFeatures`, `gl/webgpu/context.ts`).

## Design

1. **Which spheres**: per object by the renderer's facts, per sphere by size on screen.
   - Per object: spheres, WebGPU with `atomic-vec2u-min-max`, opaque (no transparency layer, alpha 1, not x-ray), no
     pixel clipping (`dClipSphere`; instance clipping is fine), illumination off (tracing reads the impostors' MRT
     outputs), the opaque pass of the color path. Everything else draws as today.
   - Per sphere: projected radius at most `R` CSS pixels (a prop, to measure; 4–8 expected). Both paths use the same
     test, from the same function in both shaders: the kernel skips spheres above `R`, the impostor vertex stage
     collapses those at or below it (as the LOD does, `z = 2w`). So each sphere is drawn by exactly one path, and large
     spheres, which the quads draw well, stay there. Why a split by size: the kernel's cost per sphere grows with its
     area (a loop over its pixels in one invocation), the hardware rasterizer's does not.
2. **The visibility buffer**: per pixel of the target, a depth and an id, cleared per frame, from `getScratchBuffer`
   (kept while used, released when not).
   - Depth: the reversed-Z float's bits (positive floats order as their bits, nearer is larger), so max wins.
   - Id: a sphere's slot in the frame's sphere table (item 4), `+1` (0 is empty).
   - **With `atomic-vec2u-min-max`: one pass.** An `array<atomic<vec2u>>`, the depth bits in the high word and the id
     in the low word, one `atomicStoreMax` per covered pixel: the nearest wins, equal depths pick one id
     deterministically. The feature is requested where the adapter has it (`OptionalFeatures`) and is a capability
     (`minimum-capabilities.md`: optional features only as optimizations behind a check). Which component of the
     `vec2u` is the high word: from the spec, checked by a Dawn test.
   - **Without it: the impostors.** Measured in step 1 for comparison, not planned as a path: two passes of 32-bit
     atomics as Gaussian density (the raster kernel writes depths with `atomicMax`, an id kernel runs the same loop and
     writes the id where its depth equals the stored one), which doubles the per-sphere work the kernel exists to cut.
     A single 32-bit packed pass (24-bit depth, 8-bit id within a workgroup tile) aliases ids or bands depths. Devices
     without the feature are mostly mobile, where the scenes this is for are rare.
   - Contention: about 48 spheres a pixel zoomed. A plain read before the atomic skips most writes (a sphere behind the
     stored depth does nothing), as point-cloud rasterizers do. With the 64-bit path the read has no atomic to use
     (`atomicStoreMax` returns nothing; whether `atomicLoad` takes `atomic<vec2u>` is to check). Binding the same
     buffer read-only as well conflicts with its read-write use in one dispatch, so the alternative is a 32-bit depth
     buffer next to it, written with `atomicMax` and read before the 64-bit write. Measured against no early-out.
3. **The raster kernel**: one invocation per sphere (the frame's sphere table), for each view.
   - The vertex stage's sphere: center by the instance transform, size by `assign_size` and the LOD scale, the LOD
     level by the pixel test (`screen-space-lod.md`), instance clipping, the frustum and near plane. `fadeLod`'s
     dither is per pixel, so it moves into the pixel loop (the same 4×4 pattern by pixel coordinate).
   - Its pixel rectangle from the same projection functions (`sphereProjection`, `quadraticProjection`, orthographic),
     then per pixel center the ray-sphere intersection of `sphereImpostor` (front hit, or the back hit inside the
     sphere when `dSolidInterior` or `uDoubleSided` needs it). The depth is exact per pixel, so intersecting spheres
     cut each other as with the impostors.
   - Shared chunks, not copies: the projection, intersection, size and LOD code becomes WGSL functions both the
     impostor shader and the kernel include (`wgsl-strategy.md`'s preprocessor).
4. **Input: the frame's sphere table**: the eligible objects' spheres as ranges `[base, base + count)` of one id space,
   with per range the object (its draw record), its instances and levels. An id maps back to (object, instance, item)
   by the ranges (few: merged objects make about one per representation).
   - Merged objects with GPU culling: the culler's compacted instances and level counts are the input, dispatched
     indirectly from `cull-args`. Others: their instances as stored.
   - The kernels read the render data pages (`tPositionGroup`, `tSize`, the LOD values) and the object blocks. Today a
     kernel can bind neither; this is the plumbing to add (the WebGPU data pages and uniform pages as kernel bindings,
     by the same offsets the render items use).
5. **The resolve**: per eligible object, a draw of its pixel rectangle (scissored to the bounds of its spheres on
   screen) in the opaque pass, before the other opaque objects, with the object's own pipeline and bindings in a
   `dVisibility` variant of the spheres program.
   - The fragment stage reads the id, discards pixels of other objects or none, recomputes the intersection for its
     sphere (the kernel keeps nothing per pixel but the id), writes `frag_depth` into `depthTextureOpaque` and shades
     with the spheres' material and lighting.
   - The theme layers are read per pixel by (group, instance), not by varyings: the fragment-stage versions of the
     `assign_*` chunks (WebGPU reads every layer from storage data pages, so it is a lookup, not a texture fetch per
     vertex).
   - Why per object, not one full-screen draw: the material, lighting switches and layer offsets stay the object's
     (uniform blocks, defines), as the impostors have them. One draw for all needs object blocks by index, which
     [batched-draws.md](batched-draws.md) built and removed. Few eligible objects (merged), each a scissored
     rectangle.
   - The other opaque objects, the impostors above `R` among them, are drawn after it and depth-tested against it, so
     the image composes as now. Hi-Z, SSAO, outlines, DOF and bloom read the depth as they do.
6. **The other variants keep the impostors** in their own passes: pick (on demand, at its own scale), marking,
   emissive, the depth passes, transparency. Their depth matches the resolve's (the same analytic intersection), so the
   marking's depth test and the outlines agree.
   - Later, by measurement: the depth passes from the visibility buffer (a resolve writing depth only), and picking
     from it at full resolution.
7. **Views**: one dispatch per view into the view's buffer: stereo eyes, XR views, the multi-sample jitter (the
   camera's jittered projection), `viewOffset` of tiled screenshots. The pick pass draws impostors.
8. **Plans and readiness**: the kernels are prepared with the scene (`ComputeContext.prepareKernels`,
   `pipeline-readiness.md`), the resolve variant is a program variant the precompile lists record. A prop
   (`sphereRaster: 'auto' | 'on' | 'off'`, `off` until measured) and the capability check; `auto` by GPU or scene size
   once measured.

## Steps

1. **Prototype and measure** (the gate): the kernels and a resolve with a flat shade, outside the passes, on the
   PetWorld model's spheres (`merging/postsyn-*`) and the globule (`render/globule-300k-spacefill`), at 640×360, 720p
   and 2560×1440, on NVIDIA and the AMD iGPU, against the opaque pass of the impostors. Measure the 64-bit pass (Chrome
   157, or Canary with the feature), the early-out read and `R`; the two 32-bit passes once, to confirm they are too
   slow. Continue only if the 64-bit kernel and the resolve beat the impostors for the small spheres by a margin worth
   the second path.
2. Shared WGSL chunks for projection, intersection, size and LOD; the impostors on them, images identical.
3. Kernel bindings for the render data pages and object blocks; the sphere table for unmerged and merged objects
   (indirect from the culler). Tests on Dawn: ids map back to (object, instance, item); every sphere drawn by exactly
   one path at the threshold.
4. The resolve variant (fragment-stage layer lookups) and the split in the opaque pass behind the prop. Visual cases
   on the WebGPU profiles: the globule, PetWorld whole and zoomed, a scene with spheres both sides of `R`, LOD fade,
   instance clipping, marking and outlines over resolved spheres, stereo and multi-sample.
5. Benchmark cases with the prop; `auto` decided from them.
6. Later, by measurement: depth passes and picking from the visibility buffer; tiles binned per workgroup for larger
   spheres; cylinders the same way.

## Open questions

- If the two 32-bit passes turn out close to the 64-bit pass, are they worth a path for devices without the feature?
- `R` fixed, by GPU, or from the measured costs per sphere and per pixel?
- Illumination: the tracing MRT outputs from the resolve too (color, normal, material, back depth from the back hit),
  or tracing keeps the impostors?
- Per-group transparency (`check_transparency` discards per fragment): eligible objects with a transparency layer
  split by group into the two paths, or the whole object stays on the impostors?
- Pixel clipping and the solid caps on clip surfaces (`clipExit`) in the kernel, or such objects stay on the impostors?
