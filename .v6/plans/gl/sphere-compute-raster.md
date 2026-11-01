# Plan: Sphere impostors rasterized by a compute kernel (WebGPU)

## Status

- **`auto` for large structures**: merged objects need GPU culling, and `gpuCulling: 'auto'` asked for 20,000 instances,
  so PDB entries never got the raster (3J3Q: about 1,400 units). It also turns on from 100,000 spheres where the device has
  the raster (`AutoGpuCullingMinRasterSpheres`). Default against culled and raster, ms, RTX 4070 Ti: crambin 0.51, 0.63;
  6k atoms 0.90, 0.88; 300k 4.39, 3.47; 4V6X 0.56, 0.63; 3J3Q 1.48, 0.63. AMD iGPU: 300k 290, 210; 4V6X 22.4, 7.29; 3J3Q
  61.9, 11.5. Decided when objects merge, where `sphereRaster` and the camera mode are not known: with the raster off or an
  orthographic camera such objects stay culled on the GPU without it (a little slower on what was measured by hand).
- **Current**: steps 1 (the prototype, the gate passed) and 2 (shared chunks) done, step 3 done, step 4 done (the resolve and the split per instance through the renderer, draw lists, visual cases, instance clipping, prepared ahead), step 5 (occlusion) and step 6 (`auto`: objects culled on the GPU) done:
  with occlusion as the impostors have it (one phase, per instance) the raster takes 9-30 % of their opaque pass zoomed
  into PetWorld and mostly 39-57 % whole, at the default levels and at higher detail (Step 1: results, Levels by size
  on screen). Small opaque spheres
  drawn without quads: a culling pass appends chunks of the visible instances' spheres, a kernel writes each covered
  pixel's nearest sphere into a visibility buffer with one 64-bit `atomicStoreMax` of depth and id
  (`atomic-vec2u-min-max`, Chrome 157), a resolve pass shades each pixel once and writes its depth. Large spheres,
  transparency and the other variants keep the impostors; devices without the feature (WebGL, WebGPU without it) the
  impostors as today.
- **Availability** (checked on the RTX 4070 Ti, D3D12): Chrome 156 exposes `atomic-vec2u-min-max` with
  `--enable-webgpu-developer-features` (the benchmark's `--browser-arg`), not without (shipping in 157); Dawn in the
  native tests exposes it. So the prototype can be measured in Chrome now, and tested on Dawn.

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
   - Per object (`isSphereRasterEligible`, `compute/sphere-raster-pass.ts`): spheres, WebGPU with `atomicMinMax64`,
     opaque (alpha 1, no transparency layer, not x-ray), clip objects only by instance without a clipping layer (its
     flags exempt locations per group or instance; pixel clipping cuts spheres the raster can't), the exact
     impostors (not `dApproximate`), no wiggle or tumble animation, no fade by the object's distance (`fadeLod`'s
     dither, `dLodFade` with a `uLod` range and no scale), at most 4 levels; merged objects only in frames culled on
     the GPU. Per view: perspective, unscaled (`camera.scale` 1), the viewport the whole target. Only the opaque pass
     of the color variant: tracing, depth, pick, marking and emissive draw the impostors. Everything else draws as
     today.
   - Per sphere: projected radius at most `R` CSS pixels (`sphereRasterRadius`, 6 by default, to measure). Both paths
     use the same test, `rasterizesSphere`: the kernel skips spheres it rejects, the impostor vertex stage (color
     variants, the object's `uRasterMaxRadius` above 0) collapses those it accepts (`z = 2w`). The impostors test with
     a margin (`R` 0.1 % smaller, the near plane 0.1 % farther): their arithmetic is not the kernel's, so a sphere at
     the limits is drawn by both, never by neither. Large spheres, which the quads draw well, stay there. Why a split
     by size: the kernel's cost per sphere grows with its area (a loop over its pixels in one invocation), the
     hardware rasterizer's does not.
2. **The visibility buffer**: per pixel of the target, a depth and an id, cleared per frame, from `getScratchBuffer`
   (kept while used, released when not).
   - Depth: the reversed-Z float's bits (positive floats order as their bits, nearer is larger), so max wins.
   - Id: the sphere's chunk in the frame's chunks (item 3) times 64 plus its lane, `+1` (0 is empty). The chunk holds
     the rest (instance, item, level, the instance's group offsets and scale), so ids are bounded by the frame's chunks
     and fit 32 bits; a slot-based id (instance times the object's spheres) overflows for merged objects (their slots
     are levels times all members' instances, their spheres all members').
   - **With `atomic-vec2u-min-max`: one pass.** An `array<atomic<vec2u>>`, the depth bits in the high word and the id
     in the low word, one `atomicStoreMax` per covered pixel: the nearest wins, equal depths pick one id
     deterministically. The feature is requested where the adapter has it (`OptionalFeatures`) and is a capability
     (`minimum-capabilities.md`: optional features only as optimizations behind a check). Which component of the
     `vec2u` is the high word: from the spec, checked by a Dawn test.
   - The high word is `.y` (checked on Dawn): the depth in `y`, the id in `x`.
   - **Without it: the impostors.** Two passes of 32-bit atomics as Gaussian density cost about twice the 64-bit pass
     (step 1), the per-sphere work the kernel exists to cut; a single 32-bit packed pass aliases ids or bands depths.
     Devices without the feature are mostly mobile, where the scenes this is for are rare.
   - **No early-out**: a 32-bit depth read before the 64-bit write was never faster (step 1); contention is handled by
     drawing fewer spheres (occlusion, item 9) and exact pixel bounds (item 3).
3. **The raster: a culling pass and an indirect raster** (step 1's compacted dispatch, exact bounds).
   - **The culling pass**: a thread per instance of an eligible object: the frustum, the object's distance range, the
     occlusion test (item 9), the instance's levels; it appends a chunk per 64 spheres of each level it is drawn at,
     and the indirect arguments. A chunk is two vec4u: (slot, first sphere, end, level) and (`aInstance`, its
     group-instance offset, its group offset, the level's scale), in the object's region of the frame's chunks (kept
     for the resolve). Not a workgroup per instance or per fixed chunk of every instance: most would find no sphere
     under the levels (step 1: 2.5× slower). For merged objects culled on the GPU the culler has culled and levelled
     the instances: a pass over its compacted slots (`sphere-raster-regions`, item 4) appends their chunks instead.
   - **The level per sphere as the impostors choose it**: an instance at each level its sphere overlaps (as the culler
     copies it), each sphere tested in the kernel against its level's range as the vertex stage does (`assign_lod_factor`:
     distance, sizes on screen converted per view by `LodView`; the fade a smaller radius, `lodAt`; culled beyond the
     range). Not one level per instance (step 1's shortcut): the two paths draw the same object and must pick the same
     spheres. A sphere's chunk carries its level: spheres in two levels' overlap are drawn at both scales, as the
     impostors draw them, and the nearer wins. Distances are from the renderer's camera plane (`uCameraPlane`: the near
     plane's normal through the eye). Objects without an instance grid are not culled and draw their spheres whole
     with the object's `uLod` (`createRenderable`), merged members without one at each level they overlap (the
     culler's records): the raster follows each.
   - **The raster kernel**: a workgroup per chunk, a thread per sphere: center by the instance transform, size by
     `assign_size` and the level's scale; the pixels whose centers it can cover from its exact perspective bounds (the
     tangents from the eye; step 1: 3× fewer pixels tested than its radius +1, where most spheres are under a pixel),
     per pixel the ray's front hit, `atomicStoreMax`. No workgroup tile (step 1: slower). Spheres through the near
     plane or above `R` are the impostors' (`rasterizesSphere`, item 1), so the eye is never inside a rasterized
     sphere and the back hit (`dSolidInterior`, `uDoubleSided`) never needed. Hits beyond the far plane are skipped:
     a negative depth's bits would win the max. Instance clipping: the instance's center (its invariant sphere's, as
     `clip_instance_wgsl`) against the clip objects in storage (`clip_storage_wgsl`, five vec4s each, packed by
     `clip-storage.ts`; kernel uniforms have no arrays): unmerged objects in the cull kernel (`clipObjects`), merged
     ones in the culler's instance copy (after the segments in `segmentData`, the binding count at the core limit's
     8), which leaves clipped instances out of all regions, the impostors' too (the same image, less vertex work).
     Object space is the clip space: `group.view` is the identity and the impostors divide the model scale out.
   - Shared chunks, not copies, free of bindings: `sphereBounds`, `spherePixelBounds`, `rasterizesSphere` and
     `raySphere` (`chunks/sphere.wgsl.ts`), `lodAt` (`chunks/lod.wgsl.ts`), `sizeIndex` and `unpackSize`
     (`chunks/size.wgsl.ts`), `unpackRGBToInt` (`pack_int_wgsl`). Each caller reads the data its own way:
     `assign_size_wgsl` from the object's blocks and `read_tSize`, the kernel from its uniforms and the bound page.
   - **A dispatch set per object** (`gl/compute/sphere-raster.ts`): its chunk pass, arguments and raster with its
     values as their uniforms and its pages bound, one after another in a compute pass, a counter and indirect
     arguments per object. The raster kernel is the same for both kinds of object (6 storage buffers: positions,
     sizes, transforms as affine rows, counter, chunks, the visibility buffer). Perspective only: orthographic views
     draw impostors.
4. **Input: the frame's objects and their chunks** (`SphereRasterTable`): per object a region of the frame's chunks,
   an id maps back to (object, instance slot, item, level) by its chunk (few regions: merged objects make about one
   per representation).
   - Merged objects culled on the GPU: the culler's compacted instances (`MergedRasterInput`,
     `Renderable.getRasterInput` once `prepareCull` dispatched): a thread per slot finds its region (level and segment)
     in the indirect records, appends the record's spheres if within its count, at the segment's scale of the level
     (`tSegment`). The culler's frustum and levels are then the raster's. A merged object's members share one material
     (one `uSizeFactor`). Others: their instances as stored.
   - **The data pages directly, the object's values packed**: the kernels bind the pages of `tPositionGroup` and
     `tSize` (`RenderItem.getDataLayer`: the page as a read-only storage buffer, the layer's first word), so the spheres
     are read where the impostors read them. The object's values (the layer offsets, the size index type, `uSize` and
     `uSizeFactor`, the model scale, the levels' ranges and counts, the instance and group counts) go into the table's
     record per object, written from the item's values when they change. Not the object blocks: their pages are
     uniform-only, their layout is generated per program variant, and a block moves to a ring slot when it changes
     after a draw read it. Transforms: the culler's compacted instances for merged objects (storage buffers), a packed
     copy for the others (their instanced attributes are in vertex-only geometry pages).
5. **The resolve**: per drawn object, one triangle covering the viewport in the opaque pass, before the other opaque
   objects, with the object's own bindings: the color program built with `dSphereResolve`
   (`RenderItem.renderFullscreen`, WebGPU; its program cached per variant and define).
   - The fragment stage reads the id, discards pixels of other objects or none, finds the sphere by its chunk and
     rebuilds what the vertex stage outputs for the impostor (`resolveSphere`: `assign_varyings` with the chunk's
     instance and group offsets, `assign_size`, the level's scale and fade, the view ray through the pixel as the
     quad's point); then the impostor's fragment stage as it is: the intersection, `frag_depth` into
     `depthTextureOpaque`, material, lighting, marker, fog. The theme layers are read per pixel (WebGPU reads every
     layer from storage data pages: a lookup, as the vertex stage's).
   - Bound with it: the visibility buffer, the chunks and the transforms (`SphereResolveSchema`: storage buffers used
     with the define, object uniforms the renderer sets per view). Scissored to the object's bounding sphere on screen
     (`projectSphere`, the hi-Z pass's; the viewport where the sphere reaches the near plane, none off screen).
   - **Revisit**: the item hook and the define stand in for a render variant of its own; whether the renderer should
     be refactored for draws like this one (a pass of an object that is not one of the variants) is open.
   - Where an object's resolve pipeline is not ready (compiling in the background), its impostors draw all its
     spheres that frame; the kernels are prepared in the background the same way, no raster until they are.
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
8. **Plans and readiness**: the kernels are prepared in the background on first use (`SphereRasterKernels`), the
   resolve's program and pipeline on first draw; and ahead with the raster on: the opaque pass's requests
   (`renderer.collect`) carry each eligible object's resolve (`PipelineRequest.define`, the program `renderFullscreen`
   draws; a precompile list's stand-ins by spheres and opacity), and `Canvas3D.precompile` and `prepareObjects`
   prepare the kernels for lists with opaque spheres. A prop
   (`sphereRaster: 'off' | 'auto' | 'on' | 'only'`, `auto` by default; `auto` objects culled on the GPU, step 6) and the
   capability check. `only`, for testing: the raster draws every sphere of its objects whatever its size, their
   impostors are not drawn in the opaque pass (spheres through the near plane are missing). The opaque pass draws
   from draw lists with the raster on too: an unmerged object's raster values change per view (set before the pass,
   reset after, its impostors' vertex stage reads them), so its draws are walked and flushed (a few rows, drawn
   directly); a merged object culled on the GPU keeps its values until a view no longer draws it (its impostors read
   none) and leaves out the raster regions by `Renderable.rasterDrawn`, set for the opaque pass only and part of its
   `drawVersion`, so the pass replays unchanged (`unplanned`). The culler's lists without the raster regions are the
   same objects as the levels' with them (`getLists`): else each call made new lists, new draw blocks and a recording
   per frame. `only`'s objects are left out of the walk as of the draws.
   - The camera's projection is a GL projection: the kernels bind it with the WebGPU programs' reversed `[0, 1]` depth
     (`ReversedZeroToOne`), so a hit's depth is its `frag_depth` and the max is the nearest. The camera plane is
     the renderer's (world space), into the objects' space by the model matrix.
9. **Occlusion**: as the impostors' (one phase, per instance): the culling pass tests each instance against the
   pyramid GPU culling already reads (`HiZPass.gpuPyramid`, the previous frame's, with its camera): built from
   `depthTextureOpaque` after the opaque pass, so the resolved spheres (their `frag_depth`) and every other opaque
   object occlude. Step 1: like for like the raster's gain; its popping about the impostors' (375-492 pixels zoomed a
   2° step), not visible in normal use. Two phases (the visible instances first, a pyramid of their buffer, the rest)
   cost the raster 3-20 % over one and are exact: a later option, all in compute (no pass split), but its pyramid
   would hold the raster's spheres only unless built from the combined depth.
10. **Which path where**: with occlusion the raster was ahead in every case measured (at worst whole at 2560×1440 at
    the default levels, 82-93 % of the impostors); no choice per object or frame at first, the prop decides (item 8).
    `auto` from the benchmark cases, by GPU or by what is drawn (spheres per covered pixel) if a case loses.

## Steps

1. **Prototype and measure** (done, the gate passed): `tests/src/benchmark/sphere-raster.ts` (below).
   - **Follow-up measurements** (to do, before or alongside step 2):
     - **The integrated GPU** (done in part): the AMD iGPU (RDNA2, 2 CUs, `--force_low_power_gpu`) exposes
       `atomic-vec2u-min-max` (Chrome 156 with `--enable-webgpu-developer-features`). PetWorld is too large for it but
       one case (zoomed, px 3, 720p: the impostors 165 ms, the raster with exact bounds 114, with one phase 24.6, 15
       %); 3j3q as spacefill instead (`sphere-raster/3j3q*`, 2.4 M atoms, a target per chain: no levels, culled per
       chain), the impostors / the raster with exact bounds / with one phase, ms: whole 720p 59.3 / 7.4 / 7.7 (13 %),
       zoomed 720p 110.7 / 33.3 / 30.5 (28-30 %), whole 1440p 130.9 / 22.6 / 24.8 (17 %); zoomed 1440p lost the device
       (the prototype's frame submits every variant at once: over the system's GPU timeout there). On the RTX 3j3q:
       1.18 / 0.21 (18 %) whole, 1.49 / 0.85, one phase 0.74 (50 %) zoomed. The raster's margin is larger on the
       iGPU, the impostors' cost per quad weighs more there.
     - **Profiles** (NVIDIA in part; AMD to do): the raster, culling and resolve kernels and the impostors' opaque
       pass under the vendors' tools. Nsight Graphics GPU Trace on the RTX (`shader-analysis.md`'s procedure, 500 ms
       at base clocks; 2 s overflowed its sampling buffer and hung), the case looped alone by the page's
       `?rasterProfile=raster` or `impostors` (`BenchmarkEnv.rasterProfile`: the one-phase raster and the resolve, or
       the impostors' frames, without end; PetWorld zoomed, px 3, 720p). The command line exports totals over the trace
       only; per shader and per line is the Nsight UI's (Shader Pipelines, Shader Source by samples; `emit_hlsl_debug_symbols`
       for readable lines).
       - The raster (the GPU 58 % active): stalled on L1TEX long scoreboard 13 % and LSU queue throttle 10 % (the
         per-pixel `atomicStoreMax` and the sphere loads), L2 30 % (about half the active time: where the atomics
         resolve), DRAM 12 %, ALU 19 %, FMA 10 %; warp launches stalled on register allocation 10 %. Bound by memory
         and atomic traffic, then registers; not by arithmetic or loop control.
       - The impostors (78 % active): SM 60 %, L2 37 %, the world pipe 32 %; pixel shader warps 26 %, vertex 11 %, 40 %
         of warp slots inactive in active SMs.
       - Unexplained: the trace's raster frame about 0.26 ms of GPU (1,107 frames in 500 ms at 58 % active), the
         benchmark's timestamps 0.71 (alone as in the mix): to read from a frame's range in the UI.
       - Radeon GPU Profiler for the iGPU: not installed (the Radeon Developer Tool Suite, captured by hand).
2. **Shared WGSL chunks** (done) for projection, intersection and the level per sphere; the impostors on them, images
   identical (the WebGPU visual cases).
3. **Kernel bindings and the sphere table** (item 4):
   - a. (done) A render item's data layers as kernel bindings (`RenderItem.getDataLayer`, `getStorageBufferView`);
     Dawn: a kernel reads the positions and sizes the item draws with.
   - b. (done) The sphere table for unmerged objects (the packed values and transforms), the culling pass appending
     chunks, the indirect raster kernel into the visibility buffer (`createSphereRaster`), the `atomicMinMax64`
     capability. Dawn (`webgpu-sphere-raster.test.ts`): ids and depths match a CPU raster of the same maths (no pixel
     differs), ids decode to (object, instance slot, item, level), spheres through the near plane and above `R` left
     out; each instance at the level its distance picks. That the impostors leave out the raster's spheres is step 4
     (their side of `rasterizesSphere`).
   - c. (done) Merged objects: the culler's compacted instances as the input (`MergedRasterInput`,
     `sphere-raster-regions`), ids by chunk (item 2). Dawn: a merged object draws the same depths as its members
     unmerged, without levels and with two faded into each other.
4. **The resolve and the opaque pass** behind the prop:
   - a. (done) The resolve (`dSphereResolve`, `RenderItem.renderFullscreen`), the impostors' skip, the renderer's
     `sphereRaster` and `sphereRasterRadius` props, `prepareSphereRaster` after `prepareCull` in the draw pass
     (`compute/sphere-raster-pass.ts`). Dawn (`webgpu-sphere-resolve.test.ts`): spacefill with the raster on draws
     the image of the impostors (multi-sample and antialiasing as by default).
   - b. (done) **The split per instance** for merged objects: the impostors' vertex stage collapsing the raster's
     spheres cost about what drawing them did (Step 4a: through the renderer). The culler's copy kernel classifies
     each instance per drawn level (`rasterized`, `cell-cull.wgsl.ts`): its segment's largest sphere (`padding`) at
     the level's scale at most `R` pixels at the instance's nearest point (the tangents from the eye), the instance in
     front of the near plane, into a raster region after the levels' (the same compacted buffers, counters and
     records: `MergedRasterParams`, from the renderer's `update` to `prepareCull`); otherwise into the impostors'
     region as before. `segmentData` holds three vec4s a segment (the sphere; the first instance and the largest
     radius; the level scales), no new binding. The impostors draw their regions without the vertex skip
     (`uRasterMaxRadius` 0 for merged objects), and the raster regions too where the raster does not draw the object
     in a view (`getLists(levels, withRaster)`: another view, a pipeline not ready); the raster draws every sphere of
     its regions (`sphere-raster-regions` from `uSlotBase` and `uRecordBase`). One path per instance. The distance
     classified at is the larger of the instance's nearest point and the level's minimum distance (the level culls its
     spheres closer: without, instances overlapping a coarse level from nearer stayed impostors). The raster draws a
     merged object where the view draws the GPU culling's lists, as the object's `cull` decides
     (`getRasterInput(frame, cameraPlane, frustum)`): a frame culled already keeps its lists in its later views (the
     samples of multi-sampling, a resting camera's temporal samples), else the view culled on the GPU; a raster of the
     culler's regions where the CPU culling's lists draw every instance would double the work, and no raster where
     the GPU's lists leave the raster regions out would draw them all as impostors (the viewer at rest, multi-sampled:
     the raster pass 0 ms most frames). `only` culls the objects it skips, so that their later views keep the GPU's
     lists. The raster
     regions exist from the first dispatch with raster parameters (a layout change). Unmerged objects keep the
     per-sphere split (few instances). Dawn (`webgpu-sphere-raster-split.test.ts`): near instances as impostors, far
     ones by the raster, both kinds of region filled, the image of the impostors.
   - c. (done) Draw lists for the opaque pass with the raster on (see 8.): its direct draws issued the merged
     object's indirect records one by one, 619 draw calls for 18 (1219 while the lists alternated); now 19, the
     frame of `on` that of `only` whole and zoomed. The resolve scissored to its object (see 5.). Dawn
     (`webgpu-sphere-raster-split.test.ts`): the same view again replays, nothing drawn directly or recorded.
   - d. Visual cases on the WebGPU profiles with the prop on (done, `renderer/sphere-raster/*`, against `webgl2`'s
     impostors, at most 7 px of 65536 differ): the globule (spheres on both sides of `R`, close up, `only`), copies
     of a globule merged and culled on the GPU (PetWorld is a file on demand) with and without two levels faded
     into each other, marking and outlines over resolved spheres, multi-sample, stereo. The canvas captured (the image
     pass's target compiles the resolve again), after frames until the raster drew (`stats.sphereRasterResolves`,
     an error where it never did on a device with it); no `fallback` (it renders the target with the raster off
     too). The runner passes `--enable-webgpu-developer-features` for the feature. Instance clipping (done, see 3.):
     `clip-instance` (unmerged, `on`) and `merged-clip-instance` (`only`: the clipped copies are the near ones, which
     `on` leaves to the impostors); each fails without its kernel's test. Prepared ahead (done, see 8.): Dawn
     (`webgpu-prepare-sphere-raster.test.ts`), the first frame after the commit compiles nothing and draws the raster,
     with what representations prepare while they build and with a recorded list; neither without the requests.
5. **Occlusion** (done): merged objects culled on the GPU had it already, their raster regions of the culler's cells
   and instances tested against `HiZPass.gpuPyramid`; unmerged objects' instances now in the cull kernel against the
   same pyramid (`SphereRasterView.pyramid`, `Renderer.setOcclusionPyramid`'s; their impostors are culled on the CPU
   by `HiZPass.isOccluded`, the CPU's pyramid: the two may differ for a frame). The test a shared chunk
   (`hiz_occlusion_wgsl`), the schema, the empty texture and the uniforms shared too (`HiZSchema`, `createEmptyHiZ`,
   `setHiZValues` in `cell-cull.ts`). Dawn (`webgpu-sphere-raster.test.ts`): a pyramid at the near plane leaves
   nothing drawn, at the far plane the image without one. `occlusion-lag/*-raster-*` (on with 6 px): the raster's
   popping is the impostors' (the same culling), pixels of 921,600 differing from no occlusion, impostors / raster:
   whole 2° 0 / 24, 10° 18 / 36, 30° 3905 / 3582, 90° 330 / 352; zoomed 2° 41 / 64, 10° 114 / 142, 30° 43 / 49, 90°
   683,640 / 683,429. The frame rotating with occlusion, zoomed, impostors / raster: 2° 5.45 / 1.14 ms, 10° 6.06 /
   1.33, 30° 6.59 / 1.45, 90° 9.10 / 1.90.
6. **Benchmark cases with the prop** (done): `sphere-raster-scene/*` (spacefill with the default props, Hi-Z on: a
   crambin, globules of 6k and 300k atoms, 60 copies of 7qpd unmerged and merged and culled on the GPU; whole and
   zoomed, 1280×720 and 2560×1440; `off`, `auto`, `on`, `only`; the cases now at the default limit, this table at 6 px) beside `sphere-raster-renderer/postsyn-*`.
   Frame medians, off / on / only, ms:

   | case                | whole              | whole 2560×1440    | zoomed             | zoomed 2560×1440   |
   | ------------------- | ------------------ | ------------------ | ------------------ | ------------------ |
   | crambin             | 0.42 / 0.59 / 0.75 | 0.58 / 0.64 / 2.43 | 0.66 / 0.69 / 19.3 | 1.43 / 1.53 / 76.0 |
   | globule 6k          | 0.83 / 0.89 / 2.50 | 2.70 / 2.87 / 9.33 | 3.18 / 3.26 / 99.1 | 12.0 / 12.2 / 395  |
   | globule 300k        | 4.34 / 4.48 / 3.47 | 14.4 / 14.7 / 12.6 | 18.0 / 18.1 / 93.0 | 68.5 / 68.9 / 372  |
   | assembly 60         | 0.98 / 1.31 / 1.31 | 2.48 / 2.31 / 2.04 | 4.46 / 4.15 / 117  | 16.2 / 16.5 / 466  |
   | assembly 60, merged | 0.93 / 0.87 / 0.88 | 2.51 / 2.22 / 1.82 | 3.92 / 3.95 / 1.40 | 14.5 / 14.6 / 5.08 |

   - Unmerged objects (the split per sphere) gain nothing with `on`, at worst 40 % slower (the crambin's fixed cost):
     their quads' vertex stage runs for every sphere, the raster's work added. `only` wins where every sphere is a
     few pixels (the 300k globule whole) and is 20-80× slower close up (a thread per sphere over all its pixels).
   - Merged objects culled on the GPU (the split per instance) never lose with `on` (within noise, 2.22 against 2.51
     at 2560×1440), PetWorld up to 7× (step 4b).
   - Hence **`auto`**: the raster draws objects culled on the GPU only (`SphereRasterPassView.culledOnly`; the
     requests of `collect` and the kernels of `Canvas3D.precompile` alike, a precompile list's stand-ins where
     `gpuCulled`); measured as `off` on unmerged objects and as `on` on merged ones (`postsyn-zoom-viewer-auto` 1.47
     ms, `off` 3.35). The default `auto` (devices without the feature draw impostors).
   - **The integrated GPU** (AMD RDNA2, 2 CUs, `--force_low_power_gpu`; PetWorld too large for it):
     `sphere-raster-igpu/3j3q-*`, 3j3q as spacefill merged and culled on the GPU, the default limit (frame / raster
     pass, ms; off / auto / only): whole 62.4 / 12.2 (7.80) / 11.3, whole at 2560×1440 141 / 34.7 (21.0) / 35.1,
     zoomed 118 / 34.2 (23.2) / 33.1, zoomed and cut by a plane 67.0 / 16.5 (8.79) / 15.7; 60 copies of 7qpd merged
     whole 54.6 / 11.8 / 11.0, zoomed 252 / 253 / 31.1. `auto` 3.4-5× the impostors; the RTX (same cases) 1.1-1.4×
     whole, zoomed within noise. The merged copies zoomed stay with the impostors at any limit (RTX, the opaque pass
     of `auto` 3.14 ms at 24 px, 3.10 at 64, `off` 3.68): their instances, chains reaching close to the camera or
     through the near plane, go whole to the impostors by the split per instance, which `only` does not draw (its
     0.98, the iGPU's 31.1, less work, not a faster path): step 8's leftover. A case at a time with a timeout (2
     min): a run of them all held the iGPU over 20 min. The unmerged scenes on the iGPU not measured (`auto` draws
     them as `off`).
7. **Pixel clipping** (done): objects clipped by pixel (`dClipVariant` 'pixel', also `dClipPrimitive`) drawn by the
   raster, now left to the impostors (`isSphereRasterEligible`). What the impostors do (`sphereImpostor`,
   `dClipSphere`): the front hit clipped (`clipTest` at its model position), with a solid interior (`dSolidInterior`,
   the spheres' default) the first point inside the sphere along the ray that no clip object removes (`clipExit`, a cap
   on the clip surface, its normal the clip's), without one the back hit where double sided, else no hit;
   `dClipPrimitive` clips a sphere whole by its center in the vertex stage; a clipping layer's flags (`dClipping`, per
   group or instance) exempt clip objects per location. The raster must decide the same per pixel: the visibility
   buffer keeps the nearest depth, so a clipped front written would hide what is behind it (another sphere, its own
   cap).
   - a. (done) **Clip functions on accessors**: `clip.wgsl.ts` reads the clip objects through functions (frame, shape, count,
     active) that the programs define from `uMaterial` and `uDraw` and the kernels from storage
     (`clip_storage_wgsl`, five vec4s each, `clip-storage.ts`): one `clipTest`, `clipExit` and `clipNormal`, not a
     copy (the cap's depth must be the resolve's). Programs unchanged: the visual suite and the clip benchmark
     (`clip/globule-300k-*`) as before. Done: `clip_core_wgsl(count)` on `clipObjectFrame`, `clipObjectShape` and
     `isClipObjectActive`, the count substituted (the programs' `oClipObjectCount`, unrolled as before; the kernels'
     uniform); `clip_storage_wgsl(buffer, base, count)` the kernels' (the culler's and the raster's instance clipping
     on `clipTest`). The clip cases on WebGPU pass; `clip/globule-300k-spacefill-large-1` / `-3` 18.3 / 25.5 ms, as
     before (18.4 / 25.5).
   - b. (done) **The raster kernel**: per sphere, `dClipPrimitive` by its center (as the vertex stage, skipped whole); per
     pixel, the front hit tested, clipped: the cap's point (`clipExit` between the front and back hits, from 0: the
     raster leaves spheres through the near plane to the impostors, so no near cap) or the back hit (double sided)
     or nothing, its depth written; the object's clip objects in a storage buffer per object (as instance clipping's,
     unmerged and merged alike: the culler does not test pixel clipping). Flags of a clipping layer: the sphere's
     from the `tClipping` data page by its group and instance (as the size's, `sizeIndex`), or such objects stay on
     the impostors at first (decided by what a layer costs to read per sphere).
   - c. (done, unchanged) **The resolve**: the impostor's fragment stage as it is (`dClipSphere` with `dSphereResolve`): the hit of the
     sphere the raster chose, its cap or back hit recomputed, `frag_depth` and the cap's normal; the raster's choice
     and the resolve's agree where both evaluate the same functions on the same ray (the raster's ray through the
     pixel center, item 5). A pixel the raster found clipped through and through is not written, so the resolve never
     sees it.
   - d. (done) **Eligibility**: pixel clipping with or without solid interior, with a clipping layer (its flag per
     sphere in the kernel, `clipFlag`: the `tClipping` data page at the impostors' index, `uInstanceGranularity`'s;
     the object left to the impostors where the page is missing); not instance clipping with a clipping layer
     (`isStorageClipSupported`: the flags per group exempt an instance's spheres apart, which the culler's test per
     instance can't). The kernels' counts by variant (`getInstanceClipCount`, `getPixelClipCount`), the clip objects
     packed for both (`getStorageClipCount`), a buffer per object, merged ones too.
   - Per sphere first (`clipSphereClass`): from the signed distance of its center to each active clip object no point
     clipped, every point (skipped), or per pixel only where it reaches a surface; the distances bounds once divided by
     their largest slope (1; the cone's `|size.xy|`; the ellipsoid's distance scaled by its smallest axis), the radius
     by the frame's largest scale. (The PetWorld clip cases measured here first had their clip objects at the origin,
     off the model: the plane clipped it whole. Their numbers, 21.0 ms without the classification, 1.28 with it, the
     ellipsoid and cone per pixel 2.00 / 1.98 and with bounds 1.46 / 1.43, are of that; the cases' objects now at the
     model's center, below.)
   - e. (done) **Tests**: Dawn (`webgpu-sphere-raster-clip.test.ts`: images of `only` against the impostors per mode;
     the backs, nothing behind single-sided fronts, a hole through a grid of atoms, whole spheres by their centers
     fail without the kernel's clipping; a convex region's caps would not, the resolve finding them from the clipped
     fronts), not against a CPU reference as first planned with a plane and a sphere clip object (front clipped, cap,
     back hit, clipped through); visual cases on the WebGPU profiles against the impostors: the globule cut by a
     tilted plane (solid interior on and off, inverted), `dClipPrimitive`, a clipping layer exempting a group; merged
     copies cut by a plane.
   - f. (done) **Measure** (`clip/globule-300k-spacefill-large-*`, `clip/postsyn-zoom-plane-*`): frame medians, ms,
     impostors / raster: the 300k globule at 2560×1440 by one clip object 18.2 / 12.7 (`only`), by three 25.5 / 12.7;
     PetWorld zoomed (`clip/postsyn-zoom-*`, the clip objects at the model's center) cut by a plane 5.12 / 2.87
     (`auto`), by an ellipsoid (its inside removed) 3.94 / 1.50, by a cone 3.97 / 1.48 (with the GPU culling of
     clipped objects and the culler's rounds by new pyramids, [clip-culling.md](clip-culling.md) step 3: 1.75 / 1.34,
     3.89 / 1.53, 3.92 / 1.53); `on` on the unmerged globule
     as `off` (step 6).
8. **What zoomed in still goes to the impostors** (done: mixed instances built, measured and dropped; large spheres
   drawn by the workgroup, the limit 24 px): `sphere-raster-zoom/postsyn-viewer-*` (the viewer's defaults; frame / raster pass / opaque pass, ms):

   | case      | off  | auto 6 px          | auto 24 px         | auto 64 px         | auto 256 px        | only               |
   | --------- | ---- | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ |
   | whole     | 0.86 | 0.93 / 0.30 / 0.03 | 1.06 / 0.39 / 0.04 | 0.91 / 0.41 / 0.04 | 0.93 / 0.41 / 0.04 | 0.83 / 0.40 / 0.02 |
   | zoom 0.3  | 3.33 | 1.45 / 1.05 / 0.05 | 1.43 / 1.03 / 0.05 | 1.41 / 1.04 / 0.06 | 1.47 / 1.05 / 0.05 | 1.48 / 1.05 / 0.04 |
   | zoom 0.1  | 5.65 | 3.40 / 1.38 / 1.62 | 3.22 / 1.47 / 1.37 | 3.23 / 1.47 / 1.36 | 3.21 / 1.50 / 1.29 | 2.03 / 1.62 / 0.04 |
   | zoom 0.03 | 5.99 | 3.71 / 1.36 / 1.96 | 3.26 / 1.60 / 1.25 | 3.25 / 1.62 / 1.20 | 3.20 / 1.63 / 1.20 | 2.05 / 1.64 / 0.04 |

   - Measured with the large spheres drawn by the workgroup and the culler counting rounds by new pyramids (below);
     runs of a few cases. Frames of about a millisecond vary up to 2× from run to run on the RTX with its clocks: the
     view-sized passes (SMAA, SSAO, Hi-Z) slower by as much tell such a run.

   - A larger limit leaves the impostors 1.2-2.5 ms at 0.1 and 0.03. `only` here is not the raster drawing everything:
     the instances the split per instance sends to the impostors (near the camera, through the near plane) `only`
     does not draw at all, so its 1.8 ms lacks the costly spheres.
   - **Mixed instances** (built and dropped): instances some of whose spheres the raster may draw (at their farthest
     point) copied into mixed regions after the raster's, their spheres split one by one (the raster's chunks
     flagged in the level word and tested against the limit, the impostors' lists of those regions leaving out per
     sphere what the raster draws by `uRasterMaxRadius` as a draw uniform). Zoom 0.03: the opaque pass 1.98 → 1.62 ms
     at 6 px, 1.25 → 1.22 at 24, the raster's pass up (1.37 → 1.64, in runs of the case alone 4.1); `only`, now
     drawing the near instances too, 6.25 ms of raster against the impostors' 5.6. The impostors' share is spheres
     large on screen near the camera, impostors whatever the split, fragment-bound; leaving out spheres in the vertex
     stage does not touch it, the mixed spheres added to the raster's work.
   - **Large spheres by the workgroup** (done, the original option in its simplest form): a sphere covering more than
     `SphereRasterLargePixels` (64) pixels is queued in the chunk's workgroup memory and drawn by all 64 threads, its
     pixels strided over them, after a barrier (no second kernel, no buffer: a chunk's spheres fit the queue; the
     per-thread early returns a `valid` flag, the per-pixel work `rasterPixel`). Dawn: the raster test against the CPU
     reference with no limit passes as before. `only`, every sphere by the raster, ms (impostors / before / after):
     the 300k globule zoomed at 2560×1440 68.7 / 372 / 33.4, at 1280×720 18.0 / 104 / 8.67; 60 copies merged
     zoomed at 2560×1440 14.4 / 5.15 / 3.12, unmerged 16.3 / 466 / 79.8; the 6k globule there 12.2 / 395 / 91.2; the
     crambin 1.43 / 76.0 / 23.0. Spheres of tens of pixels now cheaper by the raster than as impostors; of hundreds
     still not (every pixel of every sphere a ray test and an atomic, overdraw and contention, against the impostors'
     fixed-function quads), which the limit leaves to the impostors.
   - **The limit raised to 24 px** (`sphereRasterRadius`, from 6; at most 64): PetWorld as the viewer draws it (frame,
     ms) zoomed to 0.1 3.40 at 6 px, 3.23 at 24, 3.23 at 64, 3.22 at 256; to 0.03 3.77, 3.29, 3.25, 3.22; whole and
     at 0.3 unchanged (0.83, 1.43-1.46).
   - Left: the impostors' share zoomed in, the instances near the camera sent whole by the split per instance.
     Measured (`--raster-split-notes`: the split replayed on the CPU, no occlusion; `auto` at 24 px): PetWorld zoomed
     to 0.1, 7 instances to the impostors (6 through the near plane), 8.4k spheres in view, 789 of them each too large
     for the raster, 7.6k each rasterizable, their discs 8.1× the viewport, the opaque pass 1.35 ms; to 0.03 5.0k,
     425, 4.6k, 5.9×, 1.24 ms. 60 copies of 7qpd merged, zoomed: 102 instances (52 through the near plane), 16k
     spheres, 15k each too large, 791 rasterizable, 142× the viewport, 3.13 ms (at 2560×1440 137, 26k, 15k, 10k,
     144×, 12.3 ms); 3j3q none. The impostors with every sphere each rasterizable left out (`uRasterMaxRadius` for
     merged objects, the image then wrong): PetWorld 1.23 and 1.20 ms, the copies 3.11 and 12.1. So not the split: a
     split per sphere, by size or through the near plane (both built as prototypes: the instances to the raster and
     the impostors, each sphere drawn by one), leaves the opaque pass at 1.2 ms at any limit up to 256 px. Occlusion
     per sphere removes nothing ([sphere-occlusion.md](sphere-occlusion.md)).
   - **The cost: spheres reaching behind the eye** (done, the fix): the impostors' quad, the bounds of the sphere's
     projection (`sphereBounds`, for a sphere in front of the eye), is unbounded for a sphere across the camera
     plane, a quad over the view, every fragment a ray test that misses or lies before the near plane. Such spheres
     left out, the image as before but for a pixel: PetWorld 1.37 → 0.24 ms, those around the eye or wholly behind the
     near plane no change. Now the quad of a sphere through the near plane is the bounds of its part beyond it
     (`sphereNearBounds`, per axis the tangent points beyond the plane and the ends of the chord it cuts, as Mara and
     McGuire 2013 before its specialization; none where no part is beyond), at a fixed depth in the clip volume; spheres
     wholly beyond the plane as before. Frame (the opaque pass), ms, before → after, WebGPU: PetWorld zoomed to 0.1
     `auto` 3.23 (1.37) → 2.12 (0.25), `off` 5.70 → 4.55; to 0.03 `auto` 3.25 (1.25) → 2.18 (0.17), `off` 5.99 → 4.99;
     60 copies merged zoomed `auto` 3.67 (3.14) → 2.57 (2.05), at 2560×1440 13.6 (12.3) → 9.31 (7.94), `off` 3.97 →
     2.90, 14.6 → 10.2; unmerged 4.50 → 3.42, 16.3 → 11.9; the 300k globule zoomed 18.1 → 12.2, at 2560×1440 68.6 →
     45.6. WebGL2 (`off`): the globule 17.0 → 10.9, PetWorld to 0.1 10.6 → 8.25, the copies merged 4.28 → 3.12,
     unmerged 4.30 → 3.13. The visual suite as before on WebGL2, WebGL1 and WebGPU;
     `renderer/sphere-near/*` (the camera among spheres, the eye in one, solid interior, open, double sided) the images
     of the old bounds.
   - **A cellPACK model** (`sphere-raster-mg/*`, `MG_1184_cellpack_atom_instances.bcif` as the viewer draws it: 848
     units, 2.1 M instances of about 76 spheres, 161 M spheres at full detail, a scene radius of 1,709 Å; PetWorld 150,
     0.52 M, about 4,600, 2,389 M, 9,526 Å). Frame (raster pass / opaque pass), ms, `off` / `auto` / `only` / the
     impostors culled on the CPU: whole 1.93 (- / 0.66) / 2.06 (1.50 / 0.11) / 1.83 / 18.0; zoomed to 0.3 4.06 (- / 3.63)
     / 2.29 (1.58 / 0.24) / 1.97 / 15.5; to 0.1 3.53 / 2.13 / 1.85 / 8.20; to 0.03 3.56 / 2.28 / 1.82 / 8.17; at
     2560×1440 whole 3.81 / 3.45, zoomed to 0.3 14.0 / 6.60; rotating whole 2.45 / 2.33, zoomed to 0.3 4.60 / 2.76.
     Storage 886 MB, with the raster 1,437 MB. The raster's pass about the same at every zoom, whole slower than the
     impostors' opaque pass; about 1.15 ms of the frame with `off` in no timed pass (not the culler: a still camera
     does not dispatch).
     - Its raster pass by kernel (RTX, whole): appending the chunks 0.20 ms, the args 0, the raster 1.29, of 286,791
       chunks (283k instances, nearly all at the third level): 4.5 ns a chunk, as PetWorld whole (0.20 ms, 44,705 chunks,
       43.9k instances). So whole the cost is the workgroups, a chunk or so per visible instance, not pixels; closer up
       and larger the pixels add (zoomed to 0.3 153k chunks 1.38 ms, at 2560×1440 4.56).
     - How full the chunks are (the culler's raster regions, no occlusion by the kernel counted): MG_1184 whole 278k
       instances, 616k spheres in 281k chunks, 3.4 % of the threads with a sphere, 272k of the instances one sphere (the
       coarse level), 9.6k chunks if full; zoomed to 0.3 69.5 % (154k chunks, 107k if full), to 0.03 80.8 % (106k, 86k).
       PetWorld whole 23.5 % (44.7k, 10.5k), zoomed to 0.1 93 %. So whole, a workgroup per sphere; sharing chunks
       between small instances could remove most of the workgroups there, little zoomed in.
     - **Instances of one sphere** (done): a row of its own from the end of the object's region of the table, its lane 0
       the sphere (the resolve as for chunks, at most a row for a chunk), a second dispatch of the raster kernel 64
       rows a workgroup (`uSingles`). MG_1184 whole: the raster pass 1.50 → 0.47 ms, the timed GPU work 1.73 → 0.72, the
       frame 2.06 → 1.83, rotating 2.33 → 1.70; zoomed and PetWorld within noise. The images as before but for 2 pixels
       by little whole, none zoomed to 0.3. Instances of a few spheres (zoomed to 0.3: 26k of 4 or fewer, of 154k
       chunks) not yet: more rows than chunks, the capacity then to change.
     - **The chunks' buffer follows the views** (done, merged objects): sized for every instance at its finest it took
       (32 bytes a chunk) MG_1184 10.3 M chunks, PetWorld 42.9 M, where the views measured use at most 287k and 123k. Now
       as many as instances, at least 65,536, the counters read back after a dispatch doubling them where a view uses
       more than half; a view using more than all loses spheres for that frame, the read pending work (the canvas draws
       again, an image pass renders again). Storage at the end of a run: PetWorld 1,668 → 376 MB, MG_1184 1,758 → 1,509.
       Unmerged objects cull every frame, keep their bound. (The cases' byte notes after the first frame can precede the
       raster, without this buffer.)
       `renderer/sphere-raster/chunks-growing-image`: an image pass needing 192k chunks, the first render
       losing copies, renders again, the image that of the chunks for every instance from the start, to the pixel.
     - The frame beyond the timed passes (RTX, whole, about 1.1 ms; PetWorld about 0.3) is the indirect draws: a draw per
       level and unit on the CPU's side whether empty or not. The culler's lists emptied (no impostors, the image then
       wrong): `auto` 1.63 → 0.80 ms, `off` 1.43 → 0.62, though the opaque pass timed 0.09 of it; a cost of issuing,
       beside the GPU's work.
       Of 3,392 impostor draws (848 units, 4 levels) with `auto` none has instances whole, 2 and 12 zoomed to 0.3 and 0.03;
       PetWorld 0 and 6 of 600.
     - **Empty draws left out at rest** (done): a still camera's records stay (`isDispatched`), read back once; each list
       then draws only those with instances (`MultiDrawIndirect.records`), any dispatch restores them. Frame, ms: MG_1184
       whole `auto` 1.63 → 0.81, `off` 1.43 → 0.87, zoomed to 0.03 `auto` 2.28 → 2.01; PetWorld whole `auto` 0.79 → 0.67,
       `off` 0.75 → 0.68. So at rest `auto` ahead of `off` whole too. Recordings follow the records
       (`recordsVersion` among the inputs of a replay without a walk).
     - **Multi-draw indirect while moving** (done, `chromium-experimental-multi-draw-indirect`, Chrome with
       `--enable-unsafe-webgpu`, capability `multiDrawIndirect`; the loop of draws without): one call per run of
       records. On render passes only (Chrome 156: not on bundle encoders), so a chunk with an object culled on the GPU
       whose lists hold all records is drawn directly, at rest recorded as other draws. MG_1184 whole rotating, frame, ms:
       `auto` 1.73 → 1.16, `off` 2.11 → 1.21; zoomed to 0.3 2.63 → 2.60; at rest as before (CPU submission swings 0.4-0.6
       ms from run to run with the feature or without). Not for users yet.
     - The raster pass by stage, zoomed to 0.3 (a probe switching a stage off once the culler stopped dispatching: a
       pyramid missing the spheres otherwise lets more instances through; ms, at 1280×720 and at 2560×1440): the spheres
       (reading, transforming, bounds) 0.52, 0.77; the pixel loop's ray tests 0.56, 2.21; the atomics 0.48, 1.79; of all
       of these spheres above 64 pixels 0.05, 1.37. So pixels decide at 2560×1440, ray tests above atomics.
       Cheaper arithmetic per pixel, built and dropped: the ray not normalized, the sphere's term once per sphere, the depth
       from z without the matrix: the pass zoomed to 0.3 1.55 → 1.49 ms, at 2560×1440 4.73 → 4.65, whole there 2.57 → 2.31;
       and depth rounding changing which of two close spheres wins (1.4 % of the pixels by more than 24 levels whole at
       2560×1440, the raster's Dawn test failing). So the cost is not that arithmetic.
     - **Empty draws left out while moving** (done, `predictRecords`): from each unit's bounds, conservatively, per
       dispatch: outside the frustum, outside the level's range, or at a level whose instances all go to the raster.
       MG_1184, draws of 3,392: `auto` whole rotating 74, zoomed to 0.3 735; `off` 1,601, 2,310. The records change as
       the camera moves, so a chunk with such an object is drawn directly while its lists are not compact. Frame, ms,
       without the flag: whole rotating `auto` 1.73 → 1.15, `off` 2.11 → 1.58; zoomed to 0.3 `auto` 2.63 → 2.58, `off`
       4.52 → 4.58. With them multi-draw adds about nothing (1.28, 1.36, 2.79, 4.95; its runs of records many calls).
       `webgpu-gpu-culling`: the lists hold all instances the culler placed, at rest and at views moved to.
9. Later, by measurement: two-phase occlusion (item 9); depth passes and picking from the visibility buffer; tiles
   binned for the largest spheres (the workgroup's share done, step 8). Not cylinders: in practice not drawn in numbers as spheres are. The classification against clip
   objects (`clipSphereClass`) made general, culling cells, instances and primitives: [clip-culling.md](clip-culling.md).

## Step 4a: through the renderer

`sphere-raster-renderer/postsyn{-zoom}{-large}-{off,on4,on6,on8,only}` (`tests/src/benchmark/render.ts`): PetWorld as
`merging/postsyn-lod-pixel1` (levels by size on screen at 1 px, culled on the GPU with Hi-Z), the renderer's
`sphereRaster`; RTX 4070 Ti, Chrome 156 with `--enable-webgpu-developer-features`, `--gpu-timing`. Frame median / the
raster pass / the opaque pass (the impostors and the resolves), ms:

| case                 | off               | on6                | only              |
| -------------------- | ----------------- | ------------------ | ----------------- |
| whole, 1280×720      | 0.73 / - / 0.30   | 0.86 / 0.20 / 0.21 | 0.79 / 0.28 / 0.02 |
| whole, 2560×1440     | 2.24 / - / 1.78   | 3.49 / 0.77 / 2.22 | 1.30 / 0.77 / 0.08 |
| zoomed, 1280×720     | 16.1 / - / 15.9   | 10.6 / 1.91 / 8.35 | 2.25 / 1.90 / 0.04 |
| zoomed, 2560×1440    | 44.8 / - / 43.9   | 27.6 / 6.31 / 20.5 | 7.36 / 6.39 / 0.16 |

- `R` 4, 6 or 8 make no difference: at these levels almost every sphere is a few pixels across.
- **With the split per instance (step 4b)** `on` is `only`'s (on6, frame median / raster / opaque, ms): whole 0.90 /
  0.30 / 0.04, whole at 2560×1440 1.35 / 0.77 / 0.09, zoomed 2.30 / 1.91 / 0.05 (7× the impostors' 16.1), zoomed at
  2560×1440 7.34 / 6.33 / 0.16 (6.1× the impostors' 44.8). Whole at 1280×720 the frame is still slower than the
  impostors' (0.90 against 0.78): the CPU's, the opaque pass drawn directly (step 4c: with draw lists off 1.05, on6
  1.12, only 0.91, within noise; at 2560×1440 on6 1.03 against off 2.20). At 1 px about 14 % of the
  pixels zoomed differ by more than 24 levels from the impostors' (0.9 % whole): spheres under a pixel, the two paths
  pick another nearest sphere at a pixel's center; the images look the same.
- Levels at 3 px (the default), zoomed, with the level's minimum distance in the classification: off 2.80, on6 0.97,
  only 0.97 ms; at 2560×1440 15.1, 4.23, 4.17.
- Multi-sampled (`-ms`, four samples a frame, 1280×720, levels at 1 px): whole off 2.89, on6 2.38, only 1.60 ms;
  zoomed 61.9, 8.47, 8.39 (7.3×).
- As the Viewer adds them (`-viewer`: the target's default params, sized by `particle-size`, the elements by their
  physical radii), zoomed: off 3.37, on4 1.68, on6 1.47, only 1.38 ms. The split classifies by the member's
  `padding` (its largest size): one stale padding (the particle's radius, 86-391 Å, where the elements' are at most
  1.8 Å) left every instance to the impostors.
- `on` keeps the impostors' cost: their vertex stage runs for every sphere and collapses the raster's (whole at
  2560×1440 the opaque pass costs more than without the raster: the collapsed quads and a resolve covering the
  viewport). `only` draws no impostor: the raster's frame 6-7× faster zoomed. Hence the split per instance (step 4b).
- With the raster on the opaque pass drew directly (no draw lists): 619 draw calls for 18, the CPU's cost whole (draw
  lists since step 4c: 19; zoomed `-viewer` on6 1.42 ms, `only` 1.42).
- Images: zoomed, about 4.6 % of the pixels differ by more than 24 levels from the impostors', single pixels at
  sphere edges (where the raster's and the fragment stage's arithmetic decide differently, and SSAO and SMAA on them).
- Before the resolve cast the raster's ray (its point from `uInvProjection` instead) and took the tangent point on a
  grazing miss, pixels the raster claimed were discarded (white rings zoomed out) and the raster was 3× slower in the
  viewer: NaNs from the grazing rays are a suspect.

## Step 1: the prototype

The question it answers: for small opaque spheres, is a raster kernel and a resolve enough cheaper than the impostors'
opaque pass to be worth a second path. Nothing of the renderer changes; the prototype is a benchmark case with its own
code.

- **Where**: `tests/src/benchmark/sphere-raster.ts` (cases `sphere-raster/*`, on demand), on the WebGPU context's
  device directly (`WebGPUContext.device`), its own WGSL and buffers: the plumbing of steps 3 and 4 (kernel bindings of
  the data pages, the sphere table) is what the prototype is to justify, not to need. The device with the feature:
  `atomic-vec2u-min-max` added to `OptionalFeatures` (requested where the adapter has it; nothing uses it yet).
- **Input**: the spheres as the impostors draw them, from the data, instanced (PetWorld expanded is 2.39 G spheres,
  38 GB): each target's atoms once (`vec4`: center, radius as the representation sizes it: spacefill's `sizeFactor` is
  2 for coarse-grained structures, the PetWorld scene's uniform size 1), a transform per instance; the impostors'
  levels of detail mirrored (the atoms in level order, a level a prefix of them, its scale), the level picked per
  instance by its center's distance (the impostors pick per sphere), instances culled by the frustum.
- **Kernels**, per frame, timed with timestamp queries around each:
  - clear: the visibility buffer (`array<atomic<vec2u>>`, a pixel each) to 0;
  - raster: an invocation per sphere: the projected bounding rectangle (perspective, `sphereProjection`'s math), skip
    above `R` pixels (counted), per pixel center the ray's front hit, depth as reversed-Z float bits in the high word,
    the sphere's index + 1 in the low word, `atomicStoreMax`;
  - variants: with an early-out (a 32-bit `atomicMax` depth buffer beside it, read before the 64-bit write), and the
    two 32-bit passes (depth, then id where equal) as the fallback's cost.
- **Resolve**: a full-screen fragment pass reading the id per pixel, recomputing the hit, writing `frag_depth` and a
  Lambert shade (the impostors' lighting is not what is compared: the fragment cost is, by shading each pixel once).
- **Against**: the same camera's `pass opaque` GPU time of the impostors (`--gpu-timing`), with LOD as the case has it
  and off; its image as a check (the resolve's silhouette and depth agree with the impostors' within a pixel).
- **Cases**: PetWorld (`postsynstubby_2`, 3.68 M atoms) whole and zoomed (as `merging/postsyn-gpu-hiz-zoom-*`), the
  globule (`render/globule-300k-spacefill`), at 640×360, 1280×720 and 2560×1440; `R` 2, 4, 8 and unlimited.
- **Decision**: continue with step 2 if the 64-bit raster plus resolve, for the spheres at or below `R`, take at most
  about 60 % of the impostors' opaque pass for the same spheres in the PetWorld cases (a margin worth a second path
  and its upkeep); record the numbers either way, and the 32-bit passes' cost for the open question on a fallback.
- **Also measured**: the visibility buffer's memory (8 bytes a pixel: 29 MB at 2560×1440) and the clear's cost.

## Step 1: results

RTX 4070 Ti, Chrome 156 (D3D12, `--enable-webgpu-developer-features`), GPU ms per frame, medians of 15. Impostors:
`pass opaque`, culled on the GPU, the `quality` levels by distance. Raster: 64-bit raster plus resolve, `R` unlimited.

- **Two dispatch shapes**: direct (a workgroup per 1024 atoms of every instance, most empty under the levels) and
  compacted (a pass per instance culls it, picks its level and appends chunks of 64 spheres; an indirect dispatch over
  the chunks). The cull pass: 0.03-0.17 ms.

| PetWorld                 | spheres drawn | covered | impostors | direct | compacted | compacted / impostors |
| ------------------------ | ------------- | ------- | --------- | ------ | --------- | --------------------- |
| zoomed 640×360           | 40.7 M        | 88 %    | 8.11      | 10.74  | 4.87      | 60 %                  |
| zoomed 1280×720          | 40.6 M        | 88 %    | 13.14     | 13.42  | 5.91      | 45 %                  |
| zoomed 2560×1440         | 40.6 M        | 88 %    | 21.37     | 20.67  | 9.15      | 43 %                  |
| whole 640×360            | 12.2 M        | 6.5 %   | 2.53      | 6.19   | 2.78      | 110 %                 |
| whole 1280×720           | 12.2 M        | 6.5 %   | 2.45      | 6.46   | 2.95      | 120 %                 |
| whole 2560×1440          | 12.2 M        | 6.5 %   | 4.28      | 7.88   | 3.95      | 92 %                  |
| zoomed 1280×720, no LOD  | 1.92 G        | 88 %    | 377       | 136    | —         | 36 % (direct)         |

- **Gate**: passed zoomed in, failed whole (about 200 spheres per covered pixel: the atomics on the same pixels are the
  suspect). `R` changes little: almost every sphere is at most 2 pixels (+1) in radius.
- **32-bit two passes**: about twice the 64-bit raster; the early-out (a 32-bit depth read before the 64-bit write)
  never faster.
- **Globule** (300 k atoms, all above 8 pixels at 1280×720): impostors 4.04, raster 3-8 (unstable between runs): large
  spheres stay on the impostors, as designed.
- **Images**: the resolve's coverage matches the impostors' (the globule's silhouette to the pixel, PetWorld within
  1.6 % of the pixels, mostly dark pixels under the mask's threshold).
- **The impostors' occlusion**: their GPU cull already tests each cell against the previous frame's depth pyramid
  (`cell-cull.wgsl.ts`, the pyramid passed whenever objects are GPU culled); the `hiZ` prop adds only the CPU test
  (more draw calls, no gain). With the pyramid withheld the opaque pass takes 14.9 ms zoomed at 1280×720 (13.2 with
  it), 2.74 whole (2.46): it saves about 11 %; smaller cells (50, 20 Å) about 15 %. The ratios here are against
  impostors with it, as it was then. Since fixed (`merged-render-objects.md`, GPU occlusion step 7: the test rarely
  tested, a still frame stopped culling too early): the impostors' opaque pass 0.99 ms whole and 6.11 zoomed (0.71
  and 4.68 with 20 Å cells), 6.59 zoomed rotating. Against those the raster with exact bounds and no occlusion
  (1.00, 2.42) is about the same whole and 40 % zoomed; with two-phase occlusion (0.28, 0.85) 28 % and 14 %; rotating
  zoomed 1.35 against 6.59, 20 %. The gate still passes zoomed; whole, only with the raster's occlusion.

### Improvements

Compacted, the same chunks; with the cull and the resolve.

- **Exact pixel bounds**: the perspective sphere's screen bounds from its tangents, only the pixels whose centers it
  can cover (sub-pixel spheres test one pixel or none; the first bounds, the radius +1 pixel, tested about 25): the
  gain, about 3×.
- **A workgroup tile**: a chunk's depths reduced in a 32×32 shared tile at its instance's rectangle first, only each
  pixel's nearest of the chunk on the global 64-bit atomic: slower (the second loop and the barriers cost more than
  the contention saved); dropped.
- **Hi-Z**: instances culled against the farthest depth key per tile of a pyramid built from the visibility buffer
  (8 to 128 pixel tiles, the level where the instance spans at most 3 tiles a side). First single-phase against the
  previous frame (the table: a static camera only). Then two-phase, for a moving camera: phase 1 draws the instances
  visible in the previous frame (a bit per instance), the pyramid is built from their pixels, phase 2 tests every
  instance against it with the current camera, draws the newly visible and updates the bits. Checked against no
  occlusion pixel by pixel (sphere id and depth): no pixel differs, still or rotating. Per instance it culls far more
  than the impostors' per cell (78 % of the chunks zoomed in, against their 11 % of the time).

| PetWorld         | impostors | first | exact bounds | / impostors | + Hi-Z, single-phase | / impostors |
| ---------------- | --------- | ----- | ------------ | ----------- | ------ | ----------- |
| zoomed 640×360   | 8.13      | 4.87  | 1.72         | 21 %        | 0.49   | 6 %         |
| zoomed 1280×720  | 13.17     | 5.94  | 2.43         | 18 %        | 0.77   | 6 %         |
| zoomed 2560×1440 | 21.49     | 9.15  | 5.81         | 27 %        | 2.33   | 11 %        |
| whole 640×360    | 2.52      | 2.78  | 0.91         | 36 %        | 0.23   | 9 %         |
| whole 1280×720   | 2.45      | 2.94  | 0.98         | 40 %        | 0.24   | 10 %        |
| whole 2560×1440  | 4.30      | 3.95  | 1.48         | 34 %        | 0.40   | 9 %         |

Two-phase at 1280×720, the camera still or rotating 2° a frame (both sides; the second phase's raster is what the
motion adds):

| PetWorld          | impostors | exact bounds | / impostors | two-phase | / impostors | phase 1 raster | phase 2 raster |
| ----------------- | --------- | ------------ | ----------- | --------- | ----------- | -------------- | -------------- |
| zoomed, still     | 13.17     | 2.42         | 18 %        | 0.85      | 6 %         | 0.62           | 0.00           |
| zoomed, rotating  | 11.68     | 2.66         | 23 %        | 1.35      | 12 %        | 0.83           | 0.31           |
| whole, still      | 2.45      | 1.00         | 41 %        | 0.28      | 11 %        | 0.16           | 0.00           |
| whole, rotating   | 2.42      | 0.99         | 41 %        | 0.23      | 10 %        | 0.10           | 0.01           |

### Against the impostors as they are now

Re-measured after the impostors' occlusion fixes and their finer GPU cells (`merged-render-objects.md`, GPU occlusion
steps 7 and 8); the raster unchanged. 1280×720 unless noted, GPU ms, the raster with its cull and the resolve:

| PetWorld               | impostors | exact bounds | / impostors | two-phase | / impostors |
| ---------------------- | --------- | ------------ | ----------- | --------- | ----------- |
| zoomed 640×360         | 3.05      | 1.72         | 56 %        | 0.57      | 19 %        |
| zoomed                 | 5.14      | 2.42         | 47 %        | 0.83      | 16 %        |
| zoomed 2560×1440       | 8.56      | 5.79         | 68 %        | 2.41      | 28 %        |
| whole 640×360          | 0.76      | 0.91         | 120 %       | 0.27      | 36 %        |
| whole                  | 0.76      | 0.98         | 129 %       | 0.27      | 36 %        |
| whole 2560×1440        | 1.51      | 1.48         | 98 %        | 0.51      | 34 %        |
| whole, rotating 2°     | 0.81      | 0.98         | 121 %       | 0.22      | 27 %        |
| zoomed, rotating 2°    | 5.57      | 2.61         | 47 %        | 1.32      | 24 %        |

- **What is left of the margin is mostly occlusion**: the raster's spheres cost about the impostors' whole (exact
  bounds alone 98-129 %), half zoomed; with its occlusion it draws a fraction. Its occlusion is per instance (and its
  level per instance), the impostors' per cell (up to 4×4×4 of a 200 Å grid cell).
- **Per instance for the impostors too** (done, `merged-render-objects.md` GPU occlusion step 9): their copy kernel
  tests each instance of a drawn cell. Their opaque pass 0.71 ms whole, 4.64 zoomed, 5.10 zoomed rotating (from
  0.76, 5.17, 5.64): less than the raster's occlusion suggested. Against them the raster with exact bounds 2.42 zoomed
  (52 %), 0.98 whole (138 %); with two-phase occlusion 0.83 zoomed (18 %), 0.27 whole (38 %), 1.32 zoomed rotating
  (26 %). Both now cull per instance; what the raster has besides: its level picked per instance (the impostors draw a
  cell's instances at each level the instance overlaps, chosen per sphere in the vertex stage), a pyramid of the frame
  (two phases) rather than of the previous one, and cheaper spheres zoomed (exact bounds alone half of the
  impostors). Not finer pyramid tiles: the impostors' with tiles down to about 2.5×1.5 pixels draw the same
  (`merged-render-objects.md`, GPU occlusion step 10).

### Levels by size on screen

The cases above used the `quality` levels by distance; the particles' default is levels by size on screen
(`lodPixelSize` 3, the chain of `Spheres.LodPixelPresets`). `sphere-raster/postsyn{,-zoom}-px{3,1}-720p`: the impostors
with the levels and `lodPixelSize`, the raster with them as distances (a pixel's Å at a distance, 45°, 720 pixels
high; the level picked per instance, the impostors' per sphere). GPU ms, 1280×720:

| PetWorld      | spheres drawn | impostors | exact bounds | / impostors | two-phase | / impostors |
| ------------- | ------------- | --------- | ------------ | ----------- | --------- | ----------- |
| whole, px 3   | 2.6 M         | 0.30      | 0.98         | 327 %       | 0.20      | 67 %        |
| zoomed, px 3  | 7.8 M         | 2.39      | 2.51         | 105 %       | 0.74      | 31 %        |
| whole, px 1   | 2.6 M         | 0.30      | 0.98         | 327 %       | 0.20      | 67 %        |
| zoomed, px 1  | 115 M         | 15.67     | 4.09         | 26 %        | 1.68      | 11 %        |

- **At the default** (px 3) the spheres are a pixel or more across: the impostors are cheap, the raster without
  occlusion no faster (zoomed) or slower (whole); only with its occlusion (parked for the impostors as for it,
  `two-phase-culling.md`) ahead.
- **At high detail** (px 1, zoomed: 115 M spheres, mostly under a pixel) the impostors pay their cost per quad, the
  raster without occlusion takes a quarter of their time. Whole, px 1 draws as px 3 (the coarsest level either way).
- **So step 2 pays where detail is high** (`lodPixelSize` near 1, the `ultra` preset's 1.5): what the raster makes
  affordable is the detail that makes the impostors slow.

At 1.5 and at 2560×1440 (`-px{3,1.5,1}-{720p,1440p}`): the raster with exact bounds, no occlusion, against the
impostors' opaque pass, GPU ms (two-phase occlusion in parentheses):

| PetWorld     | px 3                      | px 1.5 (`ultra`)          | px 1                      |
| ------------ | ------------------------- | ------------------------- | ------------------------- |
| zoomed 720p  | 2.51 / 2.39: 105 % (0.74) | 2.75 / 9.47: 29 % (1.26)  | 4.09 / 15.67: 26 % (1.68) |
| zoomed 1440p | 6.62 / 13.93: 48 % (3.28) | 11.32 / 33.27: 34 % (4.55) | 17.33 / 44.01: 39 % (5.17) |
| whole 720p   | 0.98 / 0.30: 327 % (0.20) | 0.98 / 0.30: 327 % (0.20) | 0.98 / 0.30: 327 % (0.20) |
| whole 1440p  | 1.72 / 0.45: 382 % (0.51) | 1.57 / 1.18: 133 % (0.57) | 1.34 / 1.77: 76 % (0.62)  |

- **Where the spheres cover the screen** (zoomed): the raster without occlusion 26-48 % from `ultra` down, and at
  1440p already at the default (the impostors 13.9 ms there). The gate passes without occlusion.
- **Where few pixels take many spheres** (whole: 6.5 % of the pixels covered, about 40 spheres a covered pixel at
  720p): the raster slower, its cost the cull and the atomics on the same pixels; the impostors cheap there (their
  coarse levels). Ahead only at px 1, 1440p.
- **For step 2**: the raster not for every frame or object; a choice per object and frame by what it draws (spheres
  per covered pixel, or the impostors' time) is part of the design.

With occlusion, as the impostors had it then (one phase, per instance against the previous frame's pyramid; two since, `two-phase-culling.md`): the
raster's one phase tests each instance against the pyramid of its previous frame's visibility buffer with that frame's
camera (`cullKernel('one')`, `prevView`/`prevProj`); its two phases as above. The rotating cases (2° a frame) compare
one more step's frame with the same frame without occlusion. GPU ms against the impostors' opaque pass, the pixels
that differ (`pop`):

| PetWorld                    | impostors | exact bounds | one phase (pop)   | two phases |
| --------------------------- | --------- | ------------ | ----------------- | ---------- |
| zoomed 720p px 3            | 2.42      | 2.51 (104 %) | 0.73 (30 %)       | 0.82 (34 %) |
| zoomed 720p px 3, rotating  | 1.75      | 2.63 (150 %) | 0.53 (30 %), 375  | 0.62 (35 %), 0 |
| zoomed 1440p px 3           | 14.22     | 6.61 (46 %)  | 3.10 (22 %)       | 3.18 (22 %) |
| zoomed 1440p px 3, rotating | 8.89      | 6.96 (78 %)  | 1.79 (20 %), 492  | 2.13 (24 %), 0 |
| zoomed 720p px 1.5          | 9.63      | 2.74 (28 %)  | 1.20 (12 %)       | 1.29 (13 %) |
| zoomed 720p px 1.5, rot.    | 5.46      | 2.84 (52 %)  | 0.58 (11 %), 474  | 0.69 (13 %), 0 |
| zoomed 1440p px 1.5         | 33.92     | 11.53 (34 %) | 4.53 (13 %)       | 4.69 (14 %) |
| zoomed 1440p px 1.5, rot.   | 34.80     | 12.07 (35 %) | 3.26 (9 %), 455   | 3.67 (11 %), 0 |
| zoomed 720p px 1            | 16.00     | 4.11 (26 %)  | 1.61 (10 %)       | 1.70 (11 %) |
| zoomed 1440p px 1           | 44.57     | 17.48 (39 %) | 5.31 (12 %)       | 5.38 (12 %) |
| whole 720p px 3             | 0.30      | 0.99 (330 %) | 0.17 (57 %)       | 0.20 (67 %) |
| whole 720p px 3, rotating   | 0.32      | 1.00 (312 %) | 0.18 (56 %), 3    | 0.22 (69 %), 0 |
| whole 1440p px 3            | 0.46      | 1.72 (374 %) | 0.43 (93 %)       | 0.53 (115 %) |
| whole 1440p px 3, rotating  | 0.50      | 1.72 (344 %) | 0.41 (82 %), 10   | 0.52 (104 %), 0 |
| whole 720p px 1.5           | 0.30      | 1.00 (333 %) | 0.16 (53 %)       | 0.20 (67 %) |
| whole 1440p px 1.5          | 1.20      | 1.57 (131 %) | 0.47 (39 %)       | 0.56 (47 %) |
| whole 1440p px 1.5, rot.    | 1.09      | 1.57 (144 %) | 0.53 (49 %), 106  | 0.68 (62 %), 0 |
| whole 720p px 1             | 0.29      | 0.99 (341 %) | 0.16 (55 %)       | 0.20 (69 %) |
| whole 1440p px 1            | 1.78      | 1.34 (75 %)  | 0.53 (30 %)       | 0.62 (35 %) |

- **With one phase, like for like**: the gate passes in every case but whole 1440p at the default (82-93 %); zoomed
  9-30 %, whole mostly 39-57 %. The popping about the impostors' one phase (zoomed 375-492 pixels a 2° step, 251
  theirs at distance levels), not visible in normal use.
- **Two phases cheap here**: 3-20 % over one phase (the impostors' 38-78 %): what the second phase lets through behind
  uncovered gaps costs the raster little per sphere, and all of it is compute (no pass split). Exact in every case.
- **Step 2 with occlusion**: one phase at least, two worth considering; the choice per object and frame still for the
  whole views at the default.

## Open questions

- **The defaults `auto` and 24 px, reviewed and kept** on the RTX with the fixes of step 8. `auto` never behind `off` beyond noise:
  PetWorld and MG_1184 whole, zoomed and rotating, 1.1-2.1× ahead where it counts. Limit (frame, ms, 6 / 24 / 64 px):
  MG_1184 zoomed to 0.3 2.31 / 2.20 / 2.19, to 0.03 2.33 / 1.96 / 2.00, to 0.3 at 2560×1440 8.72 / 6.62 / 6.59; PetWorld
  whole 0.92 / 0.87 / 1.03, zoomed to 0.1 2.30 / 2.12 / 2.13, to 0.03 2.61 / 2.20 / 2.14. So 24: 6 loses zoomed, 64 gains
  nothing. The visual suite on all 7 profiles passes. Unmeasured since these fixes: the iGPU; never: other vendors' GPUs.
- `R` fixed, by GPU, or from the measured costs per sphere and per pixel? (With levels by size on screen almost every
  sphere is at most 4 pixels across; `R` matters close up, the globule's spheres above 8 cost the raster more.)
- The level per sphere in the kernel: an instance's chunks at each level it overlaps, as the culler copies it (about 4%
  duplicates at boundaries, `merged-render-objects.md`), or levels decided per sphere in one chunk set?
- Illumination: the tracing MRT outputs from the resolve too (color, normal, material, back depth from the back hit),
  or tracing keeps the impostors?
- Per-group transparency (`check_transparency` discards per fragment): eligible objects with a transparency layer
  split by group into the two paths, or the whole object stays on the impostors?
