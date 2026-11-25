# Plan: Every view of a frame culled on the GPU (stereo, XR)

## Status

- **Steps 1 to 3 done** (As built); step 4 (multiview on WebGL, view instancing on WebGPU) deferred: not now. Follows [particle-pack.md](../particle-pack.md) step 7's open item
  ("Every view culled on the GPU") and [two-phase-culling.md](two-phase-culling.md)'s "No occlusion" for the eyes.

## Why

- **One view a frame culled on the GPU**: a merged object's `prepareCull` dispatches once a frame, for its first view
  (`frame === gpuFrame`), and `cull` takes the GPU's lists only for that view's planes (`isGpuView`: the seven planes
  equal). Stereo draws its eyes one after another (`DrawPass.render`: `camera.left`, then `camera.right`) with the
  same frame, so the right eye's `cull` keeps the left eye's lists (`lastCullFrame`): drawn by the left eye's planes
  and levels, its outer side missing what the left eye does not see. XR goes the same way (`XRManager.update` fills
  the `StereoCamera`).
- **The second eye loses the sphere raster**: the raster drew only views at the target's origin (the resolve's
  fragments the visibility buffer's pixels), so the right eye, offset in its target, drew every sphere as an impostor.
- **No occlusion in either eye**: the pyramid is the main camera's; the eyes test nothing (`two-phase-culling.md`).
- **XR culls nothing on the GPU**: its camera is scaled (the world scale, `camera.scale`), and GPU culling runs only
  unscaled (`gpu-culling.ts` `updateScale`).
- Measured (the yeast cell, whole at 720p, WebGPU, RTX 4070 Ti; frame, ms): mono still 4.05 (raster 3.56, opaque
  0.20), rotating 17.4; stereo still 60.0, rotating 75.1. Stereo still by eye: left 16.0 (raster 15.7: no occlusion),
  right 43.8 (opaque 43.7: impostors, culled on the CPU); rotating left 30.3 (raster 15.8), right 44.1.

## Design

1. **One cull for the frame's views**: the views of a stereo frame culled together, once, with a frustum enclosing
   all of them (each plane of the first eye's moved out over the other eyes' corners: conservative, a few percent more
   drawn at the sides) and the main camera's plane for the levels. `cull` takes the GPU's lists for any view whose
   frustum lies inside the frustum culled (its corners within its planes) and whose camera plane is the one culled
   with, in place of the planes equal (`isGpuView`); the raster's input by the same test. Why not a dispatch an eye:
   the eyes overlap nearly whole, and a dispatch an eye doubles the compacted outputs (the yeast cell: its instance
   indices, 4 B an instance).
2. **The eyes' levels by the main camera's plane** (`renderer.update` of an eye camera: `uCameraPlane` and `LodView`
   the parent's): both eyes the same levels (on the CPU too), as the lists of step 1 choose them; no level popping
   between the eyes.
3. **Occlusion per eye** (after step 1, by measurement): a pyramid per eye from its own depth, the frame's one cull
   testing an instance against both (occluded only where occluded in both); or a dispatch an eye with its pyramid.
4. **XR's scaled camera**: GPU culling at `camera.scale` other than 1 (the planes and level distances in the scaled
   space, the kernels' model scale as the vertex stages have it).
5. **Both eyes in one draw where the device has it** (deferred). On WebGL multiview (`OVR_multiview2`, WebGL2): the passes into a
   two-layer texture array target, the vertex stages picking the eye's matrices by `gl_ViewID_OVR` (the view
   uniforms an array of two), each draw issued once for both eyes: half the draw calls and state changes of stereo,
   the CPU side of a frame (the yeast cell's second eye: 5,651 draws). Needs step 1 (both eyes the same lists). Where it
   exists: mostly mobile and headset browsers (the Quest browser; Chrome for Android), not desktop Chrome or Firefox,
   so for XR first; a pass reading the scene per eye (SSAO, outlines, postprocessing) runs per layer, or the multiview
   path draws without them (XR's `disablePostprocessing`). On WebGPU the
   [view instancing proposal](https://github.com/gpuweb/gpuweb/blob/main/proposals/view-instancing.md) (a draft, issue
   #4109): the optional feature `"view-instancing"`, a render pass's `viewCount` (dense, up to the limit
   `maxViewInstanceCount`, at least 2 with the feature), every attachment (color, depth, resolve) a `2d-array` view of
   `viewCount` layers, view `i` into layer `baseArrayLayer + i`; no view count on pipelines (a pipeline in such a pass
   draws every view); WGSL `enable view_instancing;` and `@builtin(view_index)` (`u32`, vertex and fragment stages).
   Native: Vulkan multiview (`viewMask`), D3D12 view instancing (`SV_ViewID`), Metal vertex amplification. The same
   shape as WebGL's: one path for both backends, the vertex stages' view uniforms an array indexed by the view. Its open
   questions touch us: render bundles (a compatibility bit, our draw lists), occlusion queries, indirect draws (the
   culler's records). Not in WebXR's WebGPU binding yet (an editor's draft, a pass a view). Raster and fragment work
   stay per view: it saves the CPU's draws, not the GPU's; the sphere raster (compute, a visibility buffer per view) is
   outside it.

## Steps

1. The joint cull and the containment test (items 1, 2): stereo on the yeast cell and PetWorld, still and rotating;
   images of each eye against the CPU's culling (the same pixels but at the sides, where more is drawn); the visual
   cases `renderer/camera/stereo*`.
2. XR's scale (item 4), tested with a scaled stereo camera (no headset in the tests).
3. Occlusion per eye (item 3), measured.
4. (Deferred) Multiview on WebGL and view instancing on WebGPU (item 5): the targets layered and the vertex stages'
   view arrays behind the extension or feature, the eyes' passes without it as now; tested where a browser has it
   (WebGL: Chrome for Android, the Quest browser; WebGPU: once a browser ships the feature), as the eyes drawn apart.

## As built

- **Step 1**: the eyes share one culling view (`RenderCamera.cullView`, the `StereoCamera`'s, updated with the eyes):
  the first enabled eye's frustum, each plane moved out over both eyes' clip corners, and its near plane through the
  eyes' midpoint as the camera plane of the levels (items 1, 2). `renderer.update` takes both from it where a camera has
  one, so both eyes hand `cull` and `getRasterInput` the same planes: `isGpuView`'s equality holds for both, no
  containment test needed, and the CPU's culling (the GPU's not ready, unmerged objects) culls by the same. The raster
  draws a view anywhere in its target (`uRasterOrigin`, the viewport's origin: the resolve's fragments and its view
  rays relative to it), not one starting outside it. Measured (the yeast cell as in Why): stereo still 32.3 ms (60.0),
  each eye's raster 15.8 (no occlusion: step 3); rotating 48.6 (75.1); mono unchanged (4.07, 17.5). Visual case
  `renderer/sphere-raster/merged-stereo` (both eyes by the raster, against the impostors culled per eye).
- **Step 2**: a scaled camera culls as an unscaled one, in the unscaled world: the renderer's `cullFrustum` and
  `cullPlane` (the frustum's and camera plane's constants over the scale: a point's distances to the scaled planes at
  its scaled position over the scale) go to `cull`, `prepareCull` and the raster's input test, so the level ranges
  compare unscaled distances; the raster split's near plane over the scale too. The shaders' scaled branch (the
  object's one level, no fade: `assign_lod_factor`, the raster's kernels and resolve) is gone, the distance unscaled
  already; `Renderable.cullSimple` (a level per object, the instances uncull) removed. The raster draws scaled views.
  The occlusion pyramid stays untested scaled (its spheres scaled; XR has none: its eyes, step 3). Visual case
  `renderer/sphere-raster/merged-grid-lod-scaled` (scale 0.05): the image of `merged-grid-lod` within one level at 16
  pixels; without the planes moved, 27,000 pixels differ. Not measured on a headset.
- **Step 3**: one pyramid of both eyes, not one each: the Hi-Z pass builds over the whole viewport, where the eyes lie
  side by side, and records the right eye's camera and where its half begins (`OcclusionPyramid.second`,
  `getStereoSplit`: both eyes enabled, tiling the viewport; else none, XR's views laid out another way). The test
  (`hiz_occlusion_wgsl`) places a sphere's rectangle in each eye's half, occluded only where occluded in both; texels of
  coarse levels straddling the halves hold the farther depth of both, conservative. The rectangles take the off-axis
  terms of the eyes' projections (`projectSphere` too: the raster resolve's scissor of an eye was shifted). One phase:
  the second phase builds its levels from the opaque depth so far, and the eyes draw one after the other, so what comes
  into view shows a frame late; a frame culled against levels of other views requests one more
  (`Canvas3DFrameRenderHost.requestFrame`), against levels of its own, so a still view is whole. Mono frames take no
  stereo pyramid nor stereo frames a mono one; the CPU test (one view's levels read back) stays off for the eyes; the
  culler's pyramid camera (`isPyramidCamera`) includes the second view. Measured (the yeast cell as in Why): stereo
  still 8.10 ms (32.3), each eye's raster 3.7 (15.8); rotating 15.3 (48.6), a rotating frame differing from the
  unoccluded one at 5,666 of 921,600 pixels (single spheres a frame late), a still one at none. Visual case
  `renderer/sphere-raster/merged-stereo-turned`; the catch-up frame not counted there (other pending work draws frames
  after a turn in the harness), its test a unit one (`hi-z-stereo.test.ts`).

## WebGL2 (measured, nothing stereo's left)

- The yeast cell at 720p, RTX 4070 Ti, WebGL2 (ms a frame, submit the CPU's): mono still 44.4 (submit 12, GPU wait
  32), stereo still 88.9 (submit 57, wait 32); mono rotating 129 (submit 128), stereo rotating 139 (submit 136).
  Timing mode (`--gpu-timing`) inflates the CPU side several times over: measure without it.
- **Still stereo is GPU-bound**: no sphere raster on WebGL, the impostors' vertex work per eye about a mono frame's
  (each eye 44 ms of GPU in timing mode, mono 54). The CPU profile has mono's calls (the draws, `uniform1i`) each
  about 2.3 times slower: the driver blocking behind the GPU, not more work. No occlusion: the CPU test (`hiZ.enabled`)
  is off by default and off for the eyes (one view's levels read back); a per-eye CPU test would need the eyes' pyramid
  read back.
- **Rotating is CPU-bound, stereo about mono** (one cull a frame for both eyes, step 1): 95% the CPU culling of the
  merged object's instances (`cullSegment` 44%, `emit` 18%, the frustum tests 25%, the camera distance 7%), about
  125 ms. A WebGL2 matter of large merged scenes, mono included: cells accepted or rejected whole, ranges emitted at
  once, culling less often while moving against a widened frustum; profile `cullSegment` and `emit` first.

## Open questions

- The eyes' shared levels (step 1, on the CPU too) change what an eye draws at a level boundary against each eye its
  own: the same level in both eyes, no popping between them; revisit if a stereo view shows a seam.
