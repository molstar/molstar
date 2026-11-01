# Plan: Occlusion per sphere in the impostors

## Status

- **Dropped** after step 1 (the prototype measured, then removed): the test never occludes the spheres that overdraw.
  Spheres hundreds of pixels across read pyramid texels of a few hundred pixels, whose farthest depth is the
  background through the gaps; the overdraw is of neighboring spheres at about the same depth, not hidden behind
  nearer surfaces. With the test from 0, 8 or 24 px (WebGPU, the sphere programs, `dSphereOcclusion`): every case
  within noise of without (60 copies merged zoomed 3.14 ms, at 2560×1440 12.3, PetWorld zoomed to 0.1 1.36, the 300k
  globule zoomed at 2560×1440 67.8, the 60 copies unmerged zoomed there 15.7). The same spheres dropped without the
  test (above 0 px): the globule 17.7 → 6.0 ms, the copies 3.14 → 1.17 (the path itself works).
- What the overdraw costs instead, by the raster (`sphere-compute-raster.md` step 8): every instance to the raster
  (the near plane test of the split left out, the raster skipping the spheres through it): PetWorld zoomed to 0.1
  3.24 → 2.01 ms (the opaque pass 1.36 → 0.06, the raster 1.48 → 1.54; none of its 8.4k spheres there crosses the
  near plane, its instances only by their spheres); the 60 copies 3.66 → 8.17 (the raster 0.24 → 7.68 for 15k spheres
  at 142× the view: slower than the impostors' 3.14).

- The question below blamed overdraw; much of the cost was spheres reaching behind the eye, each a quad over the view
  (fixed, `sphere-compute-raster.md` step 8: the zoomed cases below 30-35 % less).

## Question

Zoomed in, the impostors' cost is the overdraw of large spheres near the camera. The impostors write their depth, so
no early depth test: every fragment of every covered quad is shaded. Occlusion culls cells and instances (the CPU
culling, the GPU culler, against the pyramid of the previous frame), not the spheres of one instance: an instance near
the camera is visible, all its spheres drawn. Measured (`--raster-split-notes`, `auto` at 24 px; the opaque pass):

| case                               | spheres of the impostors in view | their discs   | opaque pass |
| ---------------------------------- | -------------------------------- | ------------- | ----------- |
| 60 copies of 7qpd merged, zoomed   | 16k (15k above 24 px)            | 142× the view | 3.13 ms     |
| the same at 2560×1440              | 26k                              | 144×          | 12.3 ms     |
| PetWorld (viewer) zoomed to 0.1    | 8.4k (789 above 24 px)           | 8.1×          | 1.35 ms     |
| PetWorld (viewer) zoomed to 0.03   | 5.0k (425)                       | 5.9×          | 1.24 ms     |

Unmerged objects are one instance (the raster's `auto` leaves them to the impostors): the 300k globule zoomed at
2560×1440 takes 68 ms, the 60 copies unmerged zoomed there 15.7 ms.

## Design

- The vertex stage of the sphere impostors tests the sphere against the pyramid the culler reads (`HiZPass.gpuPyramid`:
  the previous frame's, with its camera) as the culler does (`occluded` of `hiz_occlusion_wgsl`, its margins), and moves
  an occluded sphere out of [0, w]. Every corner computes the same test from the same inputs, so the quad goes whole.
  The popping is the culler's (one phase: what the previous frame lacked occludes less; accepted there, sphere raster
  item 9).
- Only spheres above a size on screen: the test costs a few texture loads per corner, small spheres' overdraw little
  (the threshold by measurement).
- One pyramid per frame: its texture and camera (`HiZSchema`) in the global bindings, not per object; the occlusion
  chunk on the block it reads (`uParams` in the kernels, the globals in the programs). Set where the renderer sets the
  culler's (`setOcclusionPyramid`): none in the passes of another camera or before a pyramid exists.
- WebGPU first; WebGL2 and cylinders after, by measurement.

## Steps

1. **Measure the potential** (WebGPU, the sphere programs only, a prototype): the opaque pass with the test, and per
   frame the impostor spheres it removes (a debug count), on `sphere-raster-scene/assembly-60-merged-zoom*` (`off`,
   `auto`), `sphere-raster-zoom/postsyn-viewer-zoom0.*`, `sphere-raster-scene/globule-300k-zoom*` and
   `sphere-raster-scene/assembly-60-zoom*` (unmerged); whole views for its cost where little is occluded. Thresholds
   0, 8, 24 px.
2. **The test**: global bindings of the pyramid, the shared chunk, the threshold. Dawn: spheres behind a pyramid at
   their depth not drawn, all without a pyramid; the visual suite's sphere cases (a static view, culled against
   pyramids of its own: the images as before).
3. **WebGL2 and cylinders**, by measurement (the pyramid texture is built on both backends).
4. **Measure** the cases of step 1 and the whole views; the iGPU's zoomed cases.

## Open questions

- The pyramid when nothing is culled on the GPU (unmerged objects): built where Hi-Z is enabled (`HiZPass.isActive`);
  with Hi-Z off, none, so no test.
- Picking and the depth passes draw the same view: the same test, or not (an occluded sphere is no pick target).
- The raster's large spheres (by the workgroup): the same test per sphere in its kernel (it tests instances now).
