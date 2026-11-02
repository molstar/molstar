# Plan: Culling by clip objects

## Status

- **Done**: what clip objects remove whole is culled, objects, cells and instances on the CPU and the GPU and
  primitives in the vertex stage, and primitives they do not cut skip the per-fragment clip tests. From the sphere
  raster's pixel clipping ([sphere-compute-raster.md](sphere-compute-raster.md) step 7), whose classification of
  spheres against the clip objects (`clipSphereClass`) it generalizes.

## Question

A clip object by pixel (`dClipVariant` 'pixel') cuts by fragment (`clip_pixel_wgsl`, `sphereImpostor`'s `dClipSphere`).
The culling (`renderer.ts` `cullObject`, the CPU culling of `MultiDrawBase`, the GPU culler `merged-culler.ts`) tests
the frustum, the levels and occlusion, not the clip objects: the clip objects only decide which solid-interior caps are
drawn (`getSolidInteriorCaps`, `Clip.canIntersectSphere`, a coarse may-intersect). So a half cut away by a plane is
culled, its vertices transformed and its fragments run and discarded; it also no longer occludes, so what is behind it
is drawn too. Measured with PetWorld zoomed and cut through its center (`clip/postsyn-zoom-*`, the clip objects at the
model's center; frame, ms): impostors culled on the GPU 5.12 (plane), 3.94 (ellipsoid, its inside removed), 3.97
(cone); culled on the CPU 8.24, 9.05, 8.90; the raster (`auto`) 2.87, 1.50, 1.48.

## Design

One test, three places: a bounding sphere (a cell's, an instance's, a primitive's) against the active clip objects,
from the signed distance of its center (`localSD`): **all** clipped (by one object, inside its clipped side by more than
the radius), **none** (outside every object's by more), or **some** (per fragment, as now). The distances divided by
their largest slope are bounds on the true ones (1; the cone's `|size.xy|`; the ellipsoid's distance already scaled by
its smallest axis), the radius scaled by the frame's largest scale: exact where it decides, conservative elsewhere.
With a clipping layer, flags per group exempt primitives of one instance apart: instances and cells are classified only
without a layer or with flags per instance (`uInstanceGranularity`), primitives always (their flag read as the vertex
stage reads it).

1. **CPU** (done, `@molstar/core` `util/clip.ts`): `Clip.classifySphere(objects, sphere, flag?)` → `SphereClass` none /
   some / all over the objects (one clipping every point decides), `classifySphereObject` per object,
   `getSignedDistance` the shaders' `localSD` on `frame` and `shape`; `canIntersectSphere` on it (some), its callers
   (the caps) unchanged: tighter than its inner and outer radii (exact for planes, cubes, cylinders, round spheres),
   and bounded for ellipsoids and cones (it was always true). `clip-classify.test.ts`: per shape (an ellipsoid and a
   wide cone among them), inverted or not, with a rotated frame and a translated and a scaled transform, 900 random
   spheres each, points within and on them sampled against the shaders' test: all and none agree, and at least 60 % of
   the spheres decided. The clip and cap visual cases on WebGL2 and WebGPU as before.
2. **Culling of objects and cells (CPU)** (done): an object whose bounding sphere its clip objects clip whole is not
   visible (`isVisible`, the draw table's rows: `DrawFact.ClipCull`, clip objects without a clipping layer, the clip
   cells among `DrawFactCells`, a clip change rewrites the row); in the CPU culling of instance grids (`cullSegment`,
   merged and unmerged) batches and cells clipped whole are dropped as the frustum drops them, a batch reaching no
   surface not testing its cells, the classes view-independent and cached per grid and clip version
   (`getGridClipClasses`: classified per frame they cost the CPU 2-9 ms, the culling runs per frame with Hi-Z), the
   clip cells among the cull cells (`getCullCells`, the cull caches follow a change); with a clipping layer nothing is
   culled. PetWorld zoomed, culled on the CPU (frame, ms): cut by a plane 8.24 → 4.47 (the opaque pass 7.34 → 3.83,
   the CPU's submission 7.11 → 4.13), by an ellipsoid 9.05 → 9.15, by a cone 8.90 → 8.63 (their pockets leave few
   cells clipped whole). `cull-segment.test.ts`: the draws as the cells give them without those clipped whole, after a
   clip change too, nothing culled with a clipping layer; the clip, cap and merged visual cases on WebGL2 and WebGPU.
3. **Culling on the GPU** (done, `cell-cull.wgsl.ts`): the cell and instance tests (`sphereLevels`) drop spheres
   clipped whole (`clipSphereClass`, moved from the raster into `clip_storage_wgsl`, the count its loop bound), in
   both kernels: the clip objects of an object without a clipping layer (`getCullClipCount`, either variant) after the
   segments in `segmentData`, now bound by the cell kernel too, the instance centers tested besides for instance
   clipping (`uClipInstance`, `syncClip` comparing it). PetWorld zoomed, culled on the GPU (frame, ms): cut by a plane
   5.12 → 1.76 (the opaque pass 4.82 → 1.46), by an ellipsoid 3.94 → 3.93, by a cone 3.97 → 3.95. `webgpu-gpu-culling`:
   members cut by a plane, the GPU drawing fewer instances than without (72 without the test, below 57.6) and the
   image the CPU culling draws; the clip, cap, merged and GPU-path visual cases on WebGL2 and WebGPU.
   - **The raster cut by a plane** (`clip/postsyn-zoom-plane-auto`): 1.25 ms (its pass 0.78), against the impostors'
     1.77; by an ellipsoid 1.48 (1.08), by a cone 1.45 (1.05). Its per-pixel clipping costs little: a 2.35 ms run
     (pass 2.0) was under-occlusion, depending on the run order: the GPU culler stopped culling a still camera after
     two rounds against pyramids of its view, but while work is pending (kernels compiling) the Hi-Z pass skips its
     builds and serves its last pyramid again, an early incomplete one, counted as a second round. A round now counts
     only a new build (`OcclusionPyramid.version`), the culler dispatching again until one arrives.
4. **Primitives in the vertex stage** (done for the impostors, spheres and cylinders; WebGL and WebGPU alike): with
   pixel clipping (spheres not by `dClipPrimitive`, which tests the center), each primitive's bounding sphere
   classified (`clipSphereClass`, now in `clip_core_wgsl` for programs and kernels, and in `common-clip.glsl.ts`):
   clipped whole, its vertices moved out of the clip volume; cut by none, a flat varying (`vClipSome`) that skips the
   fragment stage's clip tests (`clipTest` of the front, `clipExit` of the caps, the pixel test), whose results are
   then known (nothing clipped, `clipExit` returns its start); the resolve of the sphere raster tests every pixel as
   before. Exact: the clip cases' images as before on WebGL2 and WebGL1, `renderer/clipping/impostor-classes{,-interior}`
   (a grid of spheres and cylinders in all three classes, solid interiors or not) too. Frame (opaque pass), ms, before
   → after: WebGPU, the 300k globule at 2560×1440 by one clip object cutting nothing 18.2 → 15.4, by three 25.4 → 15.6
   (14.4 without clipping); PetWorld zoomed (`clip/postsyn-zoom-*-off`) by a plane 1.75 → 1.68, an ellipsoid 3.89 →
   3.70, a cone 3.92 → 3.65 (3.33 uncut: what is left is the drawing of what the cut reveals); culled on the CPU 5.31
   → 4.64, 9.07 → 8.28, 8.93 → 7.85. WebGL2: the globule by three 21.4 → 12.5 (by one 11.7 → 11.5), PetWorld's opaque
   pass by a plane 4.31 → 3.80, an ellipsoid 8.37 → 7.80, a cone 8.41 → 7.10.
   - Points and lines too: a point's fragments are all at its center (its model position the same at every corner),
     the point tested once in the vertex stage and not in the fragment stage; a line classified by the bounding sphere
     of its item's line (the start corners of a strip's line continue the line before, so not their own segment: every
     corner of the quad the same class). The images as before on WebGL2 and WebGL1 (`renderer/clipping/points-lines`,
     strips added). Meshes not (per triangle rarely pays).
5. **Measure** (done): the sphere cases with steps 3 and 4. Meshes and cartoon: `clip/assembly-60-{cartoon,gaussian}-plane{,-merged}`
   (60 copies of 7qpd, zoomed to 0.5, cut through the center by a plane; instanced and culled on the CPU, merged and
   culled on the GPU), the culling of what is clipped whole patched out → in (frame, opaque pass, ms): culled on the
   GPU, cartoon 0.88 → 0.75 (0.25 → 0.18), Gaussian surfaces 1.39 → 1.03 (1.14 → 0.78); culled on the CPU (and the
   merged ones on WebGL2) no change. The CPU culls by the instance grid's cells: at the default `cellSize` (200 Å)
   the copies (an assembly about 400 Å across) fall into 4 cells of 15, each straddling the plane (1 cell of the 28
   clipped whole); the GPU tests each instance. With cells of 40 Å 1.11 → 1.07 ms.
   - Later, by measurement: the CPU's instances classified one by one in the cells a clip surface crosses (view
     independent, cached as the cells' classes), each cell's draw split into the runs of instances not clipped whole
     (more draws where they alternate).

## Open questions

- Instance clipping (`dClipVariant` 'instance') on the CPU (settled): the grid's cells and batches classified by
  their instance centers (`getCullClipShrink`: their spheres less the invariant radius, the box of the centers expanded
  by it, the batches' by their cells' largest), as the vertex stage and the GPU culler test them; the object's sphere
  as before (conservative). `cull-segment.test.ts`: fewer instances than by the spheres, every one left out with its
  center clipped; PetWorld zoomed, cut and clipped by instance, culled on the CPU (`clip/postsyn-zoom-*-cpu-instance`):
  no measurable change (the opaque pass by a plane 3.14 → 3.13 ms on WebGPU; few cells reach a surface by their
  spheres and not by their centers, the vertex stage dropped those instances cheaply).
- The solid-interior caps (settled): culling what is clipped whole leaves the caps exact. An object clipped whole has no
  cap (its interior misses the cap object's surface, or another object clips the cap points). A mesh cap counts the
  faces behind the cap point by stencil parity, faces of the removed region beyond it too (a kept sphere's outside
  behind it); a culled cell holds whole instances, closed surfaces that cannot hold the cap point, each crossing the
  ray an even number of times. Impostor caps are per primitive. `renderer/clipping/solid-interior-cells-{plane,inverted-sphere}`
  (tori in small instance-grid cells, those outside a kept sphere culled): the images as without the culling on
  WebGL2 and WebGPU.
- Occlusion (Hi-Z) of what a cut reveals: unchanged, the pyramid is built from the frame's depth with the cut.
