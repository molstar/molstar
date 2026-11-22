# Plan: GPU culling for the views of their own

## Status

- **Done** (steps 1-5): every merged object culled on the GPU culls a view of its own there too (illumination's
  thickness pass, a ray's pick, a frame's second camera, the overlay's and the screen layer's cameras), as a page
  store's object does since
  [particle-slots.md](../particle-slots.md) step 7, and draws by the GPU culling alone: no CPU lists or grids, the
  kernels waited for (but for the sphere raster alone, `rasterOnly`). WebGL culls on the CPU as before.

## Question

A merged object culled on the GPU (`MergedCuller`, not a page store's) draws the CPU culling's lists in a view whose
frame is not the one `prepareCull` culled: the thickness pass (`tracing.ts`, `createFrame()`) and a ray's pick
(`RayHelper`). For that its members' grids are built at the first such view (`calcInstanceGrid`, lazily,
[load-time.md](load-time.md): the yeast's 1.6 s) and kept in grid order, the CPU lists rebuilt a frame of that view. A
page store's object has no CPU culling and culls those views on the GPU (`Renderable.prepareCullView`); the others
could do the same.

The image pass is not among them: it draws through `DrawPass.plan`, whose `prepareCull` culls the image's own frame
and camera on the GPU (no pyramid, `setOcclusionPyramid(null)`). Measured in step 1: no CPU cull of an object with a
culler in an image.

What is left on the CPU after: WebGL alone. While the culling's kernels compile (`compute.prepareKernels`) an object
culled on the GPU waits for them, as a page store's does, and draws no CPU fallback.

## Design

- **`prepareCullView` for every culler**: the closure in `merged.ts` without the `culler?.paged` gate, and its `cull`
  drawing the GPU lists for that view's camera (`isGpuCamera`, as the page store's), not by `frame === gpuFrame`
  alone. Into the same buffers (as decided for the page store: no second set; the next `prepareCull` culls the frame's
  view again, its dispatch key). The ghost kept, no pyramid (the thickness pass's occlusion is off: a superset).
- **The kernels compiling: waited for, always**: an object with a culler draws nothing until its kernels are ready
  (as a page store's), in every view; no CPU lists for it in between. In part already: an object with pending grids
  draws none while they compile (`merged.ts` `cull`, `areKernelsReady`), and `prepare.ts` asks for them at precompile
  and at a prepared object; what is left is the object whose grids are not pending. Why: a CPU fallback for a few frames costs the
  grids (built, held, kept in grid order) for the object's whole life. The wait kept short where it shows: the kernels
  asked for when the context is made or the first culled object is added, not at its first frame.
- **The grids where still needed**: a context without compute alone (WebGL); an object culled on the GPU builds none
  (`calcInstanceGrid` not called, the data in its stored order).
- **WebGL**: unchanged (no compute; CPU culled). Out of scope.

## Steps

1. **Measured** (done; counters at the CPU cull of a merged object with a culler, removed after): the illumination
   cases on PetWorld (`postsyn`) and MG_1184 (`illumination-mesoscale/*`, WebGPU, Chrome, 1280×720, rotating):

   | Case                 | CPU culls        | CPU cull a frame | Grids, once | Heap (fixed → auto) | Frame (fixed → auto) |
   | -------------------- | ---------------- | ---------------- | ----------- | ------------------- | -------------------- |
   | MG_1184, whole       | one a frame      | 1.6 ms           | 184 ms      | 249 → 295 MB        | 1.8 → 3.9 ms         |
   | PetWorld, whole      | one a frame      | 0.15 ms          | 44 ms       | 205 → 205 MB        | 1.0 → 2.0 ms         |
   | PetWorld, zoomed 0.3 | one a frame      | 0.3 ms           | 44 ms       | 205 → 207 MB        | 3.0 → 26.8 ms        |
   | an image (any case)  | none             |                  |             |                     |                      |

   The CPU culls are the thickness pass's alone (fixed thickness, no back depth: none). On MG_1184 its CPU lists are
   also 3,700 more draw calls a frame (3,220 → 6,959), the frame's CPU submit 0.65 → 2.7 ms. The culling's kernels
   compile in 220-270 ms from a fresh page (one batch, `prepareKernels`). PetWorld zoomed is bound by the GPU (its wait
   25.6 ms): the thickness pass draws what is behind the front without occlusion, which GPU culling does not change;
   measured again in step 4, a plan of its own if it stands.
2. **The hook for every culler** (done; `merged.ts`): `prepareCullView` for every `MergedCuller`, and its `cull`
   drawing the GPU lists of a frame other than the culled one where of the same view (a page store's by its camera
   plane, the others' by all seven planes: no other frustum's lists drawn). Tests: `webgpu-gpu-culling.test.ts` (a
   frame and camera of its own culled on the GPU, the frame's image of that camera; the same camera's view after the
   frame's, and the frame's passes after it, the GPU lists; another view not culled for it the CPU's lists). Measured
   (as step 1):

   | Case (auto thickness) | Frame        | CPU submit   | GPU wait     | Draw calls    | Heap         |
   | --------------------- | ------------ | ------------ | ------------ | ------------- | ------------ |
   | MG_1184, whole        | 3.9 → 4.1 ms | 2.7 → 0.9 ms | 1.1 → 3.1 ms | 6,959 → 4,816 | 295 → 250 MB |
   | PetWorld, whole       | 2.0 → 2.1 ms | 0.7 → 0.8 ms | 1.3 → 1.3 ms | 482 → 482     |              |
   | PetWorld, zoomed 0.3  | 27 → 6.5 ms  | 1.0 → 0.9 ms | 26 → 5.6 ms  | 2,467 → 1,700 |              |

   Zoomed, the GPU's culling by instance draws far less than the CPU's by cell. MG_1184's GPU time the same
   (`--gpu-timing`: 3.16 → 3.19 ms; the cull 0.21 → 0.44 ms, the thickness view's dispatch, the thickness pass 1.92 →
   1.75 ms): bound by the GPU before too, its CPU time (submit 4.6 → 0.8 ms in timing mode) hid it.
3. **The wait for the kernels** (done; `merged.ts`): an object with a culler (but for the sphere raster alone,
   `rasterOnly`: culled on the CPU where the raster does not take it, by design) draws by the GPU culling alone: a view
   not culled on the GPU for it, the kernels compiling among them, draws none; no CPU lists, no grids built for them.
   A second camera of a frame (a stereo or XR eye) culled on the GPU for it too (the others' early return by all
   seven planes; jittered samples cull with the camera's, `setJitterCullView`). Tests: `webgpu-gpu-culling.test.ts`
   (nothing drawn while the context tells the kernels are not compiled, a frame's second camera, a view not culled
   for it none; the CPU culling's identity region and pending grids by `rasterOnly` objects). Members placed by
   matrices still build their grids at load (`createTransform`; poses' are pending where culled on the GPU,
   `createPosesTransform`): measured in step 4, kept.
4. **Measured** (done; as step 1, after step 3): the illumination cases as after step 2 (MG_1184 4.0 ms, CPU submit
   0.7 ms, 4,816 draw calls, 249 MB; PetWorld whole 2.1 ms, zoomed 6.4 ms, its GPU wait 5.7 ms: the zoomed frame's
   GPU time no plan of its own). The grids members placed by matrices build at load (timed in `updateTransformData`,
   removed after): 100 units of 1,000 copies merged and culled on the GPU 27 ms and 7.7 MB (of 145 ms building the
   representations), of 300 copies 12 ms and 2.2 MB, 60 copies of 7qpd merged 1.8 ms, the mesoscale models under
   1 ms (particle targets: pending). Kept: not worth a pending grid for matrices.
5. **The overlay and the screen layer** (done; `overlay.ts` `prepareLayerCulls`): their scenes' objects culled on the
   GPU for their cameras before the overlay pass opens (the overlay's camera of the drawn one, each screen scene's
   own), and a ray's pick for its camera; a pick of the drawn camera draws the draw's (the same planes). Without it an
   object of many instances culled on the GPU moved to the overlay drew none (no GPU view for it, no CPU lists). Test:
   `webgpu-gpu-culling.test.ts` (an object of 40 instances culled on the GPU in the overlay drawn and picked).
