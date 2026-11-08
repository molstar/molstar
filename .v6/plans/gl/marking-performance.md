# Plan: marking performance

## Status

- **Current**: done: steps 0 to 3 and 5, step 4 dropped; what is left in the open questions. Benchmark cases
  `marking/postsyn-{whole,zoom}-{none,particle,hover,all}` (PetWorld as the viewer draws it, 1280x720, Hi-Z on; on
  demand, `postsynstubby_2-no_bonds.bcif` in `.benchmark/`) and `marking/assembly-60-hover-{redraw,rotating}`.
  - **Closing measurements**: `-hover-redraw-idle` (`RenderInput.idleMs`, 300 ms before each frame, left out of its
    time): zoomed, WebGPU / WebGL2, a frame 1.35 / 1.25 ms, the mask 0.02 / 0.06 ms (16 to 21 ms before step 5); what is
    left at idle clocks is the overlay over the screen (1.10 ms of WebGPU's 1.58). The assembly (60 instanced copies of
    7qpd, a render object per unit): the mask 0.02-0.03 -> 0.01 ms, WebGL2's hover redraw 0.54 -> 0.31 ms, the images
    the same as the culled lists'. The redraw's CPU (zoomed, profiled, about 0.5 / 0.25 ms a frame): `markedSegment`
    was the largest part, 0.19 ms, scanning the draw order through `aInstance`; it scans the counts by id with the
    inverse (`createInstancePositions`, kept per object): 0.08 ms. The rest is spread over the passes' encoding.
  - **Step 5 done**: the mask draws the marked instances alone (`Renderable.renderMarked`), unculled, in place of the
    culled lists where they are few (at most one marked instance per 4 instances, times the levels: `isMarkedFew`, else
    the culled lists). The markers keep per instance the count of its marked entries and the instances with any
    (`MarkerData.markerInstanceCounts`, `markedInstanceCount`): `Visual.mark` counts a mark's ranges before and after
    (what it marks), an action on everything by its status, or, removing only, the instances counted before. The lists
    (`SegmentedMdbList.markedSegment`) read the counts by the instance ids of the draw order (the grid reorders them),
    test each marked instance's sphere against the frustum and draw it at the levels its distance reaches, as
    `cullSegment` a cell; rebuilt where the markers (the merged object's `markerVersion`) or the view changed. With them
    the renderer skips the object's cull for the mask (`prepareMarked`): on WebGL the redraw culled every cell again
    with occlusion (a test that is never cached). No occlusion: the hidden marked instances are drawn, their ghost edges
    shown, which the culled lists left out where Hi-Z culled them (on PetWorld, a hovered particle behind others had
    none). Against step 3, WebGPU / WebGL2: zoomed, the hover redraw 1.44 -> 0.60 / 7.11 -> 0.51 ms (the mask 0.79 ->
    0.01 ms, 8 draw calls), multi-sampled 5.03 -> 1.0-1.6 / 11.8 -> 0.87 ms, rotating with a hover 2.41 -> 1.86 / 8.47
    -> 7.15 ms, one particle 2.47 -> 1.55 / 9.34 -> 7.2-8.7 ms; the whole model's redraw 0.66 -> 0.58 / 1.47 -> 0.50
    ms; everything highlighted unchanged (the culled lists). Images against step 3: one particle the same; the hovers
    differ by the ghost edges of a hovered particle hidden whole (60 to 90 px). Drawing every level of a marked instance
    instead (the shader keeps a sphere at its own) cost the whole model's mask 0.11 -> 0.44 ms: the finest level's
    vertices of a far particle.
  - **Step 3 done**: the marking depth target is filled by a quad from the scene's depth (`MarkingPass.getSceneDepth`,
    the `marking-depth` program): the opaque depth texture of the last eye (`DrawPass.getMarkingOpaqueDepth`, kept for
    the redraw) and, where unmarked transparent objects are less transparent than `marking.transparentThreshold` (0.4,
    as SSAO's), the front transparent layer where it is nearer. The mask's test allows a ten-thousandth of the distance
    and four steps of the packing's 24 bits: the scene's depth holds the marked fragments themselves, and far away a
    step of the packing is the larger (a ten-thousandth alone speckled PetWorld's far side). Multi-sampled, the scene's
    depth has the jitter of its last sample: a sample of another jitter reads it shifted by the difference, the farthest
    of the four texels around (a marked surface not hidden by itself half a pixel away, occluders' edges receding by up
    to a pixel); unshifted, the marked spheres showed ghost slices, one per sample's direction (the visual case
    `renderer/marking/multi-sample`: 201 px off the depth pass, 15017 unshifted). The depth pass stays with
    direct volumes, without an opaque depth texture (packed, path traced), with transparent objects to hide marking
    and no transparent layer drawn, and with transparent objects to dim and no transparent layer drawn (the dim pass
    takes the unmarked coverage from the marking depth: the default preset's water was not dimmed), and to dim read
    shifted (multi-sampled: the farthest of the four texels eroded the unmarked silhouettes, an undimmed rim). Zoomed, WebGPU / WebGL2: the hover redraw 3.6 -> 1.3 / 9.9 -> 8.1 ms (the depth
    0.01 ms), rotating with a hover 4.7 -> 2.7 / 16.5 -> 10.3 ms, the redraw multi-sampled 14.6 -> 5.0 / 31.5 -> 15.5
    ms; the whole model's redraw 1.55 -> 0.84 ms on WebGPU. PetWorld zoomed, everything highlighted: WebGL2's image the
    same, WebGPU's 290 single pixels at spheres' edges (the raster and the impostors pick another sphere). Visual cases:
    `renderer/marking/dim-occlusion` one pixel; `transparent-layer-*` the same with the threshold at 1 (their layer,
    alpha 0.5, hides the block).
  - **Step 2 done**: the redraw drew GPU-culled objects culled on the CPU, every cell, without the GPU's occlusion: a
    marker change updates the merged renderable, and its `update` dropped the frame's GPU culling; with the camera
    still no full frame culls again. It now keeps it where what the culling read is the same (no layout change, the
    culler's `dataVersion` unchanged: transforms, instances, clip objects). Zoomed, the hover redraw on WebGPU 10.1 ->
    3.6 ms (the depth 4.41 -> 1.88, the mask 1.84 -> 0.79, 24662 -> 738 draw calls), multi-sampled 29.8 -> 14.6 ms;
    the whole model 2.11 -> 1.55 ms. WebGL2 has no GPU culler: unchanged.
  - **Step 1 done**: `RenderInput.unforced` (frames drawn by what changed, no `requestDraw`); the cases
    `marking/postsyn-*-hover-redraw` (the camera still), `-none-rotating` and `-hover-rotating`, zoomed also with
    `-dim`, `-transparent` (an unmarked box, alpha 0.5, in the model's middle) and `-multi-sample`. Zoomed, the hover
    (frame ms, WebGPU / WebGL2):

    | | frame | marking depth | marking mask | draw calls |
    | --- | --- | --- | --- | --- |
    | nothing marked, forced | 1.50 / 11.7 | | | 19 / 12346 |
    | forced frames | 4.60 / 18.9 | 2.10 / 4.88 | 0.91 / 2.02 | 757 / 37008 |
    | redraw (camera still) | 10.1 / 14.9 | 4.41 / 4.32 | 1.84 / 1.90 | 24662 / 24662 |
    | nothing marked, rotating | 1.83 / 10.2 | | | 225 / 9408 |
    | rotating | 4.68 / 16.5 | 1.11 / 3.31 | 0.32 / 1.03 | 1399 / 28190 |
    | redraw, multi-sampled | 29.8 / 34.7 | 17.4 / 17.8 | 7.42 / 9.38 | 98634 |
    | rotating, multi-sampled | 17.1 / 47.6 | 4.57 / 13.6 | 1.27 / 4.23 | 5332 / 112692 |

    - **The redraw costs twice the marking of a full frame**: on WebGPU 10.1 ms against 4.6, the depth pass 4.41 ms
      against 2.10 and 24662 draw calls against 757; the whole model alike (2.26 against 1.18 ms). Step 2.
    - **Multi-sampling** pays the passes per sample (4 at `sampleLevel` 2): a redraw 29.8 ms on WebGPU.
    - **Dimming** costs nothing measurable here (its coverage pass 0.01 ms; `dimStrength` 0.5); the depth pass was
      drawn anyway. The transparent box changes nothing measurable.
    - WebGL2's frames are bound by the CPU (submit), noisy between runs (nothing marked, zoomed: 7.6 and 11.7 ms).
  - **Step 0 done**: a sphere's marker is the same at its four corners, so the vertex stage drops the spheres a marking
    pass does not draw (marked ones in the depth, unmarked ones in the mask) instead of shading them to be discarded.
    Zoomed, one particle highlighted: the mask 1.76 -> 0.92 ms on WebGPU, 4.14 -> 2.05 ms on WebGL2; the depth pass
    unchanged (it draws the unmarked spheres, nearly all).

## Where things stand

- **Passes per marking sample**: the depth (`markingDepth`: every object not fully marked, its unmarked fragments'
  depth and opacity, skipped when everything is marked or without depth test), the mask (`markingMask`: objects with
  any marker, its marked fragments tested against that depth), dimming's coverage (a quad over the depth and the mask)
  when dimming, the mask's antialiasing, the overlay into the layer.
- **Measured** (RTX 4070 Ti, 1280x720, after step 0), zoomed, one particle highlighted:

  | | WebGPU | WebGL2 |
  | --- | --- | --- |
  | frame, nothing marked | 1.51 ms | 7.58 ms |
  | frame, one particle | 4.57 ms | 14.3 ms |
  | marking depth | 2.10 ms | 4.98 ms |
  | marking mask | 0.92 ms | 2.05 ms |
  | draw calls | 25 | 12346 -> 37008 |

  The depth pass costs more than the scene: it draws every sphere as an impostor, the scene draws most with the sphere
  raster (1.09 ms). The mask is now its vertex stage: PetWorld's particles are merged by type, one highlighted particle
  marks an object holding nearly every sphere (`MarkerAny` is per object). WebGL2's CPU time doubles (7.3 -> 13.3 ms
  submit): both passes repeat the scene's per-cell draws.
- **What the app pays on a hover**: a marking-only change is redrawn over the last image (`MarkingPass.redraw`, from
  `Canvas3D`'s `canRedrawMarking`), the scene not drawn again: a hover costs the marking passes alone. The cases force
  full frames; they do not measure that path.

## Steps

1. **Cases for what is not measured**: marking is drawn in two ways, and each step is measured in both.
   - **Redraw**: the camera still, the hover moving (per frame another particle highlighted): a marking-only change,
     redrawn over the last image (`MarkingPass.redraw`). The frames not forced: `RenderInput` gets a flag that leaves out
     `requestDraw`, the frame drawn by the marking change alone. The marking passes are the whole frame.
   - **Full draw**: the camera rotating with the hover moving, against the camera rotating with nothing marked: the
     scene and marking drawn per frame, the scene's culling following the camera.
   - Both for the whole model and zoomed, and further: dimming on (`dimStrength` > 0); an unmarked transparent object
     over marked ones; multi-sampling. Why: the redraw is what users wait for on a hover, the full draw what marking
     adds while the view moves; dimming may cost more than the default (its coverage pass, the depth over the whole
     screen).
2. **The redraw culled as a full frame**: a merged renderable's `update` keeps the frame's GPU culling where what it
   culled by is the same (a marker change), so the redraw draws the lists of the last full frame. Why: the redraw is
   what a hover costs, and it cost more than drawing the scene and marking anew.
3. **The scene's depth instead of the marking depth pass**: the marking depth target filled from the scene's depth by a
   quad, the mask tested against it with a tolerance (marked fragments are in it themselves). Unmarked transparent
   objects hide marking where less transparent than `marking.transparentThreshold`, read from the front transparent
   layer. No depth pass in the redraw: the scene's depth is the last render's. Why: the depth pass drew the scene again,
   the largest part of marking.
4. **Dropped: WebGPU marking from the sphere raster** (the marked visible spheres from the visibility buffer). Why: it
   replaced the mask's vertex stage, which after step 5 costs 0.01 ms; the hidden marked spheres would still need
   their own draws, and the visibility buffer would have to live until marking's redraw.
5. **The mask over the marked instances alone**: per instance the count of its marked entries, kept with the markers;
   the mask draws the marked instances (in the frustum, at the levels they reach) where they are few against all, else
   the culled lists. Why it matters past the averages: a redraw after an idle pause runs at the
   GPU's idle clocks (an RTX 4070 Ti idles at 210 of 3120 MHz, memory 405 of 10501), and the mask, the vertex stage of
   every sphere of the object a marked particle is merged into, slows the most. One particle on PetWorld, zoomed,
   redrawn 300 ms apart: 0.96 ms back to back, 16 to 21 ms after the pauses (other passes 4 to 6 times slower).

## Open questions

- **Screen bounds of the marked instances** (a scissor for the overlay and the mask's antialiasing) were planned with
  step 5 and left out: the passes over the screen are what is left of a redraw, the overlay 1.10 ms at idle clocks
  (a few hundredths back to back); per instance or segment if done, never per group, the whole screen when dimming.
- **More uses of the per-instance counts**: updating visuals (an instance's markers, or other layers, only where it has
  any); not measured.

- **The pointer helper turns the scene-depth read off**: the debug and pointer helpers draw into the opaque depth
  texture (`DrawPass._renderHelpers`, `IlluminationPass.renderHelpers`), so marking draws its own depth pass while
  either is on (`DrawPass.markingOpaqueDepth` null), the cost step 3 removed. Fine for the debug helpers; the pointer
  is a user-facing option and must not cost that: draw it without writing that texture (the overlay or screen layer,
  or after the snapshot marking reads), then drop it from the check.
- Packed depth (no depth textures) keeps the depth pass: reading the packed opaque depth (`dPackedDepth`, as the
  overlay does) would drop it there too; its precision against the tolerance is untested.
