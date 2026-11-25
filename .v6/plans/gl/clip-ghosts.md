# Plan: ghosts of clipped parts

## Status

- **Built** (steps 0 to 6): the inverted clip test, the ghost pass and its composite, the props, the opt-out, the CPU
  and GPU culling, the sphere raster, Prism, illumination, the measurements on PetWorld. How it was built, where it differs from the
  design below:
  - **No new render variant**: the ghost pass draws the `depth` variant (packed depth, as the transparent depth) with a
    global uniform (`uClipGhost`) the clip code reads: no further programs, no compile. The renderer's pass
    (`Renderer.renderClipGhost`, `PassDraws.clipGhost`, no caps) draws directly, not from draw lists, with a
    visibility of its own (`ghostRowVisibilities`, `isVisible`): an object with a ghost (`DrawFact.ClipGhost`, a
    check fact: clip objects, not opted out) whose bounding sphere the clip objects reach (classified all or some),
    Hi-Z as the visible pass (removed parts in front of the cut are not occluded; behind opaque surfaces never shown).
  - **Primitives whole by their centers**: spheres, cylinders, points and lines kept where their center is clipped, not
    tested per fragment (no inverse of the impostors' `clipExit`; the ghost is coarse); meshes and direct volumes by
    fragment (`clip_pixel`, the solid interior's test), instances by their centers (`clip_instance`) as before.
  - **Resolution**: the viewport and the pixel ratio scaled (sizes in pixels kept); `resolution` 0.5 by default.
  - **Composite** (`ClipGhostPass`, program `clip-ghost`): through a framebuffer with the image's color alone (the
    image's may hold the opaque depth it reads), premultiplied, before the helpers and the overlay; the ghost forces
    the image off the drawing buffer and the transparent depth where the scene has transparency, as the overlay. The
    scene's depth as the marking depth reads it (opaque, the front transparent layer under
    `marking.transparentThreshold`). Edges: a ghost texel next to an empty one or 5 % of its distance apart in view
    space. The shader manager links the `depth` variant while a ghost is drawn.
  - **CPU culling (step 3)**: `cullSegment`'s clip classes the other way round with `clipGhost` (batches and cells the
    clip objects do not reach dropped, those clipped whole kept untested, the cut ones' cells tested); the renderables'
    `cull` take the flag. The ghost pass culls with a frame of its own per frame (`ghostFrame`: an object culls once a
    frame, its lists cached by it) and a level view of its own (`ghostLodView`, a version of its own, its scales times
    `ClipGhostLevelScale`, 4: levels as if a pixel were 4 times larger, coarser). A merged object culled on the GPU
    (its lists for its frame, `gpuFrame`) is culled on the CPU for the ghost until step 4. Tests: `cull-segment.test.ts`
    (the ghost's lists as the cells give them), `renderer/clip-ghost/grid-cells`.
  - **Measured** (`clip/postsyn-zoom-plane-{off,off-ghost,cpu,cpu-ghost}`, PetWorld zoomed, cut by a plane; frame,
    ms): WebGPU 1.27 / 4.61 / 4.10 / 9.10, WebGL2 4.95 / 8.48 / 4.95 / 8.58. The ghost's GPU cost is small (the pass
    0.5 ms against the opaque pass's 1-3, the composite 0.01); the cost is its CPU culling (WebGPU 18 draw calls ->
    5,854, the submission 0.3 -> 4.3 ms). Its levels as if a pixel were 4 times larger halve it (scale 1: 9,666 draws,
    1.16 ms). WebGL2's opaque pass reads 2.9 -> 5.3 ms with the ghost: its timer queries count the GPU idling for the
    ghost's CPU culling before the next flush (culling without drawing reads the same, drawing nothing not; WebGPU's
    timestamps show 1.03 -> 1.15), not GPU work.
  - **A list set per pass**: each renderable (plain and merged) keeps the visible passes' lists and the ghost's apart
    (`createCullListSet`, the merged one's own): list, levels' uniforms, cull cache, last frame, draw version; drawn
    from the set of the last cull; draw lists record the visible set's. Neither culls over the other's lists.
  - **GPU culling (step 4)**: one more region of the merged culler, after the levels' and the raster's, at one level
    (the coarsest: `levels × instances` for all levels was too much memory, about 220 MB for a million instances and
    four), allocated from the first dispatch with the ghost (`prepareCull`'s `clipGhost`, objects with a ghost only;
    the layout key holds it). The cell kernel marks cells in view the clip objects reach (`GhostBit`, `ghostSphere`), the
    instance kernel copies those of its instances in the first phase (global counters, not the workgroups'); the args
    kernel fills its records with the others'. The ghost's list (`getGhostLists`): a record per segment, the coarsest
    level's uniforms at every distance (`uLod` 0 to the largest float); the merged renderable's ghost set takes it
    where the frame was culled on the GPU with the ghost. The ghost culls under the frame (no frame of its own: the list
    sets keep the caches apart). PetWorld zoomed, a plane, WebGPU: the frame with the ghost 4.61 -> 1.78 ms (1.29
    without), its pass 0.47 ms, the submission 4.3 -> 0.39 ms. Test: `webgpu-gpu-culling` (the ghost's region drawing
    what the CPU culling's ghost draws).
  - **The sphere raster (step 5)**: no code. The ghost pass draws the impostors (the `depth` variant, which leaves no
    sphere to the raster: only the color variants move the raster's out), the raster-only skip is the opaque passes',
    the ghost's lists hold every instance (the CPU culling's, the GPU's region). So the raster's objects' ghosts are
    impostors at the ghost's level, not a raster of their own. PetWorld zoomed, a plane, `auto`, WebGPU: 1.11 -> 1.61
    ms with the ghost (its pass 0.47). Cases: `renderer/sphere-raster/{clip-pixel,merged-clip-pixel}-ghost`.
  - **Prism (step 6)**: Ghost in the clipping section's Settings (`clipping.settings` `ghost`: the canvas's
    `clipGhost` on or off; its look, `ghostColor`, `ghostOpacity`, `ghostEdgeOpacity`, `ghostResolution`, rows shown
    while on, kept by the area while off; refused without a canvas), in the skill's guidance; the opt-out per
    representation is its `clipGhost` param (the representation's own settings).
  - **Illumination**: the ghost composited over the traced composite, before DOF, the overlay and marking (as over the
    visible passes' image); its depth drawn with the tracing input (the first iteration and a new multi-sample offset),
    not every iteration: it changes with the input only. The tracing input's `prepareCull` takes `clipGhost` (the
    merged objects' ghost lists culled on the GPU); the transparent depth is written for the composite where the scene
    is transparent, as for the overlay.
  - Cases: `renderer/clip-ghost/{pixel,points-lines,instance,grid-cells,opt-out,full-resolution,illumination,
    illumination-transparent}` (WebGL2 and WebGPU alike).

## Question

Clip objects remove what they cut away (`clip-culling.md`): culled whole where a cell, an instance or a primitive is
clipped whole, discarded per fragment where cut. What was removed leaves no trace: the context of a cut (the rest of a
cell around a cut-open organelle, the half of a capsid in front of a plane) is gone. A ghost, a faint fill and outline
of the removed parts, keeps it. It has to cost little on mesoscale models (PetWorld, MG_1184: millions of spheres,
merged objects culled on the GPU, the sphere raster), where the clipped part can be most of the model.

Drawing the removed parts again as transparent objects is what not to do: the whole clipped content at full detail and
resolution, through order-independent transparency, uncullable by the clip objects (it is what they cut) and the
transparency passes' cost on top. A ghost needs much less: one layer (the front-most removed surface), coarse (its
shape, not its atoms), at low resolution, flat.

## Design

A **ghost pass**: the scene's objects drawn once more with the **clip test inverted** (what a clip object removes is
kept, the rest discarded), at a **coarse level**, into a **low-resolution target** (depth and a mask, no shading), then
**composited** over the frame as a faint fill and its edges, as the marking pass composites its ghost edges.

- **Inverted clipping, the same classes.** The clip culling classifies a bounding sphere against the clip objects as
  all / none / some (`Clip.classifySphere`, `clipSphereClass`). The ghost takes all and some and drops none: the
  visible lists drop all. Both come from one classification: the CPU's cached classes per grid and clip version
  (`getGridClipClasses`) serve both lists; the GPU culler's cell and instance tests (`sphereLevels`) write a second
  set of draws from the same test (one kernel, two outputs, not a second dispatch). In the vertex stage the flipped
  class moves the primitives cut by nothing out (instead of those clipped whole), and the fragment test keeps what the
  clip test removes (a variant, `dClipGhost`, in `common-clip.glsl` and `clip_core_wgsl`). Clipping layers (flags per
  group or instance) apply as they do: a group exempt from clipping has no ghost.
- **Opt-out per object.** A geometry param (`clipGhost`, on by default) leaves an object out of the ghost pass: a
  membrane whose ghost would cover everything, solvent. Its rows and culler output skipped, as an invisible object's.
- **Coarse.** Spheres draw a level coarser than the visible pass picks (a level bias on `sphereLevels` and on the CPU
  levels, the coarsest at the extreme): few spheres a particle, scaled to cover it, as the far levels draw it. Other
  geometry (meshes, cylinders) has no levels: drawn as it is, its cost bounded by the culling and the resolution.
- **Low resolution, one layer.** A target of its own at half or a quarter of the frame (depth and a group or object id
  for edges), depth-tested within itself: the front-most removed surface only, nothing blended. The sphere raster
  draws its eligible objects there too, with the inverted pixel clip, at that resolution.
- **Composite.** A screen pass over the frame (after the opaque and transparent passes and the illuminated image,
  before the overlay): where the ghost is nearer than the scene, a flat tint at a low alpha and its edges (depth
  discontinuities and id changes, as the marking pass finds its edges), upsampled; behind it, nothing. No lighting, no
  post-processing on it, not picked.
  - **The scene's depth is the marking pass's** (`MarkingPass.getSceneDepth`): the opaque depth, and the front
    transparent layer where it is less transparent than a threshold (`marking.transparentThreshold`, 0.4). So a
    near-opaque surface in front hides the ghost, a faint one (a membrane) shows it through, binary at the threshold
    (hidden by the ghost's faintness). Behind an opaque surface a ghost is never shown: what is removed there was
    hidden before the cut too (showing it would be an x-ray, another feature). Not taken: the opaque depth only (the
    ghost over nearer transparent surfaces, in front of a membrane that is in front of it); composited into the opaque
    image before the transparent pass (transparent surfaces in front blend over it, but those behind it too, the ghost
    not in their depth test; and the illuminated image composes apart, its composite taking the ghost after).
- **Views.** Stereo and XR as the visible pass: the ghost per eye or view, its target as the visible pass's.
- **Cost.** What the ghost adds is bounded by what is clipped, at a coarse level and low resolution: no transparency
  pass, no full-detail redraw. Still frames are not redrawn (the frame's caching as now); a clip change redraws the
  ghost as it redraws the visible pass.

Not taken, kept as alternatives:

- **Bounding volumes** of the clipped instances (boxes or spheres per instance): the cheapest, but a mesoscale cut is
  thousands of overlapping blobs, no shape.
- **Proxy geometry per ingredient** (a low-poly hull or a coarse Gaussian surface, instanced): the best look, but
  memory per ingredient and a build per model; a later option for the coarse level where spheres read poorly.
- **The frame before the cut**, cached as an image: wrong as soon as the camera moves.

## Steps

0. **Cases and measurements.** Visual cases: a structure and a merged instance grid cut by a plane, an inverted
   ellipsoid and a cone, ghost on (fill and edges; behind a near-opaque and a faint transparent surface; an object
   opted out), on every profile. Benchmarks: PetWorld and
   MG_1184 cut through the center (`clip/postsyn-zoom-*` with the ghost on), frame and pass times against the ghost
   off: the target is a ghost costing a fraction of the visible pass.
1. **The inverted clip test** (`dClipGhost`): fragment and vertex stages, WebGL and WebGPU, one object's ghost drawn by
   a test harness; the images against the clipped and the unclipped object (ghost plus visible = unclipped).
2. **The ghost pass and composite** in `DrawPass`: the low-resolution target, the composite (tint, alpha, edges, the
   marking pass's scene depth), canvas props (`clipGhost`: off / on, color, alpha, resolution, level bias), the
   per-object opt-out, no culling changes yet (every object with clip objects drawn inverted).
3. **CPU culling**: the ghost lists from the cached clip classes (all and some), the level bias on the CPU levels.
4. **GPU culling**: the second output of the cell and instance kernels, its indirect draws per merged object (memory:
   one more set of draw arguments per merged object), the level bias.
5. **The sphere raster**: its eligible objects into the ghost target with the inverted pixel clip.
6. **Prism**: a toggle with the clipping controls.

## Open questions

- The default level bias and resolution: from the measurements of step 0 (a ghost of PetWorld's clipped half that
  still reads as its particles, at a fraction of the visible pass).
- `ClipGhostPass.isEnabled` scans the scene's primitives twice a frame with the option on (the draw pass's plan and the
  illumination's render; `compose` takes the latter's): short scans, stopped at the first object with a ghost. Once a
  frame would mean the value passed through the passes' plan; not worth it unmeasured.
