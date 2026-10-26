# Plan: Draw submission on WebGPU

## Status

- Since the paragraph below: steps 2 (layer A) and 3 done (caps recorded, a draw block per cap,
  `setSolidInteriorCap`, not drawn directly after a bundle), 4 in part (theme layers in shared storage pages,
  [shared-draw-resources.md](shared-draw-resources.md) 3.1), 5 superseded (shared object groups with dynamic offsets
  instead of records by first instance), 6 done for every pass but illumination's inputs, data changes replayed
  without re-recording (`flushDataChanges`); `drawLists` no longer hidden, still off by default. 7 in part: GPU culling
  for merged objects ([merged-render-objects.md](merged-render-objects.md) step 6); for separate objects indirect draws
  measured a loss ([draw-table.md](draw-table.md) step 7). 8 open. Planning steps 1–5 done (5 by the draw table).
- `drawLists` as the default, measured and not done (`merging/units-1000-*`, WebGPU, RTX, the frame, two alternating
  runs each): unmerged still, hover and rotating 57%, 54% and 24% faster (cartoon hover 65%), but merged hover 19–21%
  slower (its marker changes flushed into the recordings), merged rotating 6–24% (view-dependent chunks walked), unmerged
  recolor 17–21% (464 chunks recorded again each frame), merged still within noise; `render/*` and PetWorld within a
  few percent. Merged objects (`renderMerged: auto` from 30 units) are the common case for many units, hovering them
  common. A default needs a policy: chunks recorded once stable for a few frames, drawn directly while their data or
  view changes every frame.
- That policy done (`drawChunk` in the renderer): a chunk whose inputs change in the submission after their last change
  is drawn directly where it was recorded again then (draws changing every frame, a recolor), drawn directly already,
  or holds few rows (`FewChunkRows`, 4: a merged object is one row, drawn as cheaply as walked and flushed); chunks of
  many separate objects keep being flushed and replayed (a hover sweeping across them, a rotation). It is recorded again
  after three unchanged submissions (`StableSubmissions`); `stats.drawLists.chunksDirect` counts the direct ones, the
  benchmark's `--draw-lists` turns lists on for every case. Measured (WebGPU, RTX, three tight interleaved rounds,
  lists against no lists): merged still, hover, recolor and rotating within noise (-20% to +5%; were +19-24%), unmerged
  still, hover and rotating 33-66% faster, unmerged recolor within noise (was +17-21%), `render/*` within noise but
  `boxes-1000` -72% and `units-1000-spacefill-unmerged` -76%.
- `drawLists` on by default since: also the small scenes (`render/7qpd-*`, -8% to +1%) and PetWorld
  (`merging/postsyn-lod-pixel3*`: still, zoomed, rotating, recolor, CPU culling; within ±2%) measured, three tight
  interleaved rounds. The benchmark's `--no-draw-lists` measures without; the visual suite's default pass draws
  directly, `--draw-lists` with lists. On an integrated GPU (AMD RDNA 2, `--force_low_power_gpu`, two rounds) lists are
  within noise or faster everywhere: `boxes-1000` -17%, unmerged units -3% to -7%, merged, small scenes and PetWorld
  -1% to +4% (frames there mostly GPU-bound).
- Paused after step 1: the upper bound is measured (see below) and shows the gain; steps 2–8 are not started and wait
  until draw submission becomes the priority again. Render bundles analysed against the draw path (see Render bundles:
  analysis): layer B is not needed for valid replays (per-object bind groups have fixed offsets), CPU culling can stay
  as part of the validity key; proposed order A, lean bundles (C), then B and D. Layer A done in a smaller form (a ring
  of the global blocks' own, see layer A); the measurements before C done (end of the analysis): recording costs no
  more than a direct frame, so bundles keyed by the visible set are never slower. C done for the opaque pass behind
  `RendererProps.drawLists` (see layer C): `boxes-1000` submit 6.3 to 0.7 ms; then every pass but illumination's
  inputs; data changes without re-recording next. The measurement code was removed afterwards (see below to redo it).
  Follows [webgpu-draw-overhead.md](webgpu-draw-overhead.md) (per-draw work of the render item reduced, uniform writers
  skip unchanged blocks) and replaces a plan for render bundles alone, which ran into the per-pass and per-draw state
  that bundles cannot hold. Builds on [storage-buffers.md](storage-buffers.md) (data arrays in storage buffers) and
  [indirect-draws.md](indirect-draws.md).

## Where the CPU time goes

`render/boxes-1000` (1000 separately drawn objects): about 6 ms of CPU per frame on WebGPU against about 4.5 ms on
WebGL2; `render/units-1000-spacefill-unmerged` (1000 render items of one representation, sharing the material): about
4.6 ms on WebGPU against 8 ms on WebGL2. Per object and pass, WebGPU encodes a pipeline, three bind groups (material,
object with its draw uniforms, and the global group when it changes), the vertex and index buffers and the draw; the
renderer loops over the objects (visibility, pass filters, face state, culling). Repeated state is already skipped
(`PassState`), so the cost is the state that differs per object: its own bind groups and buffers.

## Upper bound (step 1)

A benchmark case replayed the opaque pass of `render/boxes-1000` from a render bundle recorded once, through a temporary
hook of the WebGPU context: `getPassState` returned a pass state over a `GPURenderBundleEncoder` while recording
(viewport, scissor and stencil reference stayed on the pass), replay called `executeBundles` and invalidated the pass
state; the benchmark captured the renderer by wrapping `Renderer.create` and replaced `renderOpaque` (later frames ran
it on an empty group, so the pass's globals were still written and landed at the ring offset the bundle binds). The
frame equals the one drawn as usual in every pixel (checked after each run). Per frame, NVIDIA Lovelace, Chrome:

| `render/boxes-1000`          | frame  | submit (CPU) | wait (GPU) | draw calls         |
| ---------------------------- | ------ | ------------ | ---------- | ------------------ |
| WebGL2                       | 6.0 ms | 4.6 ms       | 1.4 ms     | (timing mode only) |
| WebGPU                       | 5.7 ms | 5.5 ms       | 0.2 ms     | 1008               |
| WebGPU, opaque pass replayed | 2.1 ms | 0.1 ms       | 2.0 ms     | 8                  |

The CPU time of the frame is almost all the opaque pass: replayed, the frame is bound by the GPU (2 ms, the 1000 draws
and SSAO), less than half of today's frame and of WebGL2's. Layers A–C are worth doing for scenes of many objects; the
gain needs the whole pass in a bundle (globals and per-object state without per-frame CPU work), which the replay
assumes: a static scene and camera, no culling, no change between frames.

## Goal

Frames whose CPU time does not grow with the number of objects, on WebGPU, without changing images or the WebGL backend:

1. per-object state reduced to an index, so that consecutive draws share pipeline and bind groups and each draw is about
   one call;
2. the draws of static passes replayed from render bundles;
3. culling and level of detail decided on the GPU for large instanced structures.

## Principles

- WebGL draws as today; the neutral API follows WebGPU (data arrays as `DataSpec`, object and material data as indexed
  records); WebGL lays them out as now (data textures, uniform blocks or classic uniforms).
- Every step is measured with the `render/*` benchmark (`--profile`, `--gpu-timing`) and checked with the visual suite
  (including `webgpu-core`, which requests no optional features) and the Viewer smoke test (device loss).
- WebGPU core only for every layer. Optional features (`indirect-first-instance`, multi-draw indirect, immediate data
  via the WGSL language feature `immediate_address_space`) may speed a layer up, each with a core fallback, unless a
  feature is supported virtually everywhere; none of them is today.

## Layers

### A. Per-pass and per-draw values out of the global uniforms

**Done in a smaller form** (step 2): the global blocks have a ring of their own (`WebGPUBindingContext.globalRing`).
Their snapshots of a submission take offsets in the order the frame changes the globals, so a frame of the same passes
binds the same offsets whatever the camera, time or jitter, and whatever other blocks (changed after a draw of the
submission) and compute kernels or in-pass clears write to the shared ring in between. That is what a recording needs:
it checks at replay that its pass still gets the global offset it recorded. Pipeline overrides for the per-pass values
and a persistent slot per global bindings are not needed (the ring exists because `writeBuffer` runs before the
submission's commands, so a slot rewritten in place could not hold the values of several passes of one frame);
solid-interior caps keep changing the globals between their draws and are drawn directly after a bundle. Checked: the
visual suite on WebGPU unchanged, `render/boxes-1000` submit 5.2 ms (5.5 before), a headless test that the global
offsets repeat over frames while the camera and colors change. The design below is kept for reference.


The global block holds values that change between passes or draws of one frame: `uRenderMask` (opaque or transparent
pass), `uSolidInteriorPass` and `uSolidInteriorClip` (solid interior caps, changed between the draws of one item),
`uPickType`, `uMarkingType`, `uDepthBack`. Each change means a new global snapshot (a new dynamic offset), so that a
pass cannot be recorded and caps cannot share state.

- Values with a few discrete states (render mask, cap pass, pick type, marking type, depth back) become WGSL `override`
  constants (core): part of the pipeline, which already differs per pass by its render state; more pipelines, compiled
  ahead by `linkPrograms` and cached.
- The cap's clip index (-1 for the near plane, up to `Clip.MaxObjectCount` clip objects) changes between the draws of
  one item: one draw-uniform slot per clip index in the item's uniform buffer (where its draw block already lives), so
  the cap draws differ in their dynamic offset only, which a bundle records; no pipeline per index.
- What remains global changes per frame (camera, light, fog, time): one persistent slot per global bindings, rewritten
  in place with `writeBuffer` instead of a ring snapshot.

### B. Indexed per-object data, shared bind groups

Measured upper bound (`draw-overhead/*`, see [draw-table.md](draw-table.md) step 6): 10k draws replayed from a render
bundle cost 37 ms per frame binding as Mol* does, about 1 ms with shared bindings; per distinct resource a draw binds
(geometry buffers 24 ms, bind groups 10 ms of the 37), not per command. Planned in
[shared-draw-resources.md](shared-draw-resources.md): geometry sub-allocated in shared buffers first, then the data
arrays in storage arenas with shared object groups; object records by `firstInstance` not needed for submission.

- **Object records**: the object block (and its draw uniforms) becomes a record in a storage buffer shared by all render
  items of a context (`array<ObjectRecord>`). The object index comes with the instance data: each draw's `firstInstance`
  (core for direct draws) selects the item's first entry in a shared instance section, whose entries hold the instance's
  transform, ids and object index. Allocation in a slot arena with a free list; records rewritten with `writeBuffer`
  when they change (the persistent slots of today, in one buffer).
- **Material records**: the same for materials (already shared per material id), one storage array of material records,
  indexed per draw.
- **Data arrays**: the per-item data textures (color, size, marker, overpaint, transparency, emissive, substance,
  clipping, wiggle, palette, sphere positions) become sections of shared data arenas, the
  [storage-buffers.md](storage-buffers.md) design shared per context instead of per render object; the record holds the
  arena page and the section offsets.
- **Pages, not one buffer**: theme data can be large (group-instance data of big assemblies reach the sizes of today's
  largest data textures, up to 1 GiB at 16384² texels), while core guarantees storage bindings of 128 MiB and buffers of
  256 MiB (the adapter's maxima are requested and often larger). The arenas are pages sized to the device limits; items
  are assigned to a page with room, one bind group per page and program, and the renderer orders draws by page so
  consecutive draws keep sharing the group. Arrays larger than a page get a buffer of their own with windowed bindings
  ([storage-buffers.md](storage-buffers.md) step 5) or stay data textures, with per-item bind groups as today.
- **Geometry**: vertex pulling from shared storage (spheres, cylinders, points, lines and text pull already; meshes
  follow), instance data (transforms, instance ids, segments) as sections with offsets in the record.
- Result: per program and page one bind group for everything (global, records, arenas); a draw sets its dynamic offsets
  and draws with its `firstInstance`. Textures that stay textures (grids, images, fonts, GPU-produced surfaces) keep
  per-item bind groups for those programs only.
- Page growth or a new page adds or rebuilds bind groups (it invalidates recordings, layer C).

### C. Render bundles

**Opaque pass done** (step 6, first part), behind `RendererProps.drawLists` (hidden, off by default):

- **Neutral**: `GraphicsContext.createDrawList?()` (WebGPU; WebGL has none and draws as before): `replay(globals)`,
  `record(globals, draw)`, `reset()`. The WebGPU list keeps a render bundle per target layout (`PassTarget.key`, at
  most 4: canvas and image passes), replayed while the pass binds the global slot and global textures it was recorded
  with (`WebGPUGlobalBindings.getState`, which takes the slot as a draw would, keeping later passes' offsets).
- **Renderer** (`renderPassFromDrawList`): the pass is planned as `renderPass` draws it (`isVisible`, `cullObject`,
  `getObjectState`): per draw the renderable, its render state, `programVersion` (transparency changes programs without
  changing values) and the version of its values (one `ValueCell` observer per values object, which visuals keep across
  render objects; volatile cells compared by version), plus the plan's key from `Renderable.appendDrawKey` (since
  `drawVersion`, see [draw-table.md](draw-table.md); the instance and lod ranges drawn and the level uniforms' values,
  so culled instanced items and merged items replay while their ranges hold). The same plan replays, otherwise the list
  resets (a changed plan) or records the current target; recording costs what drawing does. One list per group and time
  the group is drawn between two submissions (stereo eyes, sub-renders), see below. Objects drawing solid interior caps
  in the frame (the near plane or a clip object cuts them) keep the pass direct: the caps write the cap values into the
  object's draw block, which a replayed draw of the same submission would read.
- **Checked**: `pnpm test:visual --profiles webgpu --draw-lists` renders every case three times (compiling, recording,
  replaying) and compares the replayed frame with `webgl2` drawing directly: 145 pass, 116 of them replayed (the others:
  illumination, scenes without opaque objects, solid interiors with caps drawn); without draw lists all profiles
  unchanged. `stats.drawLists` counts replays and recordings.
- **Measured** (`render-bundles/*drawlists*`, RTX 4070 Ti): submit `boxes-1000` 6.3 to 0.7 ms (replayed while still or
  rotating), `units-1000-unmerged` 4.4 to 1.3 ms (planning culls the spheres per frame), the frames 4.6–6.5 to about 3
  ms (GPU-bound); zoomed and rotating 2.8 ms (records every frame, as direct with culling).

**Pick, depth, marking and opaque emissive passes done** (step 6, second part): the renderer's `drawPass` replays or
records every pass of these kinds, drawn directly where it can't be recorded (`planPass`: caps); the lists are per group
and pass, one per time the pass is drawn between two submissions (the pick types, stereo eyes, marking samples): the
n-th draw after a submission binds the n-th slot of the global ring each time (`GraphicsContext.submitCount`; counted
per frame token before, which marking-only redraws keep while the camera rests, so they recorded new lists every redraw
and never replayed), and per pass, so that submissions drawing other passes (a pick, a marking redraw) don't shift a
pass's lists. Checked: `--draw-lists` 146 pass, 27 cases drawn without replay (29 before); `webgpu-draw-lists`
(headless) picks the same objects, instances and groups replayed as drawn directly.
Left direct: illumination's inputs (tracing, back depth: drawn on the first iteration only).

**Transparent passes done** (step 6, third part): blended (with and without a transparent background), WBOIT, depth
peeling (a list per peel, each recorded with its layer's textures in the global set), volumes and transparent emissive
through `drawPass`. Checked: `--draw-lists` 146 pass, 15 cases drawn without replay: illumination (its inputs on the
first iteration only), solid interiors drawing caps, and a few (`geometry/merged/spacefill-mesh`, the GPU Gaussian
surfaces at low resolution, `postprocessing/outline/transparent`) whose captured frame records again, to look into.

The few that recorded again drew solid interior caps every frame: the near plane, fitted to the front of the scene's
sphere, touched their bounding sphere (an object spanning the scene, a merged item's sphere over its members'); the
near test now keeps a margin and a merged item's member spheres (`getSolidInteriorCaps`), and with `--draw-lists` only
the illumination cases draw without replay. `postprocessing/outline/transparent` replays (listed once while its pipelines compiled).

**Caps toward recording**: the cap pass and clip index are a draw state (`GlobalBindings.setSolidInteriorCap`), not a
change of the globals: on WebGPU each item binds a draw block per cap (and per multi-draw level), constant, so its
other draws keep their slot and the globals their offset; WebGL sets the uniforms as before. Images unchanged on all
profiles. Clears in recordings are draws of the bundle (`ClearDraw.drawInBundle`): under the pass's viewport (equal to the
pass's clear where the scissor lies inside it), a stencil clear to 0 with the `zero` operation (no reference, which is
the pass's), the values in a constant buffer per values; other stencil values or a scissor beyond the viewport still
fail the recording (`webgpu-bundle-clear`). **Caps recorded**: a recording draws the pass as `renderPass` does (caps and their
clears included); the plan's key holds which caps each cap-drawing object draws (`getSolidInteriorCaps`, shared with
drawing: the near plane and the clip objects cutting it, camera dependent) with its values' version, so a cut object
replays while its caps stay. A uniform arena growing (cap draw blocks are created on first use) replaces its buffer:
`UniformRing.arenaReplacements` fails a recording it happens in and drops recordings made before it (they may bind the
released buffer). Checked: `--draw-lists` 146 pass, the solid interior cases replayed; drawn without replay only
illumination (and cases whose pipelines compile while recording under load). Next: stencil bit planes per cap to avoid
most clears.

**Data changes replay** (markers on hover, colors, uniforms): a replay skips the render items' `render`, which writes
uniform slots and vertex data. `GraphicsContext.captureDraws` (WebGPU) runs draws against a capturing pass state:
their data written as by drawing, nothing encoded, what they bind returned as text (pipelines, bind groups with
offsets, buffers, draw arguments). The renderer captures each planned draw after recording; a plan that differs only
in values (`comparePlans` 'data') captures the changed draws again and replays where each binds what it bound
(`stats.drawLists.flushed`), records again where one binds otherwise (a texture or buffer replaced, a count changed, a
slot moved to the ring). Texture data (`update`, outside drawing) is written in place and seen by replays anyway.
Checked: `webgpu-draw-lists` toggles a highlight four times, replayed and equal to the direct frames (recorded six times
more without it); `--draw-lists` 146 pass. Recording captures every draw once more (its JS, no encoding).

**Planning cost** (next): a replayed frame of `render-bundles/units-1356-gaussian-drawlists` submits in 0.9 ms (7.3 ms
direct), most of it planning, redone identically on frames where nothing changed. Profiled (setup included): of
`renderPassFromDrawList` without recording, `planPass` ~87%; within it `isVisible` ~41% (frustum and lod tests,
`getAlpha`, `checkOpaque`, vector math, per object and pass), the cap decision ~12% (solid interiors are on by default),
`cullObject` ~7%, watches, versions and object states ~12%, building the plan (a draw object and key arrays per frame)
~9%; comparing and replaying ~5%; garbage collection notable. Steps:

1. Done: a plan keeps what it read (`PlanInputs`): the group's objects with their states (disposed, visible,
   pickable, color only) and program versions, an epoch of any change of their values (one observer on every watched
   cell) and the versions of volatile cells, the model view and projection, the model scale, the near plane and
   whether caps are drawn; unchanged, the list replays without planning (`stats.drawLists.unplanned`). Not with an
   occlusion test (Hi-Z changes between frames): `HiZPass.occlusionTest` is null while Hi-Z is off, it was a test
   answering false. `units-1356-gaussian-drawlists` submit 0.92 to 0.17 ms (rotating unchanged, 1.05 ms: planned);
   `webgpu-draw-lists` checks hiding, moving the camera and showing again against direct frames (failing without the
   camera and state inputs). Value changes still plan in full (then flush).
2. Dropped: visibility and caps cached per object for a view (a weak map per object, shared by the passes of a frame)
   saved ~0.03 ms; the lookup costs about what it saves.
3. Done: planning without allocating, no gain measured: two plans per list, swapped, draws, keys (`DrawKey`, written
   at a cursor: an array cleared with `length = 0` and pushed to again reallocates) and inputs written in place. Kept
   for culled objects, whose keys are long. Also done: the cap key only for objects drawing caps; an object's
   visibility once per pass, by group index, its caps reusing it (`highlighted-drawlists-rotating` planning 1.18 to
   0.99 ms, submit 1.49; plain rotating 0.76 to 0.56 ms, submit 0.83). Planning is per-object call overhead (checks,
   visibility, culling, states, about a third each for visibility and caps before), not allocation.
4. Done: the objects a pass takes (`PassMembers`: selected and checked with their draws and watches, the caps inline
   and after) per group and pass, evaluated again when the objects, their states, any value (the epoch) or whether caps
   are drawn change; the objects and states compared per plan (a frame is reused while nothing renders anew, not a
   bound). Not evaluated again while rotating; planning 0.99 to 0.85 ms (submit 1.31) highlighted, 0.56 to 0.53 ms
   (submit 0.74) plain. Replaced since by the draw table's pass predicates ([draw-table.md](draw-table.md) step 3).
5. Not done, replaced by [draw-table.md](draw-table.md): objects as flat rows maintained on change, per-view loops and
   pass masks instead of per-object planning (the steps above trim a visit, every pass still visits every object).

Measured with the `units-1356` cases (still and rotating, direct and with draw lists, also all highlighted: the marking
passes per frame); checked by `--draw-lists` and `webgpu-draw-lists`. The risk is completeness: an input not tracked
replays a stale pass.

Then: replaying zoomed views unculled by the visible fraction; then the default. The design below is kept for
reference.


With A and B, a pass's draws hold no per-frame state: they bind fixed buffers and differ only in dynamic offsets and
draw arguments. A pass is recorded into a `GPURenderBundle` and replayed with `executeBundles` while it stays valid:

- **Validity key**: pass kind and target formats, the scene version (add, remove, commit, visibility, the pass-filter
  inputs alpha, transparency average, x-ray, double-sided, flagged by `ValueCell` observers) and a draw version per
  render item (program or pipeline changes, buffer recreation, draw and instance counts). Data changes (colors, markers,
  uniforms) do not invalidate: they are buffer writes.
- **Neutral API**: `createDrawList()` with `isValid(key)`, `record(key, draw)`, `execute()`; WebGL draws in `record` and
  is never valid, so the renderer has one code path.
- **Culling**: bundled passes do not cull per object on the CPU (layer D culls on the GPU; until then the GPU discards
  off-screen geometry, and items that cull instances or switch levels of detail are drawn directly after the bundle).
- **Recording cost**: about one direct frame; scenes that change every frame (animation, trajectories) stop recording
  after repeated invalidations.

### D. GPU-driven culling and level of detail

A compute pass per frame tests instances or cells against the frustum, the previous frame's Hi-Z levels and the level of
detail distances, and writes the draw arguments and compacted instance lists (indirect draws, as for GPU surfaces). Core
only: the compacted lists hold global instance indices, which the vertex shader reads through the list's base in the
object record, so indirect draws keep `firstInstance` 0 and need no `indirect-first-instance`;

- removes the CPU culling (`isVisible`, `r.cull` rebuilding the multi-draw lists of large instanced structures), which
  grows with the instances of assemblies;
- keeps bundles valid while the camera moves (they record `drawIndirect` with fixed buffers);
- one indirect draw per object and level of detail (cheap inside a bundle); multi-draw indirect, where available, is an
  optional acceleration with this as the fallback;
- the CPU culling of today remains the fallback where compute is unavailable (WebGL) and for scenes that opt out.

## Render bundles: analysis

What a bundle needs, checked against the draw path of today (`gl/renderer.ts`, `gl/webgpu/render-item.ts`,
`global-bindings.ts`), and the options with their pros and cons.

### What a recorded pass holds, and what breaks it today

A `GPURenderBundle` holds pipelines, bind groups with their dynamic offsets, vertex and index buffers and draw calls,
for one target layout (color formats, depth-stencil format, sample count). `executeBundles` keeps the pass's viewport,
scissor, stencil reference and blend constant and resets everything else (so `PassState.invalidate` after it). Buffer
contents are read when the GPU runs it: writes between replays are seen. So a recording stays valid while the set of
draws, their pipelines, bind groups, offsets and buffers stay the same; data may change.

| per-frame state today                                                                 | breaks replay because                                       | fix                                                                                       |
| -------------------------------------------------------------------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| global block as a ring snapshot: each change (camera, time, jitter) a new offset       | the bundle records the offset of its frame                  | a persistent slot per global bindings, rewritten in place (layer A)                       |
| per-pass values in the global block (`uRenderMask`, `uPickType`, `uMarkingType`, ...)  | a new snapshot per pass                                     | one global slot per pass kind, or pipeline overrides (layer A)                            |
| global textures per pass (`tDepth`, `tDepthAlpha`)                                     | bind group 0 differs per pass                               | a bundle per pass kind (its group 0 is fixed while the targets are)                       |
| object culling (`isVisible`: frustum, lod distance, Hi-Z) per frame                    | the draw set changes with the camera                        | keep it on the CPU, the visible set in the validity key (below); or on the GPU (layer D)  |
| instance culling of instanced items (`r.cull`: multi-draw ranges, lod levels)          | draws and draw blocks change with the camera                | draw those items directly after the bundle; later layer D                                 |
| solid interior caps (stencil clears, cap pass and clip index between draws)            | clears are not bundle commands, globals change mid-item     | caps drawn directly after the bundle (as the plan's layer A for caps)                     |
| uniform writes in `render` (`materialBlock.write()`, `objectBlock.write()`)            | replay skips `render`, so changed values would not reach the GPU | dirty items flushed before replay (observers on the values), no per-draw loop        |
| pipelines compiling in the background (`getPipeline` returns undefined)                | the draw is missing from the recording                      | record only when every draw's pipeline is ready (`pending.isReady`)                       |
| arena growth (uniform slots or storage pages reallocated)                              | recorded bind groups and offsets point at old buffers       | an arena version in the validity key                                                      |

The per-object bind groups themselves do not break replay: their offsets are persistent slots. The upper bound (step 1)
replayed 1008 draws with per-object groups at 0.1 ms submit. Layer B (shared records, one group per program and page)
lowers the commands per draw and the recording cost, it is not needed for validity.

### When frames are drawn with an unchanged scene

Mol\* draws only when something changes (the idle loop sleeps otherwise), so the gain is in frames drawn repeatedly:

- **Camera moving** (rotate, zoom, spin, camera transitions): globals change every frame; with a fixed global slot and
  the visible set unchanged (the whole scene in view, the common case) every frame replays. Views from inside or close
  zooms change the visible set often: those frames record again (no gain, recording cost).
- **Static camera, accumulating frames**: temporal multi-sampling, illumination (tracing) and SSAO accumulation draw
  tens of frames with only the jitter or the sample index changing: replay throughout.
- **Hover and selection**: markers are data writes; the pick passes (asynchronous reads, several per second on hover)
  and the marking pass replay while the camera rests.
- **Animations of uniforms** (wiggle, tumble, `uTime`): replay.
- **No gain**: trajectory playback, state animations, representations being built or updated, lod switching every
  frame: the draw set changes each frame.

### Where the CPU time is today

Per frame on WebGPU (RTX 4070 Ti, `render/*`): `boxes-1000` 5.5 ms submit (1008 draws), `units-1000-spacefill-unmerged`
about 4.6 ms; scenes of one structure (`viewer-7qpd-rotating`, the 60-copy assemblies) need a few draws and little
CPU. Merged render objects (`renderMerged`) already turn many units of one representation into one draw. So bundles
matter for scenes of many distinct objects: many structures or shapes (MVS primitives, measurements, labels per
object), mesoscale scenes of many entities, and many-pass frames (illumination, pick, marking) of such scenes.

### Options

**1. As planned (A, B, C, D in order).** Pass values out of the globals, shared records and storage pages, then bundles,
then GPU culling.

- Pros: the end state has the least CPU per frame in every case (recording cheap, one call per draw, culling on the
  GPU keeps bundles valid while moving inside a scene); B also speeds frames that cannot replay (trajectories).
- Cons: B is the largest change (render items, bindings, theme data as storage pages, geometry pulling) before any
  bundle gain; storage pages need fallbacks for data beyond the binding limits; most risk sits before the first
  measurable benefit.

**2. Lean bundles first (A, then C with CPU culling as the key; B and D later).** Fixed global slots and per-pass
values, a validity key with the visible set, dirty flushing, instanced items and caps drawn directly.

- Pros: the measured gain for the cases above with the least change: no new data layout, per-object bind groups kept;
  culling semantics unchanged (the CPU decides, a different visible set records again); each piece measurable on its
  own; B and D keep their value and fit later (B makes recording cheap, D keeps instanced items in bundles).
- Cons: recording costs about a direct frame, paid whenever the visible set changes (zooming through a scene, views
  from inside); instanced items with instance culling (assemblies, mesoscale) stay direct, so the largest instanced
  scenes gain little until D; one more path in the renderer (record or replay vs. direct) and invalidation bugs show as
  stale frames, so tests must compare replayed with direct frames.

**3. No bundles: only cheaper draws (B).** One call per draw through shared records.

- Pros: helps every frame, also those that cannot replay; no validity tracking.
- Cons: still a loop and calls per object every frame; the gain per draw (fewer `setBindGroup`, vertex buffer calls) is
  a fraction of replay's (5.5 ms to an estimated 2–3 ms against 0.1 ms); the largest change of all options.

**4. GPU-driven first (D).** Culling and lod on the GPU, indirect draws.

- Pros: the largest instanced scenes (mesoscale, assemblies viewed from inside) lose their CPU culling and multi-draw
  list rebuilds, also without bundles.
- Cons: little for scenes of many separate objects (each still a draw); compute per frame; Hi-Z on the GPU side to
  build; WebGL keeps the CPU path.

### Pros and cons of render bundles in general

- Pros: CPU submit of replayed passes near zero (5.5 to 0.1 ms measured), whatever the object count; the main thread
  free for interaction while frames draw; identical images (the same commands); WebGPU core (Chrome, Firefox, Safari),
  no optional feature; several passes per frame (opaque, depth, pick, marking, emissive, the frames of accumulation)
  benefit from the same scheme.
- Cons: valid only while the draw set is unchanged, so no gain while the scene changes and a recording cost when it
  does; validity tracking across values, pipelines, arenas and targets (stale frames when missed); a bundle per pass
  kind and target layout (memory for the recorded commands, small next to the data); WebGL has no equivalent, so the
  renderer keeps its direct path (the neutral `createDrawList` records there by drawing); recording waits for
  pipelines compiling in the background; per-frame CPU work outside the draws stays (culling tests, uniform flushes,
  passes and postprocessing).

### Recommendation

Option 2: A first (measured on its own: it also removes the per-pass global snapshots of direct frames), then lean
bundles for the opaque pass behind a renderer prop, with replay-against-direct checks in the visual suite; extended to
the pick, depth, marking and emissive passes once the opaque pass holds. B and D stay planned for what option 2 leaves:
recording cost on changing views, and instanced scenes.

Measurements before C (done): `pnpm benchmark --backends webgpu --filter 'render-bundles/*'` (on-demand cases;
`WebGPUContext.recordRenderBundle` and `executeRenderBundles`, which C builds on; the benchmark wraps `renderOpaque`).
RTX 4070 Ti, Chrome, per frame:

- **Recording costs no more than drawing directly**: the opaque pass recorded and replayed every frame, submit
  `boxes-1000` 4.8 ms (direct 5.05 ms), `units-1000-unmerged` 4.1 ms (4.86 ms); 30 of 30 recordings valid. So
  re-recording when the key changes costs what a frame costs today: no case gets slower.
- **The visible set while rotating**: the whole scene in view, unchanged in 0 of 29 frames (every frame replays);
  zoomed in to 0.3 and 0.1 of the camera distance, it changed in 29 and 26 of 29 frames (those frames record).
- **The visibility key**: 0.06–0.17 ms per frame for 1000 objects (built as a string; a bit set is cheaper).
- **Drawing what culling skips** (opaque pass on the GPU, `--gpu-timing`): zoomed to 0.3 (half the boxes visible) 0.75
  ms culled, 1.06 ms unculled; to 0.1 (8% visible) 0.25 ms against 1.03 ms: under 1 µs of GPU per off-screen draw. The
  CPU of the unculled pass replayed instead of submitted (about 5 ms) outweighs it while a fair share is visible (0.3: a
  frame of about 1.8 ms replayed unculled against 2.9 ms culled), not when little is (0.1: about 2 ms against 1.5 ms).
- Not measured: flushing changed object data before a replay, design-dependent (observers mark dirty items: nothing
  to do when nothing changed).

So C with the visible set in the key is never slower than today and gains the whole CPU time when the visible set
holds (the scene in view, static-camera frames); replaying zoomed views unculled is a later option, chosen by the
visible fraction.

## Steps

1. **Upper bound** (done): in a test-only path, record the opaque pass of `render/boxes-1000` once and replay it;
   measure the CPU per frame and check the image. Tells how much layers A–C can gain: the submit time drops from 5.5 ms
   to 0.1 ms, see above.
2. **Global offsets that repeat from frame to frame** (layer A, done): a ring of the global blocks' own, see layer A.
3. **Cap clip index in per-item draw slots**: not needed for bundles (caps are drawn directly after them); only if caps
   should be recorded too.
4. **Data arrays as storage-buffer sections** ([storage-buffers.md](storage-buffers.md) steps 1–4), shared arena per
   context (layer B, data).
5. **Object and material records** in shared storage pages, object index via the instance data and `firstInstance`, one
   bind group per program and page (layer B); measure `render/boxes-1000` and the units cases (expected: about one call
   per object), and a case with theme data beyond one page.
6. **Render bundles** for the opaque, pick, depth, marking and emissive passes (layer C) behind a renderer prop, an
   opt-in visual profile rendering every case with bundles, check cases changing state between frames, then the default.
   Every pass but illumination's inputs done (`drawLists`, `--draw-lists`, see layer C); data changes without
   re-recording next.
7. **GPU-driven culling and lod** for instanced structures (layer D), first for spheres with lod levels; a large
   assembly benchmark case with the camera inside.
8. **Optional accelerations** with core fallbacks (checked on this machine, Chrome 156: `chromium-experimental-multi-draw-indirect`
   only with `--enable-unsafe-webgpu`, so not for users yet): multi-draw indirect where available (one call per pass), immediate
   data if it helps beyond `firstInstance`.

## Open questions

- Page size: the device limits, or smaller pages so that unloading frees memory and growth stays cheap.
- Instance data in shared sections: transforms of instanced items (assemblies) are already per instance; non-instanced
  items get one entry each, an indirection per vertex to measure against today's per-item vertex buffers.
- How many override-constant pipeline variants the passes need, and their compile time (pipelines compile in the
  background, `linkPrograms` creates known ones ahead).
- Shared arena fragmentation and growth on scenes that load and unload often (trajectories, many structures).
- Whether frustum culling without layer D costs more GPU time than it saves CPU time on large scenes viewed from inside.
