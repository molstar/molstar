# Plan: Draw submission on WebGPU

## Status

- Paused after step 1: the upper bound is measured (see below) and shows the gain; steps 2–8 are not started and wait
  until draw submission becomes the priority again. The measurement code was removed afterwards (see below to redo it).
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

## Steps

1. **Upper bound** (done): in a test-only path, record the opaque pass of `render/boxes-1000` once and replay it;
   measure the CPU per frame and check the image. Tells how much layers A–C can gain: the submit time drops from 5.5 ms
   to 0.1 ms, see above.
2. **Pass values as pipeline overrides**, the global block in a persistent slot per global bindings (layer A); visual
   suite unchanged; measure `render/*`.
3. **Cap clip index in per-item draw slots** (one per clip index, dynamic offsets per cap draw) instead of global
   changes between draws.
4. **Data arrays as storage-buffer sections** ([storage-buffers.md](storage-buffers.md) steps 1–4), shared arena per
   context (layer B, data).
5. **Object and material records** in shared storage pages, object index via the instance data and `firstInstance`, one
   bind group per program and page (layer B); measure `render/boxes-1000` and the units cases (expected: about one call
   per object), and a case with theme data beyond one page.
6. **Render bundles** for the opaque, pick, depth, marking and emissive passes (layer C) behind a renderer prop, an
   opt-in visual profile rendering every case with bundles, check cases changing state between frames, then the default.
7. **GPU-driven culling and lod** for instanced structures (layer D), first for spheres with lod levels; a large
   assembly benchmark case with the camera inside.
8. **Optional accelerations** with core fallbacks: multi-draw indirect where available (one call per pass), immediate
   data if it helps beyond `firstInstance`.

## Open questions

- Page size: the device limits, or smaller pages so that unloading frees memory and growth stays cheap.
- Instance data in shared sections: transforms of instanced items (assemblies) are already per instance; non-instanced
  items get one entry each, an indirection per vertex to measure against today's per-item vertex buffers.
- How many override-constant pipeline variants the passes need, and their compile time (pipelines compile in the
  background, `linkPrograms` creates known ones ahead).
- Shared arena fragmentation and growth on scenes that load and unload often (trajectories, many structures).
- Whether frustum culling without layer D costs more GPU time than it saves CPU time on large scenes viewed from inside.
