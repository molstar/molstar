# Plan: Draw table, objects as flat data

## Status

- Steps 1 to 5 done (4 with chunked recordings): the table is the renderer's path. The GPU side went elsewhere: GPU
  culling tests the cells of merged objects with their own buffers ([merged-render-objects.md](merged-render-objects.md)
  step 6, `Renderer.prepareCull` over the renderables), not table rows; per-row indirect draws were measured a loss
  (step 7: an indirect draw replays at about 0.3 µs whether culled or not). So step 6's row buffer had no consumer
  (its test only) and was removed; the design below records it. Layer B of
  [webgpu-draw-submission.md](webgpu-draw-submission.md) was dropped for submission (shared object groups instead).
  Replaced the per-object planning steps of that plan (Planning cost, steps 1–5): those trimmed the cost of an object
  visit (1.18 to 0.85 ms planning for `units-1356-gaussian-highlighted-drawlists-rotating`), but every pass still
  visited every object through calls into its values.
- Left: the setup cost of 100k objects (a plan of its own, not written).

## Why

- Scenes from user data reach 10k–100k small objects. Per frame and pass the renderer walks all of them: state checks,
  bounding sphere scaled and tested, lod, occlusion, pass checks (`checkOpaque`, ... reading `r.values.*.ref.value`,
  four dependent loads into scattered objects each), caps, culling, render state, value versions. About 200 ns per
  object and pass with draw lists, more on the direct path; it grows linearly and is redone while nothing changed.
- What decides the draws are few facts per object that change rarely (on commits, value and state changes) and the view
  that changes per frame. Kept as flat columns, maintained on change, the per-frame work becomes tight loops over typed
  arrays: about 10 ns per row for the view test, bit tests per pass.
- The same rows, uploaded, were to be what GPU-driven culling reads (layer D); GPU culling came per merged object
  instead (Status).

## Design

Render objects, renderables and visuals change where the table is better served that way: the facts and state changes
come from their owners instead of being derived from the outside.

### Scene: the facts (renderer independent)

- **Slots**: each renderable gets a stable slot (a free list, reused after removal), not its index in `renderables`; the
  draw order is a separate `order` array of slots, rebuilt from the group's order when objects were added or removed.
  One table per group (`Scene.Group.table`, `gl/draw-table.ts`): the renderer draws `primitives` and `volumes`
  separately, and the scene itself as a group (`Scene.table`, all renderables, `renderBlended` of helper scenes; told of
  the scene's sort, `invalidateOrder`); groups only planned with (precompile) have none and take the object path. Stable
  slots let rows be uploaded as they change and become the object record index of layer B.
- **Rows** in one `ArrayBuffer` read through a `Float32Array` and a `Uint32Array` view (no `DataView`), laid out as a
  WGSL struct (std430, 32 bytes): bounding sphere (center, radius; unscaled), lod distances (min, max), fact bits, state
  bits. Interleaved because the view test reads nearly all of a row and rows upload as they are; CPU-only columns
  (versions, the renderable) are separate arrays by slot.
- **Fact bits**: what the renderer's checks read, as bits: draw count > 0, alpha 1, alpha 0, transparency any and full,
  marker any and full, emissive, direct volume, text or image, fuzzy points, x-ray, transparent back faces opaque and
  off, double sided, single sided (`uDoubleSided` given and false), flip sided, reflection, solid interior, instance
  grid, several instances, lod levels (`DrawFact`, `getDrawFacts`). Clip
  caps stay with the renderable (`getClipInfo`), read for solid-interior rows only.
- **State bits**: disposed, visible, pickable, color only, opaque, write depth (`DrawState`).
- **Maintenance**:
  - Facts: derived from the values (`getDrawFacts`, one place for what a fact means), refreshed when a cell it reads
    changes (`DrawFactCells`): the table observes those cells per slot with an observer whose `changed` setter queues
    the slot (the `ValueCell.Observer` interface as is). Queued rows are written by `refresh`; a version per row and one
    for the table advance.
  - State: observable (`RenderableState.create`, accessors notifying observers on change), so that every writer keeps
    writing fields: about 50 places write state directly (visuals, geometry props, helpers, extensions), moving them to
    setters would be churn with stale rows for any missed or new writer. The geometry utils create observable states
    (`createRenderableState`); plain states (built by hand) stay as they are and are compared on `refresh`.
  - Values: any change of an object's values (all its cells, not only the facts') advances its values version and the
    table's values epoch at once (`valuesVersions`, `valuesEpoch`); volatile cells are compared on `refresh`. The
    renderer reads them for its plans and draw lists (data changes, flushed), it keeps no watches of its own.
  - Adds and removes: slots taken and freed with the scene's adds, removes and clear; observers removed with them.
- The scene's aggregates (marker, emissive, wiggle, opacity averages, transparency minimum, `hasOpaque`,
  `hasAnimation`) come from the columns instead of walks over the values: what each object adds is written with its
  row (`DrawTable.contributions`, the `Animated` fact), the scene sums over the visible rows of its primitives.
- Changed objects only: the table lists the slots whose values changed (`takeChanged`); `Scene.updateChanged` updates
  those objects (marking, which changes a few of many per hover) instead of all. The scene's visibility follows the
  table's `visibleVersion` (rows whose visibility changed, added or removed) instead of a hash over all objects.

### Renderer: what the facts mean for its passes

- **Per view**: one loop over the rows writes a visible byte per slot (`getRowVisibility`): state, draw count, frustum,
  lod, then the occlusion test for those passing; inlined plane tests on the row's floats, nothing per object called.
  Keyed by the table's version, the view's planes by value (the frustum's and the camera plane, `viewPlanes`: passes
  of several cameras alternate within a frame, e.g. the camera helper's, a counter of view changes would never hold),
  the model scale and the occlusion test with the frame (Hi-Z changes between frames); shared by all passes and caps
  of the view.
  Pickable is the pass's (`isVariantDrawn`). The near-plane cap test stays per drawn row (`getSolidInteriorCaps`).
- **Pass checks**: predicates over the fact and state bits (`FactTest`: `checkOpaque`, `checkTransparent`, ...) are the
  definition; objects outside a table (pipeline requests with stand-ins, groups without a table) get their facts from
  their values (`getCheckFacts`, the stand-ins hold those values) and take the same predicates, so there is nothing to
  keep equal. A drawn row's sphere and culling facts come from the row (`setRowBoundingSphere`, `cullFacts`), not
  from the values: with every object in view, reading `r.values.boundingSphere` per object was half of planning.
- **Culling**: `r.cull` only for rows that cull instances. The renderable's `drawVersion` advances when the ranges
  it draws change: `cull` and `cullSimple` compare their key after writing it (`createDrawVersion`), `uncull` only when
  culling was on; read as a number (compared on every read it cost five times the old key).
- **Render states** from the facts (`getObjectState(facts, ...)`: geometry type, flip sided, single sided,
  reflection), cached per pass state, cull and front face (`getFaceState`): no values read per draw; a change that
  alters a render state changes the row.
- **Replay decision**: a recording keeps a snapshot (`Snapshot`, written in place): per draw its slot, object, render
  state and program, values and draw versions in typed arrays, what it binds (`digests`), the caps as a flat list
  (slot, versions, caps), and what the walk read besides the rows (table version and values epoch, view planes, model
  scale, near plane, whether caps are drawn). Those unchanged (`isSnapshotCurrent`, a few scalars and 28 planes):
  replay without walking. Otherwise one walk (`walkPass`) compares the pass's draws with the snapshot, stopping at the
  first difference: equal replays, only values versions differ flushes those draws (`flushDataChanges`), else records
  again and the walk writes the snapshot. Programs change without the rows (`setTransparency`): the scene tells its
  tables (`DrawTable.invalidate`). Plans, plan keys, `PlanInputs` and `PassMembers` are gone.
- **Chunks**: a pass's recordings are chunks of 64 objects of the table's order (`PassDrawList`, `PassPart`), then
  one of the caps drawn after all objects; each chunk replays, flushes or records on its own (a snapshot each), so that
  a change of what a pass draws records the chunks it touches: a hover sweeping across objects (the marking passes draw
  the marked objects, the depth pass skips fully marked ones), objects entering and leaving a zoomed view while
  rotating, visibility toggles. Chunks drawing nothing keep no recording; unchanged inputs since the pass was drawn
  replay all chunks without a check of their own; otherwise a chunk of objects whose own rows are unchanged (slots in
  order, row and values versions per position, the table's `invalidations`) replays without walking, while other
  objects change; after a change of the view only (a camera rotating), a chunk whose rows are as visible as when drawn
  replays without walking, unless a visible row's draws follow the view otherwise (instances culled, levels of detail,
  solid interior caps: `isViewDependent`); WebGPU executes the bundles of consecutive replays in one
  `executeBundles` (`executeRenderBundles` queues them until the next command of the pass). A few whole-pass recordings
  kept per list instead did not help: sweeps across thousands of objects never return to a recorded state.
- The direct path (`renderPass`, WebGL) draws from the same per-pass slot lists.

### GPU: the rows uploaded

Built (step 6), then removed without a consumer (Status).

- Rows to a storage buffer, changed ranges written per frame (`DrawTableBuffer`, `draw-table-buffer.ts`): the rows
  by slot as the table lays them out, the order of slots; rows written since the last update merged into a few ranges
  (`TextureChanges.mergeRanges`: gaps under 256 rows merged, all past 8 ranges or half the rows), the order where it
  changed, buffers outgrown replaced. Through `ComputeContext` (storage buffers): WebGL has none and keeps the CPU
  loops. The slot is the object record index.
- Layer D reads them: a compute pass tests rows (and instances of rows that cull) against the frustum, lod and the
  previous Hi-Z, writes per-pass indirect arguments and compacted lists; bundles record indirect draws and stay valid
  while the camera moves. The CPU loops stay the fallback (WebGL, no compute).
- Many small objects need fewer draws, not only cheaper planning: one indirect draw per object is still 100k draws.
  Batching per program over shared records (layer B, multi-draw indirect where available) follows.

## Steps

1. **Benchmarks** (done): `draw-table/boxes-10k*` (on demand), 10k boxes each its own representation, still and
   rotating, direct and with draw lists, one highlighted, with the setup times. WebGPU, NVIDIA, per frame: direct 75.7
   ms submit (155 highlighted), draw lists 0.56 ms unchanged, 14.4 ms rotating (34.2 highlighted); `isVisible` is three
   quarters of planning, its own loads into the objects, not calls. Replaying 10k draws waits about 50 ms (38 rotating),
   of which the GPU executes 5.2 ms (`--gpu-timing`): the rest is Dawn in the GPU process, replaying the bundle's
   commands into the D3D12 command list, about 4.5 µs per draw (hidden behind the 75 ms submit on the direct path). The
   draw count needs batching (layer B) as much as planning needs the table. Setup: representations 0.8 s, add and commit
   2.1 s (render items with their textures and bind groups), first frame 2.2 s. 30k (setup 8.6 s, first frame 13.4 s,
   both superlinear, replay 146 ms) and 100k take minutes to set up: not cases.
2. **Scene table**: slots, `order`, rows and CPU columns, facts and observable state, unit tests (`draw-table.test.ts`)
   done; no measurable cost on the benchmarks (`units-1356` direct and with draw lists, `boxes-10k` setup and frames).
   Aggregates from the columns done (see the sweep in step 4). That the facts equal the renderer's checks is
   tested with the predicates (step 3).
3. **Per-view visibility and pass predicates** (done), used by `planPass` and `renderPass`. Measured back to back
   against step 2 (submit per frame, the machine slower than in step 1): `boxes-10k-drawlists-rotating` 19–21 to 11.2
   ms, `boxes-10k-rotating` (direct) 100 to 88 ms, `units-1356-gaussian-highlighted-rotating` (direct) 25.5 to 21 ms,
   with draw lists about equal (1356 objects). What remains per drawn row: render states (`getObjectState`, step 5),
   the plan inputs per replay (`setPlanInputs`, step 4), `uncull`, the values watches and versions.
4. **Snapshots** in place of plans (done); the draw version of culling renderables. Measured back to back against
   step 3 (submit per frame): `boxes-10k-drawlists` (still) 0.6 to 0.15 ms, `-rotating` 7.0 to 6.1 ms,
   `-highlighted-drawlists-rotating` 14.2 to 11.3 ms, `units-1356-gaussian-highlighted-drawlists-rotating` 1.1 to 0.88
   ms; `boxes-1000-zoom0.3-drawlists-rotating` (records every frame) 4.2 to 4.0 ms.
   Then **chunked recordings** (done, see Chunks), measured back to back, per frame:
   `draw-table/boxes-10k-drawlists-sweep` (a whole box highlighted per frame, another each time) 187 to 61 ms (about 2
   chunks recorded per pass instead of the whole pass), `render-bundles/boxes-1000-zoom0.3-drawlists-rotating` 4.0 to
   1.2 ms; full replays of 10k objects 0.1–0.8 ms more CPU (the chunks' walks), as much less GPU wait. What remains of
   the sweep's 60 ms is not the draw lists: `Canvas3D.resolveMarking` updates every render object of the scene per hover
   (`scene.update(undefined, true)`, 60%), the scene's averages are computed again over all objects (about 10%), the
   shader requirements check the scene. Updates of the changed objects only: 61 to 28.6 ms (direct 257 to 159 ms); the
   aggregates from the table then: 18.5 ms (direct 112 ms); chunks checked by their own rows (not walking when others
   changed) then 17.9 ms, the submit's CPU 5.2 to 2.9 ms: the rest is GPU wait (10k draws, step 6 and layer B);
   shared object groups (layer B) then 8.9 ms.
5. **Render states** from facts (done): `boxes-10k-drawlists-rotating` 4.6 to 2.3 ms,
   `-highlighted-drawlists-rotating` 11 to 4.3 ms, `units-1356-gaussian-highlighted-drawlists-rotating` 1.6 to 1.4 ms;
   the direct path unchanged (encoding dominates).
6. **Rows on the GPU** (done, then removed without a consumer, see Status): `createDrawTableBuffer(compute, table)`,
   unit tested (rows and order as the table's,
   only changed rows written, outgrown buffers replaced); created and updated by its consumer (step 7), not per frame
   without one. The object records of layer B are not needed for submission (shared object groups, layer B step 4).
   Its upper bound measured
   (`draw-overhead/*`, raw WebGPU, 10k tiny cubes in one render bundle, NVIDIA, D3D12): binding as Mol* does (per draw
   an object and a draw block bind group, five vertex buffers and an index buffer of its own) 37 ms per frame, only the
   buffers of their own 24 ms, only the bind groups 10 ms, shared bindings with a dynamic offset per draw 1.1 ms,
   objects by `firstInstance` 1.0 ms, one instanced draw 0.2 ms. The replay costs per distinct resource a draw binds
   (about 0.4–0.5 µs each), not per command: layer B removes nearly all of it, geometry in shared storage (vertex
   pulling, a shared index buffer) first, then the object records. Planned with layer B in
   [shared-draw-resources.md](shared-draw-resources.md): shared buffers bound per draw at an offset cost 2.2 ms, so
   geometry needs sub-allocation, not vertex pulling.
7. **GPU culling** from the rows (layer D), indirect draws in bundles. Measured first (`draw-overhead/*`, 10k draws in
   a bundle, shared bind group with dynamic offsets, NVIDIA, D3D12): direct 1.1 ms per frame, indexed indirect 4.1 ms
   (about 0.3 µs more per draw), half of them culled (instance count 0) 4.1 ms as well: an indirect draw costs its
   replay whether it draws or not. Against today's CPU culling (`boxes-10k-drawlists-rotating` 5.5 ms, submit 2.5 ms;
   `boxes-1000-zoom0.3-drawlists-rotating` 1.3 ms; `assembly-60-*-rotating` 0.1 ms submit), per-object indirect
   draws would add GPU time for every object, culled or not, to save CPU time of the walks: a loss without fewer
   draws (compaction with multi-draw indirect, not in WebGPU; or batching per program). Batching was dropped
   ([batched-draws.md](batched-draws.md)); GPU culling is planned for merged render objects, where one object holds
   many units and their instances ([merged-render-objects.md](merged-render-objects.md)).

Each step is checked by the unit tests, `--draw-lists` and `webgpu-draw-lists`, and the direct visual suite where the
direct path changes.

## Later: membership as GPU data

Where the features are available, chunks could give way to recordings that never change with what a pass draws: every
candidate object of a pass recorded as an indirect draw, whether it takes part (marked, visible, in the pass) written
as its instance count (0 or n) into a per-object arguments buffer, a marker or visibility change then a few bytes of
data. The same arguments GPU culling (step 7) writes. Needs indexed indirect draws in the render items
(`drawIndexedIndirect`, only `drawIndirect` of GPU surfaces today), `indirect-first-instance` (optional) or another way
for culled ranges starting at an instance offset, and the per-draw cost of indirect draws in bundles measured
(`draw-overhead/*`) against direct draws.

## Open questions

- Volatile cells (versions without notification): polled per row that has them.
- Whether helper scenes (debug, handles, pointer) use the table or stay on the per-object path.
- Creating render objects and geometry for 100k objects costs more than drawing them; a plan of its own.
