# Plan: Geometry residency, CPU copies released after upload and playback without allocations

## Status

- **Current**: proposed, nothing done. One notion of where geometry lives, on CPU and GPU:
  - the CPU copy of geometry that can be computed again is released once it is uploaded, and computed again (or read
    back) when something needs it: a lost context, a backend switch, export, smoothing, a props change that rebuilds
    from it;
  - geometry rebuilt every frame (trajectory playback) allocates nothing once warm: each geometry keeps its buffers
    with slack and hands out views, and each computation's intermediates come from a scratch arena reset after it.

  Step 1 (accounting and measurements) decides how far the rest goes.

## Today

- **Every CPU copy is kept.** The geometry object and the render object share the same `ValueCell`s (`Mesh`
  `vertexBuffer` is `aPosition`, `Spheres.createValues`'s `tPositionGroup` next to `centerBuffer` and `groupBuffer`).
  Visuals keep the geometry, the render object and the position iterator, whose closures capture the arrays
  (`spheres.ts`, `mesh.ts`, `units-visual.ts`). Renderables spread the values. Merged objects' `SegmentedArray` parts
  reference the members' arrays. Nothing ever drops a copy.
  - Spheres hold their data twice on the CPU: centers (12 B) and groups (4 B) as `ValueSpec`s that are never uploaded,
    plus `tPositionGroup` (16 B) built from them, against 16 B on the GPU. Points work the same way.
  - `Mesh.meta.originalData` keeps the index buffer from before `uniformTriangleGroup`.
  - Textures keep their load source (`getSource`) for `transferTexture`, so volume grids, images and font atlases
    persist by design.
- **Restore is a rebuild that reads the CPU copies** (`backend-switch-followups.md`, step 4). `Scene.setContext` adds
  the same render objects to the new context and creates renderables from `values`. Objects with foreign resources
  (GPU surfaces: textures or arena chunks of the lost context) are dropped. Their visuals compute the geometry again:
  `Visual.foreignResourcesCheck` in `prepareUpdate` forces `createGeometry`. That is the one path today for geometry
  that has no CPU copy.
- **Readback on demand exists for texture meshes only**: `GeometryUtils.prepare` and `GeometryRequest`
  (`geometry-preparation.md`) for position iterators, and geo-export's own `readTextureData`.
- **CPU readers after upload**:
  - bounding spheres (cached by a version hash, `mesh.ts`, `center-group.ts`) and shape group bounds;
  - the spheres' LOD and size-factor changes (`setPositionGroup` from `centerBuffer` and `groupBuffer`);
  - position iterators for vertex and volume themes;
  - CPU color smoothing (overpaint, transparency, emissive, substance) and `processValues`;
  - geo-export (`values` of mesh, lines, points, spheres and cylinders);
  - `getGroupCount`, clip volumes, the meshes and debug-helpers extensions;
  - merged objects' layout changes (every part written again from the members).

  Picking, culling and Hi-Z read no geometry.
- **Trajectory playback**: a conformation change sets `createGeometry`. Builders take the previous geometry as their
  first chunk (`MeshBuilder.createState(…, mesh)`, `SpheresBuilder.create(…, spheres)`, about 70 call sites).
  `ChunkedArray.compactTrimmed` returns that chunk in place only when it is a single chunk at least 80% full; otherwise
  it allocates and copies. So a frame of the same count reuses its arrays, and a frame that grows or shrinks by more
  than 20% allocates (surfaces, selections that change with the frame). A shrinking frame leaves a trimmed copy as the
  next frame's first chunk, so sizes going up and down allocate (a chunk, then the compaction copy) every frame.
  - The builders overwrite the previous frame's arrays while it is still drawn and referenced (its geometry, its
    position iterator): in-place reuse is already the contract, safe because the GPU holds the drawn copy.
  - Upload is by version: the whole array is written again (`updateSubData` while it fits).
  - Not counted anywhere: the intermediates of a computation (splines, lookups, density grids, element lists),
    allocated per call.
  - GPU surfaces keep their arena and scratch buffers while playing (`gpu-surface-buffers.md`).
- **Accounting is GPU only**: `stats.resourceCounts`, `ResourceByteCounts` and the renderables' `getByteCount`. On
  WebGPU nothing adds to `attribute`, so the shared geometry pages go unreported. No CPU bytes are counted anywhere.
- **Related**: the geometry store of [shared-draw-resources.md](shared-draw-resources.md) (CPU pages mirroring GPU
  pages, owner pages reused by their own updates) was deferred because it measured little gain for draw overhead. Its
  CPU half is what item 6 below takes up, for allocations rather than draws. `GraphicsContext.release` and `restore`
  (Prism's canvas presence, at most 12 live contexts) are the GPU side of residency.

## Goal

- A scene's CPU memory for geometry shrinks to what is in use: what is drawn and idle lives on the GPU only.
- Every reader gets its data, through one asynchronous request.
- Playback allocates nothing per frame once it has reached its high-water mark.
- CPU and GPU bytes are reported side by side.
- Geometry and render objects keep working without a graphics context, as today (item 10).

## Design

1. **Residency per geometry**: `resident` (CPU and GPU, as today) or `released` (GPU only, plus a way to get the data
   back). It is a visual's state, because the visual holds the recipe: its `createGeometry` with the same inputs
   (structure, unit, props, theme sizes), the same path `foreignResourcesCheck` takes today.
   - Why recompute, not keep a compressed copy or spill to IndexedDB: the inputs (the model) are kept anyway, and the
     recompute path exists and is tested by every context loss.
   - **Theme layers too**: color, size, overpaint, transparency, emissive, substance and clipping are released and
     acquired like the geometry. Their recipe is the visual's theme update (`createColors`, `createSizes`, the layers'
     producers from the representation's state). Why: per group they are small next to a surface's vertices, but not
     next to a sphere's 16 B (color alone is 3 B a group), so for sphere and point scenes they are a large part of
     what is held.
   - **Except the marker layer**: marking actions change it in place (`applyMarkerAction`), so it is the canvas's
     highlight and selection state, not derived data. It stays resident, at 1 B a group.
2. **Releasing without a change**: the arrays in the geometry's cells are swapped for a shared empty array of the same
   type without a version bump (`ValueCell.releaseRef`, or the geometry kind's own `release`). So render items upload
   nothing, version-hashed caches (bounding spheres) keep their values, and draw counts stay.
   - The geometry is marked released (`Geometry.isResident`). The position iterator and anything else holding the arrays
     is dropped with it and created again on acquire.
   - A dev check: reading a released geometry outside an acquire throws. A test-mode proxy, not in production.
3. **Getting it back: `acquire`**: `Visual.acquireGeometry(): Promise<Geometry>` generalises
   `GeometryUtils.prepare` and `GeometryRequest`, for the geometry and the released layers. Two strategies; which one
   on-demand readers use while the context is alive is measured per kind and backend (step 3), not assumed:
   - **Read back** from the GPU (WebGPU `mapAsync`, WebGL2 `getBufferSubData` with a fence): exact, no CPU work, but
     a staging buffer and a wait for the GPU. Only while the context is alive.
   - **Recompute** by the recipe: the only way after a loss or a switch; also for props changes that rebuild from the
     data (the spheres' LOD, which then just recomputes).
   - Recomputed geometry must match what the GPU holds: same vertex and group counts, checked, and re-uploaded on a
     mismatch (non-deterministic producers are not eligible).
   - Acquired geometry is resident again until the policy releases it.
4. **Readers through acquire**: geo-export, CPU color smoothing and `processValues`, position iterators for vertex
   themes, shape bounds, `getGroupCount`, clip volumes, the extensions' mesh utilities. Each is already asynchronous
   or runs in a task; the synchronous ones (`getGroupCount`, shape bounds) keep what they need as metadata at
   release (counts, group bounds), so they never acquire.
5. **Context loss and switches**: a released object is treated like one with foreign resources. `Scene.setContext`
   drops it and the visual recomputes in `prepareUpdate`. Nothing new on the restore path; it widens the existing one.
   A context released on purpose (`GraphicsContext.release`) with released geometry restores by recomputing,
   proportional to the scene (measure: Prism's canvas presence restores a background canvas this way).
6. **Rebuilt every frame without allocating**: two lifetimes, two granularities.
   - **Outputs live until replaced: per geometry.** A geometry's arrays are its own buffers with slack. Builders write
     into them as today, and compaction returns a view (`subarray(0, count)`) instead of a trimmed copy. The buffer
     grows to the high-water mark × 1.25 when a frame does not fit, and never shrinks while frames come (a shrinking
     frame is a shorter view). So sizes going up and down allocate only at a new high-water mark.
     - Why per geometry, not one arena shared by the frame's visuals: lifetimes are mixed. Visuals a frame does not
       touch (an unchanged ligand, units off-screen or of another model) keep their geometry, so a bump arena reset
       per frame would free live data. The geometry's lifetime is already the visual's, and replacement in place is
       already the contract (Today). Nothing is double-buffered, because nothing reads the old frame's CPU arrays
       while the new one is built that does not read them today.
     - Why not one shared paged store (the deferred geometry store of `shared-draw-resources.md`): it buys contiguous
       uploads for draws, which measured little, at the cost of an allocator; per-geometry slack gets the allocation
       win alone.
     - Readers take the geometry's counts (`vertexCount`, `sphereCount`, …), never an array's length; one that takes
       the length is a bug, slack or not. The view's length is the count anyway.
     - Uploads by range: the view, written as today (`updateSubData` while it fits). The GPU buffer keeps its size with
       the CPU buffer's.
   - **Intermediates live within one computation: a scratch arena per thread.** A bump allocator of typed arrays, taken
     by `createGeometry` and its helpers for their temporaries and reset when the call returns. Always on, not only
     during playback; one per worker as well. Only for the producers step 1 finds allocating per frame.
   - **Knowing a geometry is rebuilt every frame** matters only for when to give the slack back, not for correctness:
     the slack costs at most 25% and reuse needs no detection. The signal is time alone: two rebuilds of a geometry
     within about a second mark it as streaming, a second without one clears it.
     - Why time, not the reason for the rebuild or a hint from the app: any rebuild that comes again soon gains from
       reuse (trajectory frames, a slider dragged through a parameter, coordinates set by an embedder), whatever drives
       it. The plugin's `isAnimating` says neither which geometry rebuilds (it is set for camera animations and
       snapshot transitions that rebuild nothing) nor anything the timing does not.
     - While streaming the slack does not shrink, even for frames well below the high-water mark: the next frame may
       grow again.

     While streaming, the geometry is neither released (item 7) nor trimmed. When it clears, it is trimmed to its
     count, then idle like any other.
7. **Policy**: a Canvas3D prop `geometryResidency: 'keep' | 'release'`, plus a CPU budget in MiB (under `release`).
   `release` becomes the default once measured (step 5); until then `keep`.
   - Released: geometry idle (not rebuilt or acquired) for some seconds, largest first, while over the budget (0
     releases everything idle).
   - Never released: streaming geometry (item 6), the marker layer (item 1), objects that are not eligible (no
     recipe: shapes from data the visual does not own; non-deterministic producers), and v1 merged members (their
     parts are rewritten from the members on layout changes).
   - The plugin and Prism set it (Prism's many canvases first).
8. **Accounting**: CPU bytes per render object (geometry and theme layers, from the values' arrays) and per context,
   next to the GPU counts in `stats`, the debug console and the benchmark. The WebGPU geometry pages are counted under
   `attribute`. Released and resident are counted separately, and so are the slack and the scratch arenas.
9. **Quick wins found on the way, done regardless of the policy**:
   - spheres and points without their duplicate: `centerBuffer` and `groupBuffer` as views of, or derived on demand
     from, `tPositionGroup` (or the other way round), half of their CPU bytes;
   - `Mesh.meta.originalData` dropped where nothing switches back.
10. **Without a graphics context, as today**: representations, geometry and render objects are built and read with no
    context and no Canvas3D (export to glTF and the other geo-export formats in Node or a worker, tests, servers).
    Nothing here may change that:
    - release needs a context: there is no GPU copy to keep, so without one everything stays resident and `acquire`
      returns the geometry at once;
    - nothing in a visual, a geometry or a render object starts to require a context or a residency manager (the policy
      and the accounting live in Canvas3D and the context, and are absent without them);
    - item 6 (views with slack, the scratch arena) is CPU only and works the same without a context.

## Steps

1. **Accounting and measurements** (the gate): CPU bytes per render object and per context in `stats` (item 8), the
   geometry pages counted. Measure:
   - CPU geometry and layer bytes against the JS heap for typical scenes (1cbs, 7qpd surfaces and cartoon, an assembly,
     the PetWorld model);
   - allocations and GC time per frame of trajectory playback (cartoon, spacefill, a molecular surface; a trajectory of
     a few thousand atoms and one of 100k), split into output arrays and intermediates by producer (allocation
     sampling in the DevTools heap profiler).

   Continue with release if geometry is a large share of the heap; with item 6 if playback allocates per frame
   measurably, its outputs part if the output arrays dominate, the scratch arena for the producers whose intermediates
   do.
2. The quick wins (item 9). Tests: spheres and points values unchanged for their readers; heap per sphere measured.
3. Release and acquire for meshes and spheres and their theme layers (items 1-4): the released state, the
   version-neutral swap, the dev check, read back and recompute, the readers through acquire, metadata for the
   synchronous ones. Tests: every reader on a released geometry or layer gets the same data as before, recompute equals
   the GPU counts, nothing is uploaded on release, marking works on released objects. Measured: read back against
   recompute per kind (mesh, spheres, layers), size (7qpd's surfaces, PetWorld's spheres) and backend (WebGPU, WebGL2),
   for the time to the data and the peak memory (the staging buffer, the intermediates); the faster one becomes the
   on-demand default per kind and backend.
4. Context loss and switches with released geometry (item 5). Tests: `renderer/context/release-restore` and
   `webgpu-recreate` with everything released, `pnpm test:viewer-smoke --switch` and `--fallback`. Item 10: geo-export's
   `prism-viewer.test.ts` (1cbs from the headless plugin's render objects, no WebGL, written as GLB, STL, OBJ and USDZ)
   keeps passing through every step; a case with `geometryResidency: 'release'` set and no context added to it.
5. The policy and prop (item 7); Prism's canvas presence uses it. Measure memory of a Prism session with many canvases
   and the step 1 scenes under `release`, and the cost of acquires in normal use (export, theme changes, restore); with
   the gain shown and no regressions, `release` becomes the default.
6. Item 6: views with slack for meshes, spheres and cylinders (`ChunkedArray` compaction, the builders), the streaming
   mark by time, the trim when it clears; the scratch arena for the producers step 1 names. Measured with the playback
   cases of step 1: allocations per frame at zero after the high-water mark, frame time unchanged or better. Tests:
   sizes going up and down reuse the buffer, readers see the count (not the capacity), a cleared mark trims.
7. Later: cylinders, lines and text released; merged members; textures without their source (volumes and images
   recomputed or read back instead of `transferTexture`).
