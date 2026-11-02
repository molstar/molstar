# Plan: Geometry residency, CPU copies released after upload and playback without allocations

## Status

- **Current**: proposed, deferred for now (nothing done). One notion of where geometry lives, on CPU and GPU:
  - the CPU copy of geometry is released as soon as nothing immediately needs it: once the purpose it was created for
    is served (drawn: uploaded, with the readers that run at creation; exported: written). It is computed again (or
    read back) when something needs it later: a lost context, a backend switch, export, smoothing, a props change that
    rebuilds from it. Release is configurable (keep everything, as today, or release);
  - geometry rebuilt every frame (trajectory playback) allocates nothing once warm: each geometry keeps its buffers
    with slack and hands out views, and each computation's intermediates come from a scratch arena reset after it.

  Step 1 (accounting and measurements) decides how far the rest goes. Its first part measured (below): CPU geometry is
  most of what CPU meshes retain (surfaces as CPU meshes, cartoon: 89 % of 3j3q's at `medium`), little of a sphere
  scene's or of texture-mesh surfaces'; playback allocations not yet.

- **Measured** (`pnpm benchmark --byte-notes`: the render objects' CPU bytes, GPU bytes and the JS heap the scene
  retains; a heap snapshot of the Viewer for what holds the array buffers):
  - 7qpd as cartoon, ball-and-stick, molecular and Gaussian surface (`render/viewer-7qpd-rotating`): CPU geometry
    42.6 MB, layers 2.0 MB;
  - 4V6X as Gaussian surfaces, CPU meshes (`merging/4v6x-gaussian-*`): CPU geometry 80.1 MB, layers 1.1 MB; the heap
    snapshot after loading it in the Viewer (WebGPU): 131 MB of array buffers, of them 45 MB of `Float32Array` values
    (positions, normals, groups), 39 MB of `Uint32Array` (elements) and 6.4 MB of `Int16Array` (the encoded normals the
    render item derives and uploads), the rest the file (12 MB) and the model; JS objects about 23 MB. Geometry is about
    two thirds of what the scene holds;
  - spheres (`render/7qpd-spacefill`, `render/assembly-60-spacefill-rotating`, `render/globule-300k-spacefill`): the
    positions (`tPositionGroup`, 16 B a sphere) and layers 0.3, 0.9 and 7.2 MB, the scenes' heaps 2.6, 4.5 and 29.3 MB:
    small against the visuals and model;
  - the collected heap the benchmark notes for a scene varies between runs for surfaces (98 to 266 MB for 4V6X): the
    snapshot is the reliable measure there.
  - PetWorld (`postsynstubby_2`, particles as spacefill with the mesoscale settings, `merging/postsyn-gpu-hiz`,
    WebGPU, `--heap-snapshot`): 823 MB of array buffers in the page; the render objects' CPU bytes geometry 87.7 MB,
    instances 33.5 MB (and 41.4 MB of merged copies), layers 34.0 MB, about 200 MB. The rest, by retainer paths
    (3.69 M atoms, about 517 k particles):
    - the parsed BinaryCIF, held by the models' source data: its encoded bytes 47 MB, kept for lazy decoding; and the
      columns it decoded, kept by their fields (`binary/field.ts`), 191 MB: `Cartn_x/y/z` as `Float64Array`, 88.5 MB,
      decoded only for the particle bounds (`getAtomSiteBounds` reads `value(i)`; the conformation reads the Float32
      decode, 3 × 14.7 MB, that the model shares); `pdbx_struct_oper_list`'s `matrix`/`vector` columns (about 41 MB)
      are raw `Float64` bytes, decoded as views of the file's buffer, no copy;
    - the model, 260 MB: columns (`auth_seq_id`, `label_seq_id`, `atomId` 3 × 14.7 MB as `Int32Array`, the atom names'
      dictionary indices), the conformation (75 MB, its coordinates 44 MB), the ingredients' structures and
      hierarchies (81 MB);
    - the particles (the mesoscale property), 62 MB: their transforms 33 MB, rotations, coordinates, keys, labels;
    - transforms per instance, five times: the particles' (33 MB), the render objects' `transform` (33 MB), their
      `extraTransform` (33 MB, all zeros: allocated per instance but used only by `setTransform`'s instance
      transforms), `aTransform` (27 MB, a copy of `transform` without a matrix), the merged objects' `aTransform`
      (33 MB);
    - the culler's lists 15.6 MB, picking and upload rings a few MB.

    So for sphere scenes of this size the geometry is not what to release. Released instead: the bounds read the
    Float32 coordinates (88.5 MB) and `extraTransform` is allocated on first use (33 MB): 122 MB less, 676 MB of
    arrays attributed;
  - 3j3q as Gaussian surface with the defaults (`merging/3j3q-gaussian-defaults`, 2.44 M atoms, 1356 chains; texture
    meshes, the small units CPU meshes merged; WebGPU and WebGL2 alike; 36ZA is too large for the benchmark page):
    264 MB of array buffers, the model 120 MB, the file 114 MB (encoded 53 MB; decoded 61 MB, its raw `Float64`
    coordinates the conformation's views), the render objects 27 MB (CPU geometry 9 MB, the encoded normals 1 MB,
    the merged colors' parts 8.7 MB); 88 MB of JS objects.
  - 3j3q as cartoon at medium quality (`merging/3j3q-cartoon-medium`, CPU meshes, merged, WebGPU): 2163 MB of array
    buffers, the geometry 1929 MB (89 %): 38.2 M vertices, positions 459 MB, normals 459 MB and their encoding 135 MB,
    groups 153 MB, elements 724 MB; the merged object's parts are the members' arrays (no copy); the GPU holds 1329 MB
    (attributes 644, elements 686). The model and the file the rest.

  So release (items 1-5) pays for CPU meshes (their positions, normals, groups, elements and the encoded normals: two
  thirds of 4V6X's surfaces, 89 % of 3j3q's cartoon, 1.9 GB), not for the texture meshes surfaces get by default (3 %
  of 3j3q's) nor for sphere scenes; the normals are kept twice, as `Float32` and encoded for upload (3j3q's cartoon:
  459 and 135 MB).

  The default quality (`auto`, `getStructureQuality`) lowers the detail with the element count (`lowest` above 1 M
  elements, `low` above 100 k, `medium` above 20 k), so large scenes at the defaults understate what is held: the
  geometry grows with the quality a user picks (for a figure, a closer look), and a mid-size structure at `auto` gets
  `medium` or `high`. Measure large scenes at `medium` or `high` (3j3q's cartoon at `medium` above), not only at the
  defaults.

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
  CPU half is what item 7 below takes up, for allocations rather than draws. `GraphicsContext.release` and `restore`
  (Prism's canvas presence, at most 12 live contexts) are the GPU side of residency.

## Goal

- A scene's CPU memory for geometry shrinks to what is in use: geometry lives on the CPU only while the purpose it
  was created for needs it; what is drawn lives on the GPU only after that.
- One rule for every kind, no heuristics: released as soon as it is not immediately needed (not after an idle time,
  not by a budget), configurable.
- Every reader gets its data, through one asynchronous request.
- Playback allocates nothing per frame once it has reached its high-water mark.
- CPU and GPU bytes are reported side by side.
- Geometry and render objects keep working without a graphics context, as today (item 11).

## Design

1. **Created for a purpose**: geometry is created by a visual's `createGeometry` for one of a few purposes, and is
   needed on the CPU only until that purpose is served:
   - **draw** (creating and updating a visual with a context): needed until it is uploaded and the readers that run at
     creation are done (the bounding sphere, `processValues` and CPU color smoothing, the position iterator a vertex
     or volume theme reads while the theme is computed, the counts kept as metadata);
   - **export** (geo-export, with a context or without one): needed until the exporter has written it. Geo-export
     walks the canvas's render objects (`controls.ts`), so the way to get released geometry back must be reachable
     from a render object (the visual registers it with the object), not only from the visual;
   - **another reader later** (a theme change reading positions, smoothing again, a props change that rebuilds from the
     data): it acquires (item 4) and is served like a creation for that purpose;
   - **restore** (a lost context, a backend switch): a creation for drawing.

   Each purpose ends at a known point (the upload, the end of the creating task, the end of the export or the
   reader's request), and the release happens there: no idle timers, no budget, no ordering by size. Step 1 lists the
   creation sites and the readers at creation for each geometry kind, so none is missed.
2. **Residency per geometry**: `resident` (CPU and GPU, as today) or `released` (GPU only, plus a way to get the data
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
3. **Releasing without a change**: the arrays in the geometry's cells are swapped for a shared empty array of the same
   type without a version bump (`ValueCell.releaseRef`, or the geometry kind's own `release`). So render items upload
   nothing, version-hashed caches (bounding spheres) keep their values, and draw counts stay.
   - The geometry is marked released (`Geometry.isResident`). The position iterator and anything else holding the arrays
     is dropped with it and created again on acquire.
   - A dev check: reading a released geometry outside an acquire throws. A test-mode proxy, not in production.
4. **Getting it back: `acquire`**: `Visual.acquireGeometry(): Promise<Geometry>` generalises
   `GeometryUtils.prepare` and `GeometryRequest`, for the geometry and the released layers. Two strategies; which one
   on-demand readers use while the context is alive is measured per kind and backend (step 3), not assumed:
   - **Read back** from the GPU (WebGPU `mapAsync`, WebGL2 `getBufferSubData` with a fence): exact, no CPU work, but
     a staging buffer and a wait for the GPU. Only while the context is alive.
   - **Recompute** by the recipe: the only way after a loss or a switch; also for props changes that rebuild from the
     data (the spheres' LOD, which then just recomputes).
   - Recomputed geometry must match what the GPU holds: same vertex and group counts, checked, and re-uploaded on a
     mismatch (non-deterministic producers are not eligible).
   - Acquired geometry is resident while the request that acquired it runs, and released when it ends (item 1).
5. **Readers through acquire**: geo-export, CPU color smoothing and `processValues`, position iterators for vertex
   themes, shape bounds, `getGroupCount`, clip volumes, the extensions' mesh utilities. Each is already asynchronous
   or runs in a task; the synchronous ones (`getGroupCount`, shape bounds) keep what they need as metadata at
   release (counts, group bounds), so they never acquire.
6. **Context loss and switches**: a released object is treated like one with foreign resources. `Scene.setContext`
   drops it and the visual recomputes in `prepareUpdate`. Nothing new on the restore path; it widens the existing one.
   A context released on purpose (`GraphicsContext.release`) with released geometry restores by recomputing,
   proportional to the scene (measure: Prism's canvas presence restores a background canvas this way).
7. **Rebuilt every frame without allocating**: two lifetimes, two granularities.
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

     While streaming, the geometry is immediately needed (the next frame builds into its buffers), so it is neither
     released (item 8) nor trimmed. When it clears, it is trimmed to its count and released like any other.
8. **Configurable**: `geometryResidency: 'keep' | 'release'`. `keep` is today's behaviour; `release` releases at the
   end of each purpose (item 1). No budget, no idle time, no ordering: the one rule is simpler to reason about and to
   test, and holds for every kind and scene size.
   - Configured where geometry is created, not only in Canvas3D: the visual context carries it (Canvas3D sets it from
     its prop, the plugin and Prism from their settings), so a caller without Canvas3D (an export in Node or a worker,
     item 11) sets it too.
   - Why configurable: an app that reads geometry often (exports repeatedly, smooths colors on every theme change, an
     editor) trades memory for not recomputing; a page with many canvases (Prism) or a large scene wants the memory.
   - `release` becomes the default once measured (step 5); until then `keep`.
   - Never released: streaming geometry (item 7), the marker layer (item 2), objects that are not eligible (no
     recipe: shapes from data the visual does not own; non-deterministic producers), and v1 merged members (their
     parts are rewritten from the members on layout changes).
9. **Accounting**: CPU bytes per render object (geometry and theme layers, from the values' arrays) and per context,
   next to the GPU counts in `stats`, the debug console and the benchmark. The WebGPU geometry pages are counted under
   `attribute`. Released and resident are counted separately, and so are the slack and the scratch arenas.
10. **Quick wins found on the way, done regardless of the policy**:
   - done: spheres and points keep center and group in one array (`positionGroup`, `center-group.ts`), which the
     shaders read (`tPositionGroup`): 16 instead of 32 CPU bytes per item. With levels of detail it is in the levels'
     order, reordered in place from the strides applied before (`orderPositionGroup`), so the order depends on the
     strides only; the bounding sphere is computed in the items' order. Shared with
     [primitive-expansion.md](primitive-expansion.md) (per-item data of the quad primitives);
   - left as it is: `Mesh.meta.originalData` (the index buffer before `uniformTriangleGroup`, for geometry export).
     The split only runs without flat interpolation (WebGL1 and a switch to it, `makeCompatible`), so it is absent
     on WebGL2 and WebGPU.
11. **Without a graphics context, as today**: representations, geometry and render objects are built and read with no
    context and no Canvas3D (export to glTF and the other geo-export formats in Node or a worker, tests, servers).
    Nothing here may change that:
    - without a context there is no GPU copy, so released geometry can only be recomputed: under `keep` (the default
      there) everything stays resident and `acquire` returns it at once; under `release`, set by the caller, geometry
      created for an export is released once it is written, and the next export recomputes;
    - nothing in a visual, a geometry or a render object starts to require a context or a residency manager (the policy
      and the accounting live in Canvas3D and the context, and are absent without them);
    - item 7 (views with slack, the scratch arena) is CPU only and works the same without a context.

## Steps

1. **Accounting and measurements** (the gate): CPU bytes per render object and per context in `stats` (item 9), the
   geometry pages counted. The creation sites per geometry kind with their purpose and the readers at creation (item
   1), so the end of each purpose is known. Measure:
   - CPU geometry and layer bytes against the JS heap for typical scenes (1cbs, 7qpd surfaces and cartoon, an assembly,
     the PetWorld model), large ones at `medium` or `high` quality as well as the defaults;
   - allocations and GC time per frame of trajectory playback (cartoon, spacefill, a molecular surface; a trajectory of
     a few thousand atoms and one of 100k), split into output arrays and intermediates by producer (allocation
     sampling in the DevTools heap profiler).

   Continue with release if geometry is a large share of the heap; with item 7 if playback allocates per frame
   measurably, its outputs part if the output arrays dominate, the scratch arena for the producers whose intermediates
   do.
2. The quick wins (item 10; done for spheres and points). Tests: builders give the same content (`builders.test.ts`),
   the levels' order depends on the strides only and the bounding sphere is the items' (`spheres.test.ts`), images
   unchanged on all profiles.
3. Release and acquire for meshes and spheres and their theme layers (items 1-5): the release at the end of each
   purpose, the released state, the version-neutral swap, the dev check, read back and recompute, the readers through
   acquire (geo-export from the render objects), metadata for the synchronous ones. Tests: released after the upload
   and the readers at creation, after an export, after an acquiring request; every reader on a released geometry or
   layer gets the same data as before, recompute equals the GPU counts, nothing is uploaded on release, marking works
   on released objects; `keep` changes nothing. Measured: read back against recompute per kind (mesh, spheres,
   layers), size (7qpd's surfaces, 3j3q's cartoon at `medium`, PetWorld's spheres) and backend (WebGPU, WebGL2),
   for the time to the data and the peak memory (the staging buffer, the intermediates); the faster one becomes the
   on-demand default per kind and backend.
4. Context loss and switches with released geometry (item 6). Tests: `renderer/context/release-restore` and
   `webgpu-recreate` with everything released, `pnpm test:viewer-smoke --switch` and `--fallback`. Item 11: geo-export's
   `prism-viewer.test.ts` (1cbs from the headless plugin's render objects, no WebGL, written as GLB, STL, OBJ and USDZ)
   keeps passing through every step; a case with `geometryResidency: 'release'` set and no context added to it (two
   exports, the second recomputing).
5. The setting (item 8) in the visual context, Canvas3D's prop and the plugin's and Prism's settings; Prism's canvas
   presence uses it. Measure memory of a Prism session with many canvases
   and the step 1 scenes under `release`, and the cost of acquires in normal use (export, theme changes, restore); with
   the gain shown and no regressions, `release` becomes the default.
6. Item 7: views with slack for meshes, spheres and cylinders (`ChunkedArray` compaction, the builders), the streaming
   mark by time, the trim when it clears; the scratch arena for the producers step 1 names. Measured with the playback
   cases of step 1: allocations per frame at zero after the high-water mark, frame time unchanged or better. Tests:
   sizes going up and down reuse the buffer, readers see the count (not the capacity), a cleared mark trims.
7. Later: cylinders, lines and text released; merged members; textures without their source (volumes and images
   recomputed or read back instead of `transferTexture`).
