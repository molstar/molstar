# Plan: GPU surface buffers (WebGPU)

## Status

- **Current**: steps 1–3 done; decisions 1A, 2A, 4A (4 K vertex chunks) and 5C (several arenas, a surface in one) made, 3
  superseded by chunked storage. Step 4 postponed: multi-draw indirect only replaces the per-range draw loop of
  render items with multi-draw lists (culling, level of detail of instanced structures), needs
  `indirect-first-instance` for its fallback, and few adapters offer it (as `WEBGL_multi_draw` on WebGL, which still
  helps real cases); to measure first (draws of multi-draw lists per frame in an assembly scene), exact GPU counts per
  range for chunked surfaces being the part independent of the feature.
- After step 3: the WebGPU texture output removed (`computeIsosurface`, the kernels `mc-count`, `mc-emit` and
  `mc-finish`, the texture target of `computeTextureMeshIsosurface`, the draw arguments of `TextureMesh.DoubleBuffer`,
  which keeps its textures for WebGL); compute color smoothing reads only the arena (its texture variant had no other
  user); the benchmark's volume and Gaussian surface cases extract into the arena (until the count is read back), its
  WebGPU color smoothing cases smooth a surface extracted into the arena from the GPU density (not timed; triangles as
  the GPU extracts them); the Dawn test compares the count with the CPU marching cubes. The prefix sum (`gl/compute/prefix-sum.ts`) is no longer
  used by marching cubes, kept as a general kernel. Stats: `ResourceByteCounts.storageBuffer` (WebGPU, the bytes of the
  storage buffers created and not destroyed, `ComputeContext.getStorageBufferByteCount`; WebGPU still reports 0 for the
  other resources), the byte count of a chunked texture mesh's renderable includes its pages and page table (read back,
  `meta.chunkByteCount`), and the debug console stats show the arena (`getVertexArenaStats`: used and allocated pages).
- Step 3 done: render schemas declare read-only storage buffers (`StorageBufferSpec` with `usedIf`, WebGL ignores
  them; WebGPU numbers them after the draw uniforms of the object group, bound when referenced); the vertex arena
  (`gl/compute/marching-cubes/chunks.ts`: one growable arena of 4 K vertex pages per context, 12 B per vertex since
  the compact format follow-up, a free page stack) with the kernels `mc-chunk-count`
  (one atomic reservation per workgroup), `mc-chunk-allocate`, `mc-chunk-emit`, `mc-chunk-free` and `mc-chunk-push`;
  texture meshes in the arena (`TextureMesh.createChunked`, `dIsoChunks`, the WGSL mesh program reading through the
  page table, vertices beyond the drawn count clipped), the color smoothing kernel and the readback (position
  iterator, geo export) reading through it. `computeTextureMeshIsosurface` extracts into the arena where the kernels
  are supported (`IsosurfaceTarget`, `extractIsosurface`); the count read back afterwards raises the mesh's upper bound
  (`vertexCount`, for draws without the GPU count and the samples of derived data, whose handlers then run again) and a
  surface missing pages is extracted again after the arena grew. Tests on Dawn (`webgpu-isosurface-chunks.test.ts`:
  count as the exact path, vertices on the surface, pages reused and returned, missing pages after growth); visual tests
  (all profiles), smoke and unit tests pass. WebGPU, median of 25 runs, after step 2 → after step 3: 7qpd created 15 →
  6.6 ms, updated 9.8 → 7.8 ms; globule-10k created 8.6 → 5.3 ms, updated 5.4 → 4.5 ms; blobs-128 created 8.7 → 6.0 ms
  (no readback awaited for new volume surfaces), updated 3.4 ms (unchanged). Memory of 7qpd (7 polymer units, 3.69 M
  vertices): 902 pages, 89 MB used (0.3% page slack), in an arena of 1520 pages (143 MB, grown ×1.5 from 64 pages); the
  textures before were 28 B per vertex of a 1.25× capacity, about 130 MB, twice that after an update (double buffer).
- Step 2 done: `getScratchBuffer` (`gl/compute/scratch.ts`) keeps a storage buffer per context and key, grown,
  released after 3 s unused; the density scratch buffers use it (cleared per slab) and the marching cubes `counts` and
  `offsets` with a capacity (without one, the count is awaited between the kernels, when another extraction could use
  them); the shared density texture of the Gaussian surface only grows (`growTexture`, read by grid coordinates, reads
  clamped to the grid). Keeping output textures up to 2× was dropped: no gain, as the double buffer already redefines
  each set to the size it had. WebGPU, median of 25 runs: 7qpd update 34 → 9.8 ms, created about 40 → 15 ms,
  globule-10k update unchanged, blobs-128 update 3.0 → 3.4–3.6 ms (outside the GPU passes, which got faster; not
  explained). Visual tests (all profiles), smoke and unit tests pass.
- Step 1 done: `pnpm benchmark --gpu-timing` times the passes of the representation cases (GPU time per pass label,
  summed over the run); the GPU passes are a fraction of the wall time, the rest is resource creation in the GPU process
  (see Measured), so the pool moved before the storage buffers (step 2).
- Follows [indirect-draws.md](indirect-draws.md) (capacity, indirect draws,
  overflow re-extraction); builds on [storage-buffers.md](storage-buffers.md) (WebGPU core, vertex pulling from storage
  buffers) and shares multi-draw indirect with layer C of [webgpu-draw-submission.md](webgpu-draw-submission.md)
  (step 8).

## Question

Marching cubes on WebGPU write their vertices into three textures (position `rgba32float`, normal `rgba16float`, group
`rgba8unorm`) sized for a capacity: an estimate for new Gaussian surfaces, exact for new volume surfaces, 1.25× the last
count on updates; the count goes into draw arguments and is read back without waiting, a count above the capacity
extracts again. The output size varies with the surface (iso value, smoothness, resolution, the units of a structure).
Should the output instead go into storage buffers allocated in chunks on the GPU, so that a surface never has to be
extracted again for lack of space, never holds much more than it uses, and the CPU never waits for (or needs) its
vertex count; drawn with multi-draw indirect where Chrome offers it?

## Measured (7qpd, 21 Gaussian surface units, NVIDIA, Windows)

- GPU passes are a fraction of the time: of the 32 ms update of 7qpd (7 polymer units extracted), the passes take 6.7
  ms (density 1.8, count with the prefix sum 3.9, emit 0.9, finish 0.1); the main thread records all units in about 7
  ms and is idle otherwise. The rest is spent in the GPU process, outside the timestamps and the page's CPU profile.
- It is resource creation: each extraction creates large resources and drops them again (the density's three grid-sized
  scratch buffers, the shared 3d density texture defined to each unit's grid, the `counts` and `offsets` of marching
  cubes, 4 B per cell each, and output textures whenever the capacity changes). Kept and grown instead (an experiment),
  in steps:

  | Kept                                     | 7qpd update | 7qpd created | blobs-128 update | globule-10k update |
  | ---------------------------------------- | ----------- | ------------ | ---------------- | ------------------ |
  | none (today)                             | 32 ms       | 40 ms        | 3.0 ms           | 5.5 ms             |
  | density scratch buffers                  | 23 ms       | 36 ms        |                  |                    |
  | + marching cubes scratch, output ≤ 2×    | 15 ms       | 26 ms        | 4.0 ms           | 5.7 ms             |
  | + shared density texture, grown only     | 12 ms       | 19 ms        | 4.0 ms           | 6.0 ms             |

  Disabling Dawn's lazy clears (`--disable-dawn-features=lazy_clear_resource_on_first_use`) changed nothing, which points
  at the allocation itself (D3D12 resources of tens of MB) rather than its zero-initialization (not verified that Chrome
  applied the toggle). Single-surface cases gain nothing (blobs-128 slightly slower, to be checked: clears of larger
  kept buffers). - Memory: capacities are 1.25× the count (13.8 M vertices for 11.1 M, at 28 B per vertex about 390 MB
  of textures, 77 MB of it headroom); 3 of 21 units overflowed their estimate (small units), none of the updates.
  Texture meshes keep two texture sets (`TextureMesh.DoubleBuffer`), so an updated surface holds about twice its size
  (about 650 MB for 7qpd).

What this plan can gain: time (resources kept instead of created per extraction), memory (packed vertices, one set
instead of two, pooled capacity), one emit pass instead of count + prefix sum + emit (GPU time, the scratch buffers),
and fewer draw calls with multi-draw indirect. Mol\* extracts on updates, not per frame, so a ring buffer of per-frame
output does not fit: surfaces keep their vertices until the next update.

## Design

- **Vertex buffers on WebGPU**: the extraction writes storage buffers, read in the vertex stage by vertex id (as the
  textures are read today). WebGL keeps its textures (no storage buffers); the texture mesh holds either, behind the
  neutral geometry (`TextureMesh` with a buffer variant, or a separate GPU mesh kind, see decision 1).
- **Layout**: position and group in one `vec4<f32>` (group bitcast into `w`) and the normal in 8 B
  (`pack2x16snorm` × 2, or 4 B octahedral): 24 or 20 B per vertex instead of 28. One vertex buffer of 16 B entries at
  the default `maxStorageBufferBindingSize` (128 MiB) holds 8.4 M vertices, the cap of the Gaussian estimate today;
  larger surfaces request the adapter's limit and otherwise split (as the textures are limited by
  `maxTextureDimension2D` today).
- **Chunked vertex storage**: vertices live in fixed-size chunks (pages) of one arena of storage buffers per graphics
  context (decision 5); a surface owns a page table (a small buffer: its vertex count, page count and page ids) and its
  vertex `v` is at `page[v / ChunkSize] * ChunkSize + v % ChunkSize` of the arena. Growing replaces the arena's
  buffers by larger ones with their contents (page ids stay valid), the meshes bind the new ones (arena listeners);
  the arena is at most `maxStorageBufferBindingSize`. Waste per surface is less than one chunk; no capacity, no size
  classes.
- **Allocation on the GPU, without readback**: the arenas' free pages are a stack in a pool buffer (count and page ids).
  An extraction runs three kernels in one submission:
  1. **count**: per workgroup, the vertex count of its cells (workgroup memory), reserved with one `atomicAdd` on the
     surface's vertex count; the workgroup's base is stored (4 B per workgroup, not per cell);
  2. **allocate** (one workgroup): the pages the count needs, popped from the free stack into the page table (surplus
     pages of the previous extraction pushed back, so updates reuse their pages); writes the draw arguments (the count,
     clamped to the pages it got);
  3. **emit**: per workgroup again, its vertices from its base, written through the page table.

  The vertex order varies between runs (decision 2). Destroying a surface runs a free kernel pushing its pages back.
  Kernels touching the free stack run in queue order (single workgroup), so it needs no locking.
- **Pool growth, off the critical path**: the free page count is read back asynchronously after extractions (never
  awaited); the CPU adds an arena (its pages pushed by a kernel) when the free pages fall below a headroom, and before
  an extraction whose estimate (atoms and resolution for Gaussian surfaces, a fraction of the cells for volumes)
  exceeds the free pages last known minus the estimates in flight. Should the stack still run dry during an
  extraction, the allocate kernel clamps the count and sets a status word; read back asynchronously, it makes the CPU
  grow the pool and extract that surface again: the only re-extraction left, a safety net rather than the path.
- **Readbacks are bookkeeping**: the vertex count (stats, position iterator, export) and the pool's free count arrive
  asynchronously; nothing waits for them, and no allocation size depends on them.
- **Scratch resources**: the per-extraction buffers (density accumulation; marching cubes counts and offsets until the
  chunked path replaces them) and the
  shared 3d density texture are kept per context and grown, not created per call (shrunk when unused for a while); the
  density texture only grows where its readers address it by the grid size (the marching cubes input), not where it is
  sampled by texture coordinates (the density volume visual keeps an exact texture).
- **One storage per surface on WebGPU**: an update writes after the frames that drew the previous vertices (queue
  order), so the double buffer is not needed (WebGL keeps it).
- **Multi-draw indirect** (`chromium-experimental-multi-draw-indirect`, optional, requested where the adapter has it):
  draws the multi-draw lists (culled and level-of-detail ranges) from an arguments buffer with one `multiDrawIndirect`,
  instead of one draw per range; the fallback is the loop of draws (core). For GPU surfaces in multi-draw lists the
  per-range counts can then come from the GPU (decision 2B of indirect-draws.md, today the capacity is drawn). With all
  GPU surfaces in the same arenas, surfaces sharing a pipeline could later be drawn by one call (per-draw data indexed
  by draw). Not in the TypeScript DOM types: declared locally. Few adapters report it so far (web3dsurvey: about 0.25%
  on Windows), so the core path stays the measured one.

## Decisions

### 1. Where the buffer variant lives

- **A. `TextureMesh` with buffers on WebGPU**: the same geometry kind and visuals; values and shader chunk pick
  textures or buffers by backend (a define). Consumers (color smoothing, the readback for the position iterator and
  export) branch on the storage.
- **B. A GPU mesh kind**: a new geometry for compute-produced vertices on WebGPU, `TextureMesh` stays WebGL's. Cleaner
  types, but every visual creating texture meshes chooses between two geometries.

Chosen: A (one geometry for GPU surfaces; the storage is a backend detail).

### 2. Emission order

- **A. Atomic counter**: one pass, no scratch buffers; the vertex order varies between runs (workgroups reserve in the
  order they run). Images are the same where order does not matter (depth-tested opaque surfaces, order-independent
  transparency, fixed-point color smoothing sums); blended transparency without OIT and exact z-fighting ties may
  differ between runs.
- **B. Count and prefix sum, into buffers**: today's order (deterministic), only the output changes.

Chosen: A, directly (no prefix-sum step in between); blended transparency without OIT may look different between
runs.

### 3. Capacity classes

Superseded by the chunked storage (earlier choice: ×1.5 steps).

### 4. Chunk size

- **A. 4 K vertices** (48 KB at 12 B): little waste for small units (waters, ligands), page tables of a few thousand
  entries for the largest surfaces (11 M vertices: 2.7 K pages).
- **B. 16 K vertices** (384 KB): smaller tables, up to 384 KB unused per surface (21 units of 7qpd: up to 8 MB).
- **C. 64 K vertices** (1.5 MB): fewest pages, up to 1.5 MB unused per surface.

Chosen: A (the table lookup per vertex is the same for all; small units are common).

### 5. Arenas

- **A. One growable arena**: one buffer per vertex attribute holding all pages; growing allocates larger buffers and
  copies on the GPU (rare, off the critical path). Fewest bindings (page table and the arena buffer), simplest
  shaders; the total is limited by the storage buffer binding size (11 M vertices at the 128 MiB default).
- **B. Several fixed arenas**: growing adds a buffer without copying; shaders select the arena with a `switch` (no
  binding arrays in WebGPU), two bindings per arena, about 3 arenas within the default 8 storage buffers per stage.

Chosen first: A (more arenas only if the binding size is ever reached). Then C, when 3j3q at medium quality needed
more than a 2 GiB buffer (see Follow-ups):

- **C. Several arenas, a surface in one**: as B, but a surface's pages are all in one arena and what draws it binds
  that one buffer (no `switch`, one binding as A). The first arena grows (copying) to 256 MiB, more are added at
  256 MiB (larger for a surface of more pages, up to the storage buffer limit) without copying; a surface whose arena
  can't take its pages moves (its pages freed, taken from the other arena).

## Steps

1. Measure the GPU passes of the extraction (density, count, prefix sum, emit, finish) per unit with timestamp
   queries. Done (results in Measured).
2. Keep the scratch resources: density scratch buffers and the marching cubes counts and offsets per context, grown;
   the shared density texture grown only. Done (results in Status).
3. Chunked vertex storage on WebGPU: the arena, the free stack and page tables; the count, allocate and emit kernels
   (replacing the prefix sum and the marching cubes scratch) and the free kernel; the WGSL texture-mesh chunk reads
   through the page table, color smoothing and the readback too; one storage per surface (no double buffer on WebGPU);
   pool growth from asynchronous free counts and estimates, the status word as safety net. Visual tests, benchmark and
   memory on 7qpd (blended transparency may differ between runs). Done (results in Status).
4. Multi-draw indirect where available: multi-draw lists from an arguments buffer, the loop as fallback; GPU counts for
   GPU surfaces in multi-draw lists. Benchmark on an instanced scene with culling (`render/*`).

## Follow-ups

- **Compact vertices, the arena's limit reported** (done): 3j3q (1356 chains, 2.44 M atoms) at medium quality
  (resolution 0.8) has about 124 vertices per atom, 303 M vertices: 6.8 GiB at 24 B, the arena capped at a 2 GiB
  storage buffer (2 GiB of positions, 1 GiB of normals: 134 M vertices) and surfaces past it drawn incomplete
  without notice. Vertices are 12 B now, one buffer (`iso_vertex_wgsl`): the position quantized to 16 bits per axis in
  the box of the grid's vertices (its min and step in the page table header, written by `mc-chunk-allocate` in queue
  order; extent / 65535, about 0.002 Å for a 100 Å unit), the normal octahedral in 24 bits, the group as packed rgb
  (24 bits, as before); 179 M vertices fit a 2 GiB arena. Page tables are sized to the estimate (twice it and 16
  pages, not the arena's pages: 128 KiB per surface at 2 GiB), a surface beyond its table is extracted again with its
  count. A surface still missing pages with the arena at its limit is not extracted again but reported, all of them
  in one message (`GraphicsResources.outOfMemory`, `Canvas3D.outOfMemory`, logged as an error by the plugin: Prism's
  notifications), as WebGPU out-of-memory errors (once per 5 s).
- **Several arenas** (done, decision 5C): the pages of an extraction are reserved before it
  (`reserveVertexArenaPages`, in turn per context): an arena with the free pages, else one grown, else a new one, its
  buffers created in an out-of-memory error scope and awaited (`ComputeContext.tryCreateBuffer`): a buffer that failed
  invalidates the command buffer using it, the other surfaces' work and the frame with it. So extractions await before
  recording and compute shared inputs themselves: the Gaussian surface computes its density in each extraction
  (`IsosurfaceExtraction`). Without memory for another arena, the arenas are exhausted until a surface returns its pages;
  the surfaces that did not fit are reported.
- **Indexed surfaces** (done): a vertex per crossed edge of the grid, owned by the edge's lower point and shared by
  the cells around it (as the CPU marching cubes, whose vertex counts the Dawn test matches), a triangle the 3 indices
  of its vertices; a page holds 4 K vertices or 4 K triangles (12 B each), the table lists vertex pages, then triangle
  pages. `mc-chunk-count` runs per grid point: its crossed edges (a bit per axis) and first vertex in a scratch word
  per point (`pointBases`, 4 B a point, as the old per-cell scratch), the cell's triangles per workgroup;
  `mc-chunk-emit` writes a point's vertices and its cell's triangles, a corner's vertex the base of the edge's point
  and the crossed edges before it (`countOneBits`, no value reads). Draws pull through the indices (`isoSlot`, 3 per
  triangle, non-indexed draws), so the mesh and color smoothing programs are unchanged; readbacks expand to the
  unindexed layout. A vertex's group is of the nearer edge end (else the other), a triangle's of its first vertex
  (flat), as the CPU path; before, a group per cell. 18 B a triangle against 36 B: 3j3q at medium about 1.7 GiB.
  The kernels compile in about 150–250 ms on Dawn against 90 ms (the emit kernel's edges in a loop of the set bits,
  not unrolled per axis: 1 s unrolled).
- **Empty arenas released** (done): an arena without surfaces (counted, not by its free pages, which a readback still
  pending leaves too high) and without reservations is released after 3 s (`ArenaIdleMs`), the first too (a 64-page
  one is created again when needed); a reservation or a surface in between keeps it. Longer than a frame of a
  trajectory played back, whose surfaces may be created again each frame (the new one's pages taken before the old
  one's are freed): its arenas are kept while playing (a Dawn test of 5 frames) and released when it stops. New
  surfaces fill the first arenas that fit, so the later ones drain.

- **Size estimates for volume isosurfaces** (done): volume isosurfaces and segments pass an estimate from the grid
  (`estimateIsosurfaceVertexCount`): the triangles of the cells by `TriTable`, classified as the compute kernels'
  `cubeCase`, so exact where the grid is visited whole (up to 2^20 cells, about 6 ms; 12.5 ms at 2^21), larger grids by
  every n-th z-slice, scaled, with 20% headroom; computed only for a mesh without a last count. `1tqn.ccp4` (80x100x144)
  at its default level: before, the first extraction missed 204 pages (445 in the browser) and was extracted again;
  with the estimate, once. Earlier, that path bound a destroyed page table, fixed (`TextureMesh.bindChunkedTable`).
  Tests: the estimate against the CPU marching cubes (exact, strided, sampled), on Dawn the two spheres that missed
  pages without an estimate get all with it.
- **CPU or GPU for small Gaussian surfaces** (done): measured (`mc/gaussian-cutoff/*`, globules of 100 to 30k atoms at
  resolutions 0.5 to 3; CPU, WebGL and WebGPU paths) on an RTX 4070 Ti and on the integrated AMD Radeon of a Ryzen
  7900X (2 CUs, `--browser-arg --force_low_power_gpu`). WebGL: the cutoff at resolution 1 holds on the integrated GPU
  (above it, the histogram pyramid is about as fast as the CPU or slower: 30k atoms at 1.5 in 70 ms, the CPU 55 ms), kept.
  WebGPU: about 6–7 ms whatever the size (dispatches, the count's readback), faster than the CPU from about 100k
  density cells at any resolution, on the integrated GPU too (30k atoms at resolution 2: 11 ms, the CPU 30 ms; at 1:
  21 ms, the CPU 124 ms); below, the CPU takes 1–6 ms. So the compute path decides by size instead of the resolution:
  the CPU for density grids under 100k cells of structures under 5000 atoms (the atoms bound the CPU's density work on
  coarse grids: 10k atoms at resolution 3, 43k cells, 10 ms). The surfaces of both paths have the same triangles
  within 0.1% at all resolutions; the visual cases `gaussian-surface-gpu-res2` and `-res3` compare WebGPU's GPU surface
  of a 6k-atom ball with WebGL's CPU surface.
