# Plan: Indirect draws for GPU-produced geometry (WebGPU)

## Status

- **Current**: done, steps 1–5.
- Steps 1–5 done (decisions 1C, 2A, 3A). A neutral `drawIndirect` render item value (WebGPU `drawIndirect` without
  multi-draws); marching cubes with a capacity (`mc-finish`: the count into the draw arguments, a degenerate tail);
  `computeTextureMeshIsosurface` picks the capacity (exact for new volume surfaces, an estimate for new Gaussian
  surfaces, 1.25× the last count on updates) and extracts again on overflow, running the re-extraction handlers of the
  mesh (texture mesh smoothing registers one per kind, as it samples the vertices); the vertex count of texture meshes
  is one cell shared by `drawCount` and `uVertexCount`. Visual cases `draw-indirect`, `isosurface-gpu-update` and
  `isosurface-gpu-update-overflow`; the smoothed Gaussian cases overflow their estimate (thin helices) and cover the
  handlers. Results below. Step 5 of [webgpu-compute.md](webgpu-compute.md) and step 6 of
  [webgpu-compute-shaders.md](webgpu-compute-shaders.md) (counts read back first, indirect draws as the follow-up).
- Motivation, from [marching-cubes-benchmark.md](marching-cubes-benchmark.md): WebGPU marching cubes pay a fixed ~6 ms
  for the round trip of the vertex count (`readAsync` of the prefix sum total) before the emit kernel can run into
  textures of the right size. It dominates small and medium surfaces (the 10 k atom Gaussian surface representation: 14
  ms on WebGPU, 8.5 ms on WebGL) and every interactive update (e.g. dragging the iso value).

## Where the CPU uses the count today

| Use                                                                 | Where                                         | Needs the exact count?                                   |
| ------------------------------------------------------------------- | --------------------------------------------- | -------------------------------------------------------- |
| Size of the output textures (position, normal, group)               | `computeIsosurface`                           | an upper bound suffices                                  |
| `drawCount` of the render item                                      | `TextureMesh.createValues`, render item       | no, with indirect draws                                  |
| Culled / level-of-detail multi-draws (per-draw counts from the CPU) | `renderable.ts` mdb lists, WebGPU render item | no, if the tail beyond the count is degenerate           |
| `uVertexCount` (vertex-instance color index)                        | `varyings.wgsl.ts`                            | an upper bound suffices (a stride, consistent with data) |
| Color smoothing sample count                                        | `calcTextureMeshColorSmoothing`               | no, if the kernel reads the count from the GPU buffer    |
| Vertex-granularity colors and the position iterator                 | `createColors`, `createPositionIterator`      | WebGL only today (reads the textures back)               |
| Stats, export, bounding data                                        | stats; bounding spheres come from the input   | stats only                                               |

## Design

- **Neutral render item value**: an optional `drawIndirect` value (a `StorageBuffer` holding draw arguments) next to
  `drawCount`. WebGPU draws with `drawIndirect` when it is set and there is no multi-draw list; WebGL never sets it (no
  storage buffers). The arguments buffer holds `vertexCount, instanceCount, firstVertex, firstInstance`; the count is
  written by a kernel (copied from the prefix sum total, clamped to the capacity), the instance count by the CPU when it
  changes.
- **Capacity**: the output textures have a capacity of at least the vertex count; the emit kernel skips vertices beyond
  it and writes degenerate (zero) positions into the tail, so drawing the full capacity is always safe (multi-draw
  paths, and any consumer that only knows the capacity). `TextureMesh` gets `vertexCount` (exact once known) and a
  `vertexCapacity`, and an optional GPU count buffer.
- **The count still comes back, asynchronously**: `readAsync` of the total is started but not awaited by the geometry.
  When it resolves, the exact count updates the texture mesh (stats, CPU-side consumers); if it exceeds the capacity the
  surface is computed again with the exact capacity (tracked as pending work, so screenshots wait for it).
- **Direct draws once the count is back** (`TextureMesh.setDrawnVertexCount`): the count read back goes into the draw
  arguments (`DrawIndirect.count`), and the WebGPU render item draws it directly until the next extraction; `drawCount`
  and `uVertexCount` stay the upper bound (the stride of vertex-instance data). An indirect draw costs far more GPU time
  than a direct one in Chrome (Apple M4, without render bundles): 3j3q in the Prism viewer, 1356 Gaussian surfaces drawn
  about 3× a frame, took 105–540 ms a frame indirectly and 20–25 ms directly, with the same CPU time (WebGL2: about 40
  ms, its CPU surfaces merged). Indirect draws stay for the frames between an extraction and its readback. Render
  bundles (layer C of [webgpu-draw-submission.md](webgpu-draw-submission.md)) recorded with such a draw must be
  recorded again when the `drawIndirect` cell changes: the count changes twice per extraction (cleared, then read back)
  without a change of `drawCount` or the instance count.

## Decisions

### 1. Capacity without knowing the count

- **A. First computation reads back, updates go indirect**: a new texture mesh is sized exactly (as today); updates of
  an existing one reuse its capacity (with headroom, e.g. 1.25× the last count) and draw indirectly, growing on
  overflow. Interactive updates lose the round trip, creating a surface does not.
- **B. Always indirect, with an estimate**: the first capacity comes from an estimate (Gaussian surfaces: atoms and
  resolution; volumes: a fraction of the cells), overflow recomputes. Removes the round trip everywhere, but a bad
  estimate costs a second computation (dense volumes such as the gyroid) or memory.
- **C. A, plus estimates where they are reliable**: A for volumes, B for Gaussian surfaces (whose vertex count scales
  with atoms / resolution²).

Chosen: C. Volumes as A; Gaussian surfaces estimate their first capacity (see step 4).

### 2. Multi-draws (culling, level of detail)

- **A. Draw the capacity**: per-draw counts are the capacity; the degenerate tail draws nothing.
- **B. Indirect arguments per draw**: a GPU copy of the count into each draw's arguments, every frame the list changes.

Chosen: A (B only if the tail shows up in profiles).

### 3. Color smoothing on indirectly drawn meshes

- **A. Count from the GPU buffer**: the accumulate kernel reads the vertex count from the count buffer (dispatched for
  the capacity).
- **B. Wait for the readback**: smoothing runs once the exact count is known.

Chosen: A.

## Steps

1. Remove the stale wait in `createComputeGaussianSurface` (`pending.whenReady()` before the extraction, from when the
   density rendered again after its pipelines compiled).
2. Neutral `drawIndirect` render item value; WebGPU `drawIndirect`; a compute visual case drawing a kernel-written
   count.
3. Marching cubes with a capacity: degenerate tail, the count copied into draw arguments, the total read back without
   being awaited; `TextureMesh` capacity, exact count updated when it arrives, recompute on overflow.
4. Volume isosurface and Gaussian surface visuals on WebGPU use it: volumes sized exactly when created and reused (with
   headroom) on updates, Gaussian surfaces from an estimate (atoms and resolution, calibrated on the benchmark
   globules); color smoothing reads the GPU count.
5. Benchmark: the update path (a second extraction into the same textures) next to the creation path.

## Results

Same machine as [marching-cubes-benchmark.md](marching-cubes-benchmark.md), median ms; "visible" until the GPU work
completed, the rest until pending work settled (the count read back, extracted again on overflow).

| Case                                                  | WebGL | WebGPU before | WebGPU | WebGPU visible |
| ----------------------------------------------------- | ----- | ------------- | ------ | -------------- |
| repr/isosurface-blobs-128 (created, exact)            | 17    | 9.4           | 9.2    | 9.1            |
| repr/gaussian-surface-globule-10k (created, estimate) | 6.6   | 14            | 5.9    | 3.0            |
| repr/isosurface-blobs-128-update                      | 2.8   |               | 3.1    | 2.9            |
| repr/gaussian-surface-globule-10k-update              | 20    |               | 4.9    | 4.8            |

The triangle counts of indirectly drawn surfaces in the benchmark output are their capacity (the drawn count).

## Gaussian surface estimate

Calibrated on the CPU Gaussian surfaces of the units of PDB entries in `data/examples` (1crn, 1cbs, 1grm, 7qpd, 1bna) at
resolutions 0.5 and 1 Å, unindexed vertices (three per triangle) · resolution² / atoms:

| Units                                                    | Vertices · res² / atom  |
| -------------------------------------------------------- | ----------------------- |
| Protein and DNA chains (80–3365 atoms)                   | 73–86                   |
| Ligands, short peptides                                  | 98–133                  |
| Waters (isolated atoms)                                  | about 240               |
| Benchmark globules (dense jittered lattice, no cavities) | 5–39, falling with size |

Real chains are linear in the atoms (internal cavities), unlike the compact globules the first estimate (∝ atoms^⅔) was
fitted to, which underestimated chains about 3×. The estimate is 100 · atoms / resolution² (headroom over the polymers),
at most 8 M vertices (about 224 MB of textures); waters and ligands overflow it and are extracted again (small units).
The `gpu-paths/representations` check asserts that the GPU paths are taken with support as reported without strict mode
(the strict backend of the visual tests reports every program as supported, which hid a fallback to the CPU).
