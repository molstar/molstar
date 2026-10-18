# Plan: Marching cubes benchmark (CPU, WebGL, WebGPU)

## Status

- **Current**: done, steps 1–3; the benchmark also has color smoothing and render cases (`pnpm benchmark`).
- Step 1 done: `pnpm benchmark` (`scripts/benchmark.mjs`, page `tests/src/benchmark`, built with
  `scripts/esbuild/visual-tests.mjs --suite benchmark`); volume cases with the CPU and WebGL paths, baseline below. The
  WebGL stage times wait with a 1-pixel readback: the fence of `waitForGpuCommandsComplete` is polled and added about 90
  ms per wait.
- Step 2 done: Gaussian surface cases (protein-like balls of 1 k, 10 k and 100 k atoms at resolutions 1 and 0.5, the
  structure visual's density and marching cubes on each path) and end-to-end cases (the volume isosurface and Gaussian
  surface representations created from scratch, `tryUseGpu` off and on; a fresh volume per run, which caches its GPU
  texture). A full garbage collection before every run (`--expose-gc`): without it, CPU times depended on the garbage of
  earlier cases (up to 3× slower after the large volumes).
- Step 3 done: the WebGPU path (compute marching cubes on a storage buffer of the volume values, and on the 3d gaussian
  density) in all cases; results below.
- Color smoothing cases (`smoothing/gaussian/*`): the smoothing of the colors of a Gaussian surface on each path,
  results below.
- Render cases (`render/*`): frame times of scenes with a `Canvas3D`, on the backends `webgl2`, `webgl2-no-ubo`,
  `webgl1` and `webgpu` (`--backends`); results in [backend-selection.md](backend-selection.md#default-backend).
- Measures the isosurface extraction of the three paths, so the WebGPU port of [webgpu-compute.md](webgpu-compute.md)
  (step 3, marching cubes with compute shaders) can be compared with what it replaces: the CPU path
  (`computeMarchingCubesMesh`) and the WebGL path (`extractIsosurface`: active voxels, histogram pyramid, isosurface
  buffers). Related: [gpu-timer.md](gpu-timer.md) (GPU times), [resource-readiness.md](resource-readiness.md)
  (asynchronous readback).

## Goal

Numbers per path for the same inputs, on the same machine, reproducible with one command:

- time of the extraction alone (the input grid already on the device that runs it), and end to end (the representation
  from volume or structure to drawable geometry);
- first run (shader compilation, pipeline creation, allocations) separately from the steady state;
- triangle counts as a correctness guard (equal across paths for the same grid and iso level, up to vertex order).

## What is measured

| Path   | Extraction (stage times)                                                                 | Done when                                       |
| ------ | ---------------------------------------------------------------------------------------- | ----------------------------------------------- |
| CPU    | `computeMarchingCubesMesh` on the grid (one stage)                                       | the `Mesh` is returned                          |
| WebGL  | active voxels, histogram pyramid (with the synchronous sum readback), isosurface buffers | `waitForGpuCommandsComplete` after the last one |
| WebGPU | count, prefix sum, emit (compute kernels), then the asynchronous count readback          | the count is read back and the emit completed   |

- **Wall time** (`performance.now`) per stage and in total, with the GPU work completed before a stage's time is taken
  (`waitForGpuCommandsComplete`); the CPU work of the GPU paths (uploads, buffer setup) is part of it.
- **GPU time** per stage from timer groups where timestamps are available (WebGPU `timestamp-query`; WebGL
  `EXT_disjoint_timer_query_webgl2`, often unavailable in Chrome); reported, not required.
- **Repetitions**: 2 warm-up runs after the first, then the median and min of 7 runs; the first run reported separately.
- **Output**: vertex and triangle counts, the byte count of the result (texture mesh textures or mesh buffers).
- **End to end**: `createOrUpdate` of the volume isosurface representation and of the Gaussian surface representation
  with `tryUseGpu` on and off, until the render object is ready (and, on WebGPU, the count readback settled).

Not measured: rendering of the result (the same `mesh` program on all paths), and the Gaussian density on its own (a
separate, optional case below).

## Inputs

Deterministic, generated in the page or read from `data/examples`, no network:

- **Volumes**: the `blobVolume` generator of the visual tests (sum of Gaussian blobs) at 64³, 128³, 192³ and 256³ cells
  (16.8 M), at two iso levels each: a sparse surface (few active cells) and a dense one (many small surfaces, e.g. added
  high-frequency noise), since the compaction cost depends on the active fraction.
- **Gaussian surfaces**: synthetic structures (the `helixStructure` generator, repeated chains) of about 1 k, 10 k and
  100 k atoms at resolutions 1.0 and 0.5; optionally `data/examples/7qpd.fw2.cif` as a real structure.
- Grid sizes are capped per path where a path can't handle them (e.g. WebGL texture size limits for the 2d atlas): the
  case is reported as skipped for that path, not failed.

## Harness

- A benchmark page next to the visual tests (`tests/src/benchmark`, built like `tests/build/visual`) with the cases as
  functions per path, and a runner `scripts/benchmark.mjs` (Playwright, headless Chrome with the GPU, the same launch
  flags as the visual tests):
  - one page per backend (`cpu` runs in the WebGL page without the GPU paths), cases run **sequentially** in one page,
    one page at a time: no parallel pages, which would distort the times;
  - output: a table on the console (per case and path: first run, median, min, GPU time, triangles, speedup against CPU)
    and a JSON file (`--out <file>`) with the machine info (GPU renderer string or adapter info, Chrome version);
  - `--compare <file>` prints the change against an earlier JSON (for following the WebGPU port and regressions),
    `--filter` and `--paths cpu,webgl,webgpu` to select.
- Not part of `pnpm test` or the visual tests: timings depend on the machine and its load; run on demand.
- The existing manual page `tests/src/browser/marching-cubes.ts` (console timings of the WebGL stages) is superseded by
  it.

## Fairness

- The WebGL path stalls on its synchronous pyramid sum readback; that stall is part of its time. The WebGPU path's count
  readback is asynchronous; its time runs until the count is available, so both include the round trip.
- Inputs are uploaded before the extraction is timed (the grid texture or storage buffer), uploads timed as their own
  stage; the CPU path reads the `Float32Array` directly.
- Same iso level and grid; the paths' results differ only in vertex order (and float rounding of positions).
- The machine is otherwise idle; the runner warns when the first-run and median times of the CPU path differ by more
  than a set factor between two runs (noise check).

## Steps

1. Benchmark page and runner with the CPU and WebGL paths and the volume cases; baseline numbers recorded in this plan.
2. Gaussian surface and end-to-end cases (representations with `tryUseGpu` on and off).
3. The WebGPU path, as the marching cubes port lands ([webgpu-compute.md](webgpu-compute.md) step 3); the same cases,
   compared against the baseline.
4. Optional: the Gaussian density alone (CPU, WebGL splatting, WebGPU splatting), as input to step 2 numbers.

## Baseline

NVIDIA GeForce RTX 4070 Ti, ANGLE D3D11, Chrome; median of 7 runs in ms (in parentheses the first run where it is much
slower: shader compilation); WebGL with the representations' defaults (volume `gpuDataType` 'byte', density in half
floats). Triangle counts of the paths agree within 0.5% (equal on the gyroid). GPU times of a few ms vary by up to ±40%
between runs, and single larger cases occasionally by more (the WebGL `gyroid-192` isosurface stage measured 15 and 80
ms in two runs); compare WebGPU against several runs.

| Case                              | CPU  | WebGL     | WebGL vs CPU | CPU stages | WebGL stages          | Triangles |
| --------------------------------- | ---- | --------- | ------------ | ---------- | --------------------- | --------- |
| volume/blobs-64                   | 5.3  | 4.9 (296) | ×1.1         | 0.5 / 4.8  | 3.2 / 0.3 / 0.8 / 0.4 | 17 k      |
| volume/gyroid-64                  | 15   | 3.0       | ×5.1         | 0.9 / 15   | 1.4 / 0.4 / 0.8 / 0.5 | 157 k     |
| volume/blobs-128                  | 32   | 17        | ×1.8         | 2.0 / 30   | 16 / 0.4 / 0.7 / 0.5  | 68 k      |
| volume/gyroid-128                 | 89   | 12        | ×7.3         | 2.2 / 87   | 6.1 / 0.6 / 1.2 / 4.1 | 1.3 M     |
| volume/blobs-192                  | 87   | 56        | ×1.6         | 5.8 / 81   | 45 / 3.1 / 6.1 / 1.6  | 153 k     |
| volume/gyroid-192                 | 292  | 108       | ×2.7         | 6.8 / 285  | 22 / 2.1 / 4.3 / 80   | 4.4 M     |
| volume/blobs-256                  | 189  | 121       | ×1.6         | 15 / 173   | 90 / 7.3 / 21 / 2.6   | 272 k     |
| volume/gyroid-256                 | 681  | 86 (310)  | ×7.9         | 27 / 654   | 44 / 0.8 / 2.0 / 38   | 10.4 M    |
| gaussian/globule-1k-res1          | 7.9  | 2.6 (258) | ×3.0         | 4.3 / 4.0  | 1.5 / 1.0             | 8 k       |
| gaussian/globule-1k-res0.5        | 30   | 3.3       | ×9.1         | 18 / 14    | 2.4 / 1.0             | 33 k      |
| gaussian/globule-10k-res1         | 41   | 6.6       | ×6.3         | 30 / 12    | 5.5 / 1.1             | 36 k      |
| gaussian/globule-10k-res0.5       | 209  | 23        | ×9.1         | 145 / 65   | 21 / 2.6              | 146 k     |
| gaussian/globule-100k-res1        | 338  | 34 (267)  | ×10.0        | 272 / 65   | 32 / 1.3              | 163 k     |
| gaussian/globule-100k-res0.5      | 2028 | 76 (288)  | ×26.5        | 1564 / 464 | 73 / 3.2              | 660 k     |
| repr/isosurface-blobs-128         | 133  | 16        | ×8.3         | 133        | 16                    | 68 k      |
| repr/gaussian-surface-globule-10k | 138  | 14        | ×9.8         | 138        | 14                    | 130 k     |

Stages: volumes on the CPU setup / marching cubes, on WebGL upload / active voxels / pyramid / isosurface; Gaussian
surfaces density / marching cubes on both; representations one stage. The WebGL volume path is dominated by the upload
(building the 2d atlas on the CPU, then the texture upload), the CPU Gaussian surface by the density; the GPU density
makes the Gaussian surface ×6–26 faster on WebGL, the sparse volume surfaces gain least (×1.6–1.8).

## WebGPU (step 3)

Same machine, one run of all paths; median ms. WebGPU stages: volumes upload / marching cubes, Gaussian surfaces density
/ marching cubes (each including the asynchronous count readback). Triangle counts match the CPU's (volumes in f32:
equal on the blobs, within 0.004% on the gyroid, whose zero crossings sit exactly on grid points).

| Case                              | CPU  | WebGL | WebGPU | WebGPU vs WebGL | WebGPU stages |
| --------------------------------- | ---- | ----- | ------ | --------------- | ------------- |
| volume/blobs-64                   | 7.6  | 4.3   | 6.6    | ×0.7            | 0.5 / 6.1     |
| volume/gyroid-64                  | 15   | 2.5   | 7.8    | ×0.3            | 0.6 / 7.2     |
| volume/blobs-128                  | 30   | 17    | 10     | ×1.7            | 2.1 / 7.8     |
| volume/gyroid-128                 | 84   | 14    | 14     | ×1.0            | 2.1 / 12      |
| volume/blobs-192                  | 82   | 54    | 18     | ×3.0            | 8.8 / 9.4     |
| volume/gyroid-192                 | 283  | 109   | 32     | ×3.4            | 7.7 / 24      |
| volume/blobs-256                  | 187  | 113   | 57     | ×2.0            | 23 / 34       |
| volume/gyroid-256                 | 653  | 87    | 62     | ×1.4            | 15 / 47       |
| gaussian/globule-1k-res1          | 8.3  | 2.6   | 11     | ×0.2            | 4.4 / 6.5     |
| gaussian/globule-1k-res0.5        | 31   | 3.4   | 11     | ×0.3            | 4.9 / 6.0     |
| gaussian/globule-10k-res1         | 40   | 7.3   | 12     | ×0.6            | 5.8 / 6.0     |
| gaussian/globule-10k-res0.5       | 206  | 25    | 15     | ×1.6            | 7.4 / 7.5     |
| gaussian/globule-100k-res1        | 343  | 33    | 44     | ×0.8            | 36 / 8.3      |
| gaussian/globule-100k-res0.5      | 2058 | 76    | 113    | ×0.7            | 98 / 14       |
| repr/isosurface-blobs-128         | 128  | 16    | 9.4    | ×1.7            | 9.4           |
| repr/gaussian-surface-globule-10k | 137  | 15    | 20     | ×0.7            | 20            |

- Volumes: WebGPU is ×1.4–3.4 faster than WebGL from 128³ on, mostly by uploading the values directly into a storage
  buffer instead of building the 2d atlas on the CPU; small grids pay a fixed ~6 ms, the round trip of the count
  readback (indirect draws would remove it).
- Gaussian surfaces: the marching cubes are fast (6–14 ms), the 3d density is the cost on WebGPU (98 ms at 100 k atoms
  and 0.5 Å, against about 73 ms for WebGL's 2d density and marching cubes together): the slice-by-slice splatting into
  the 3d texture, a render pass per slice and pass type, is the next thing to measure and improve.
- The density with compute kernels instead (per atom and column of its splatting quad, the slices near the atom:
  fixed-point atomic sums of the density, the nearest atom by atomic max, then its group): Gaussian surfaces on WebGPU
  in median ms (density / marching cubes), triangle counts within 0.005% of WebGL's.

| Case                              | WebGL | WebGPU (render passes) | WebGPU (compute) | WebGPU vs WebGL |
| --------------------------------- | ----- | ---------------------- | ---------------- | --------------- |
| gaussian/globule-1k-res1          | 2.7   | 11                     | 8.8 (3.3 / 5.5)  | ×0.3            |
| gaussian/globule-1k-res0.5        | 3.2   | 11                     | 9.1 (3.5 / 5.9)  | ×0.4            |
| gaussian/globule-10k-res1         | 4.1   | 12                     | 10 (4.5 / 6.0)   | ×0.4            |
| gaussian/globule-10k-res0.5       | 22    | 15                     | 14 (5.8 / 8.0)   | ×1.5            |
| gaussian/globule-100k-res1        | 33    | 44                     | 22 (13 / 9.1)    | ×1.5            |
| gaussian/globule-100k-res0.5      | 76    | 113                    | 36 (22 / 15)     | ×2.1            |
| repr/gaussian-surface-globule-10k | 8.5   | 20                     | 14               | ×0.6            |

The small surfaces are bound by the fixed cost of the marching cubes count readback (about 6 ms).

## Current numbers

Measured again after the WebGPU work (compute density, indirect draws, render benchmark changes), same machine, median
ms; stages density / marching cubes.

| Case                                     | WebGL     | WebGPU    |
| ---------------------------------------- | --------- | --------- |
| gaussian/globule-1k-res1                 | 1.7 / 0.9 | 3.2 / 6.0 |
| gaussian/globule-10k-res0.5              | 6.7 / 1.2 | 5.8 / 8.5 |
| gaussian/globule-100k-res1               | 34 / 1.3  | 13 / 8.3  |
| gaussian/globule-100k-res0.5             | 95 / 3.5  | 22 / 14   |
| gaussian/1cbs-res1                       | 1.8 / 0.9 | 3.4 / 5.5 |
| gaussian/7qpd-res0.5                     | 12 / 2.2  | 10 / 16   |
| repr/gaussian-surface-globule-10k        | 7.3       | 8.6       |
| repr/gaussian-surface-globule-10k-update | 6.1       | 4.8       |
| repr/gaussian-surface-7qpd               | 40        | 38        |
| repr/gaussian-surface-7qpd-update        | 33        | 30        |

- The density is no longer the cost on WebGPU: faster than WebGL's from 10 k atoms at 0.5 Å on (×4 at 100 k atoms), 2–3
  ms on small inputs (the GPU round trip of the stage timing).
- The marching cubes stage of these cases extracts directly and waits for the count (about 6 ms at least, the readback
  round trip); the visuals draw indirectly into a capacity and don't wait, so the representations are at parity with
  WebGL or faster.

## Real structures

PDB entries of `data/examples` (`exampleStructure`), as the globules lack the cavities of real chains (2–15× fewer
vertices per atom, see [indirect-draws.md](indirect-draws.md)): 1cbs (about 1.2 k atoms, protein, ligand, waters) and
7qpd (about 11.5 k atoms in 7 units). Structure Gaussian surfaces (density / marching cubes, the count read back) and
the unit-wise representation (indirect draws), median ms:

| Case                              | CPU | WebGL | WebGPU | WebGPU stages |
| --------------------------------- | --- | ----- | ------ | ------------- |
| gaussian/1cbs-res1                | 14  | 3.0   | 9.1    | 3.2 / 5.6     |
| gaussian/1cbs-res0.5              | 63  | 4.8   | 12     | 4.2 / 7.9     |
| gaussian/7qpd-res1                | 108 | 6.0   | 16     | 6.7 / 8.9     |
| gaussian/7qpd-res0.5              | 601 | 12    | 24     | 9.9 / 14      |
| repr/gaussian-surface-7qpd        | 729 | 35    | 39     | 37 visible    |
| repr/gaussian-surface-7qpd-update | 713 | 31    | 29     | 28 visible    |

- The representation is at parity with WebGL; the structure surfaces pay the count readback of the direct extraction
  (6–14 ms of marching cubes against WebGL's 1–2 ms).
- Running these found two WebGPU paths falling back or coarsening silently, both from the WebGL 2d atlas criterion sized
  by the smaller WebGPU 2d texture limit: the representation took the CPU path for 7qpd at 0.5 Å (`suitableForCompute`
  now checks the 3d grid), and units coarsened their resolution (11% fewer vertices for a 2552 atom chain;
  `ensureComputeResolution` limits all cells instead). The check case `gpu-paths/gaussian-surface-units` compares the
  GPU surface of 7qpd per unit with the CPU one.

## Color smoothing

The Gaussian surface (resolution 0.5 Å) of the globules, a color per atom, smoothed on a 1 Å grid (the default
`resolutionFactor` 2) and a 0.5 Å grid, every 3rd vertex sampled. CPU: `calcMeshColorSmoothing` on the indexed mesh,
including the upload of the grid; WebGL: point sprites per slice into float targets, then a normalize pass; WebGPU: the
accumulate kernel (fixed-point atomic sums), then the normalize kernel into the atlas, in slabs of slices when the sums
exceed the storage buffer limit (the 100 k atom surface on the 0.5 Å grid: 22 M cells, 445 MB of sums). The GPU paths
sample every 9th vertex of the unindexed texture mesh, the same sample count. Same machine; median ms.

| Case                | CPU | WebGL | WebGPU | WebGPU vs WebGL | Grid texture |
| ------------------- | --- | ----- | ------ | --------------- | ------------ |
| globule-1k-res1     | 15  | 1.1   | 3.0    | ×0.4            | 0.3 MiB      |
| globule-1k-res0.5   | 21  | 1.5   | 3.5    | ×0.4            | 1.7 MiB      |
| globule-10k-res1    | 37  | 3.1   | 3.9    | ×0.8            | 1.7 MiB      |
| globule-10k-res0.5  | 91  | 17    | 5.0    | ×3.4            | 11.8 MiB     |
| globule-100k-res1   | 192 | 21    | 12     | ×1.7            | 11.7 MiB     |
| globule-100k-res0.5 | 577 | 37    | 8.9    | ×4.2            | 84.9 MiB     |

- Both GPU paths are ×5–65 faster than the CPU. WebGPU has a floor of about 3 ms (buffers, bind groups and the
  submission per call) and wins on the larger grids, where WebGL draws every sample once per slice.
- The 1 Å grids are slower than the 0.5 Å ones at 100 k atoms on WebGPU: more samples per cell contend for the same
  atomics.

## Decisions

1. In the repository (`tests/src/benchmark`, `scripts/benchmark.mjs`), run on demand.
2. Sizes up to 256³ cells on all paths; larger only on WebGPU, later if needed.
3. Step 1 lands first and records the baseline the WebGPU port is compared with.
