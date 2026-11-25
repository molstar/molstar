# Plan: The sphere raster's work shared in subgroups and workgroups

## Status

- **Done**: steps 1, 2 and 4 (measured; chunks packed; pixels balanced: Results); steps 3 and 5 measured and dropped
  (persistent workgroups, kernel variants). Left: none planned; the zoomed cases' cost is the pixels' ray tests and
  atomics, the yeast cell's its 73k workgroups' reads.
  Follows [sphere-compute-raster.md](sphere-compute-raster.md) (the raster as built: a workgroup of 64 a chunk, a
  thread a sphere over its pixels, one 64-bit `atomicStoreMax` a covered pixel, spheres above 64 pixels strided over
  the workgroup, instances of one sphere in rows of their own).

## Why

- The yeast cell's raster takes 5.9 ms of its 6.4 ms frame (RTX 4070 Ti, whole), its workgroups two thirds of it (Results); MG_1184
  whole costs about 4.5 ns a chunk, a workgroup per visible instance (`sphere-compute-raster.md` step 8).
- What the raster is bound by, measured: the dependent reads of a chunk (its words, the pose, the sphere, its size)
  in the latency of every chunk (decoding quantized poses there: 5.9 -> 7.7-7.9 ms, `particle-pack.md` step 8); and
  zoomed in at 2560×1440 the pixels (ray tests 2.21 ms, atomics 1.79). Not the pose built per lane: lanes run in
  lockstep and one address is one read (the rows once a chunk: 5.9 -> 6.2 ms). Not the arithmetic per pixel (built,
  1.55 -> 1.49 ms).
- What others did, for primitives of a few pixels:
  - **Schütz et al. 2021** (points by compute, 64-bit `atomicMin` of depth and color): the atomics on one pixel
    reduced in the warp first (`subgroupMin` where all lanes share a pixel, pairs by clustered ops, or partitions):
    an overview of high overdraw 21.45 -> 4.29 ms, close up 5.07 -> 3.91; an early depth read before the atomic
    about 34 % on average, most with overdraw; batches of 128 points in Morton order, the batches shuffled, 2-3.5×
    over the original order (clustered writes contend, scattered ones miss the cache).
  - **Schütz et al. 2022** (2 billion points): a workgroup a batch of about 10k points, each thread several;
    positions per batch at the precision its size on screen needs (often 10 bits an axis), a bandwidth gain of
    brute force over unique points.
  - **Nanite** (Karis 2021): triangles up to about 32 pixels by a compute rasterizer, a thread a triangle, about 3×
    the hardware's; larger by the hardware; a 64-bit visibility buffer; materials resolved by screen tiles.
  - **CuRast** (Schütz et al. 2026): three paths by size: a thread a triangle up to 128 pixels (a persistent kernel,
    256 threads taking the next 256 triangles by a global counter), a warp a triangle up to 4096, a workgroup a
    64×64 tile beyond; no shared memory, no warp intrinsics, no contention scheme; whole bounding boxes traversed
    (inside only cost more); register pressure the measured lever (two kernel variants, with and without instancing:
    the instancing code cost occupancy).

## Design

1. **Measured first** (the yeast cell, MG_1184, PetWorld; whole, zoomed to 0.3 and 0.03, 1280×720 and 2560×1440): the
   raster's stages as step 8 of the raster's plan probed them (the spheres, the ray tests, the atomics), and counts:
   chunks and singles, lanes with a sphere, spheres of one pixel or none, pixels tested, atomics, atomics a covered
   pixel. Each step below is kept or dropped by these cases.
2. **Chunks packed** (done): a workgroup takes 4 consecutive chunks, or 16 where the object's chunks hold 16 spheres or
   fewer on average, their spheres dealt over its threads 64 a round (a prefix of the chunks' counts; thread `t` of
   round `r` the sphere `64 r + t`), ids as before (chunk * 64 + its place in the chunk), the append and the resolve
   unchanged. The append counts the chunks' spheres (a counter word); the arguments kernel picks the kernel and gives
   the other none (two indirect dispatches an object). Two kernels, the count a constant of each.
   - Built first by width (`w` = the next power of two, `64 / w` instances a workgroup): not needed, a prefix over
     whole chunks packs as well and keeps the rows.
   - Measured, the raster pass (ms; chunks a workgroup fixed, then as built):

     | chunks a workgroup | yeast | MG whole | MG 0.3 | MG 0.03 | MG 0.3, 2560×1440 | PetWorld whole | PetWorld 0.1 |
     | ------------------ | ----- | -------- | ------ | ------- | ----------------- | -------------- | ------------ |
     | before             | 5.88  | 0.52     | 1.71   | 1.32    | 5.67              | 0.40-0.53      | 1.45         |
     | 1                  | 6.36  | 0.54     | 1.79   | 1.39    | 5.97              | 0.35           | 1.50         |
     | 2                  | 4.79  | 0.54     | 1.58   | 1.27    | 5.37              | 0.44           | 1.51         |
     | 4                  | 4.05  | 0.54     | 1.50   | 1.29    | 5.37              | 0.33           | 1.45         |
     | 8                  | 3.66  | 0.58     | 1.66   | 1.51    | 6.28              | 0.33           | 2.05         |
     | 16                 | 3.52  | 0.57     | 2.05   | 2.01    | 7.79              | 0.26           | 3.01         |
     | 4 or 16 (as built) | 3.53  | 0.55     | 1.50   | 1.28    | 5.40              | 0.20-0.26      | 1.45         |

     The yeast cell's frame 6.37 -> 4.04 ms. More chunks a workgroup lose where they are full (rounds in series where
     workgroups ran side by side). The count chosen on the GPU per object by the spheres a chunk and read at run time
     (one kernel, its loops bound by a uniform count): 5-10 % slower than the fixed counts (the count pinned with the
     rest kept: as fixed), so two kernels.
3. **Persistent workgroups** (built, measured, dropped): at most `P` workgroups dispatched, each drawing the groups
   `w`, `w + P`, ... (the counts read only, so uniform; a barrier between groups for the large spheres' queue).
   Raster pass, ms (yeast, MG whole, 0.3, 0.03, 0.3 at 2560×1440, PetWorld whole, 0.1): one workgroup a group 3.59,
   0.62, 1.56, 1.45, 5.74, 0.20, 1.63; `P` 512 3.79, 0.62, 1.93, 1.75, 7.19, 0.38, 2.10; 2048 3.64, 0.58, 1.70, 1.50,
   6.23, 0.29, 1.79; 8192 3.53, 0.60, 1.66, 1.48, 5.76, 0.31, 1.81. After packing (73k workgroups for the yeast cell)
   launching is not what is left, and a fixed share of groups unequal in work (pixels) leaves workgroups idle at the
   end. A queue by an atomic counter would balance, but reading ahead (the reason) has nothing to gain then. Not
   pursued.
4. **Pixels balanced over the workgroup** (done): each thread looped over its own sphere's pixels, a warp as long as
   its largest sphere's (spheres up to 64 pixels, 19 tested on average zoomed to 0.3). Now a round's spheres of more
   than one pixel are queued in workgroup memory, their pixels in a row (one atomic add of the count and the pixels
   packed, 7 and 25 bits, so the first pixels rise with the slots), strided over all threads; each thread finds a
   pixel's sphere by walking on from the last (its pixels rise). Spheres of one pixel stay with their thread, those
   above `SphereRasterLargePixels` (4096, not to overflow the row) are drawn one after another as the large ones were.
   Built without counting first.
   - Raster pass, ms (yeast, MG whole, 0.3, 0.03, 0.3 at 2560×1440, PetWorld whole, 0.1): before 3.53, 0.55, 1.50,
     1.28, 5.40, 0.20-0.26, 1.45; balanced 3.54, 0.45, 1.27, 1.14-1.19, 4.66, 0.33-0.40 (frame unchanged, 1.02-1.05;
     the case swings 0.16-0.44), 1.24. Frames: MG zoomed to 0.3 at 2560×1440 7.34 -> 6.67 ms, PetWorld zoomed to 0.1
     2.16 -> 1.88.
   - The threshold for sharing (more than 1, 2, 4 or 8 pixels) changes nothing beyond noise. A binary search of the
     first pixels per pixel instead of the walk: slower than before (MG 0.03 1.55, PetWorld 0.1 1.98).
5. **Fewer registers** (measured, dropped): kernel variants with the branches an object does not take compiled out, as
   CuRast's variant without instancing. Tried with the largest, clipping (`clipped` a constant false: the clip flags,
   the per-pixel clip tests and exits gone): raster pass, ms (yeast, MG whole, 0.3, 0.03, 0.3 at 2560×1440, PetWorld
   whole, 0.1) 3.51, 0.44, 1.28, 1.23-1.27, 4.66-4.69, 0.20-0.31, 1.26-1.29 against 3.54, 0.45, 1.27, 1.14-1.19, 4.66,
   0.33-0.40, 1.24: none beyond noise, so the others (singles, levels) not tried.

Not planned, with why:

- A workgroup tile of depths (built, step 1 of the raster: slower) and bins by screen tile (a sort for primitives of
  a pixel costs more than drawing them).
- Pixels inside the disc only (CuRast: the bounding box faster; the arithmetic per pixel not the cost here).
- Positions quantized per chunk (2022's adaptive precision): bandwidth for unique points; here instanced, the spheres
  in the cache, and decoding measured slower (the poses).
- **Atomics reduced in the subgroup** (Schütz 2021) and the early depth read: measured out (Results). Where lanes share
  pixels (whole: 42 % of the yeast cell's one-pixel lanes share one with an earlier lane of their 32) the atomics cost
  0.2-0.8 ms; where they cost (zoomed, 5.2 of 18 ms at 2560×1440) spheres cover many pixels and no two lanes of a
  chunk share one. So also not a shuffled order of chunks.
- Materials by screen tile for the resolve (Nanite): few eligible objects, each a scissored draw (the opaque pass
  0.04-0.24 ms with the raster).

## Steps

1. The measurements (item 1; done).
2. Chunks packed (item 2; done).
3. Persistent workgroups (item 3; measured, dropped).
4. Pixels balanced (item 4; done).
5. Kernel variants (item 5; measured, dropped).

Each step: `renderer/sphere-raster/**` on WebGPU, the Dawn raster test, the cases of step 1; kept only where the yeast
cell or MG_1184 gains and none loses beyond noise.

## Open questions

- Unmeasured: other GPUs than the RTX 4070 Ti (the iGPU, AMD, Apple), where packing and the shared pixels may weigh
  differently (subgroup sizes 8-64, workgroup memory).

## Results

RTX 4070 Ti, Chrome (D3D12), the raster pass (`--gpu-timing`), ms. Measured by probes (a uniform switching the kernel
to return after a stage, set after the warm-up frames, the stage's values kept alive) and counters (workgroup atomics
summed per frame). The yeast cell in a page of its own: after it, the other cases' raster passes up to 3× slower
(MG_1184 whole 0.53 -> 1.65 ms).

- **Occlusion off** (the culler's `occluded` false: probes write no depth, so with occlusion the pyramid empties and
  the culler draws more): the same work in every mode, each stage the difference to the row above.

  | stage reached                | yeast | MG whole | MG 0.3 | MG 0.03 | MG 0.3, 2560×1440 | PetWorld whole | PetWorld 0.1 |
  | ---------------------------- | ----- | -------- | ------ | ------- | ----------------- | -------------- | ------------ |
  | kernel empty (and appending) | 37.1  | 0.45     | 1.45   | 0.88    | 2.58              | 0.62           | 0.66         |
  | chunk words                  | 38.3  | 0.56     | 1.48   | 0.88    | 2.62              | 0.66           | 0.78         |
  | sphere read                  | 39.0  | 0.58     | 1.59   | 0.92    | 2.82              | 0.85           | 0.77         |
  | pose, size, level            | 43.1  | 0.79     | 1.77   | 1.00    | 3.09              | 0.95           | 0.90         |
  | pixel bounds                 | 48.4  | 0.79     | 1.93   | 1.08    | 3.34              | 1.10           | 0.97         |
  | ray tests                    | 55.0  | 0.96     | 5.14   | 2.87    | 12.8              | 1.30           | 2.77         |
  | atomics (the kernel whole)   | 55.8  | 1.33     | 6.35   | 3.50    | 18.0              | 1.35           | 3.79         |

  The yeast cell: 18.8 M chunk workgroups and 0.26 M of singles for 128.8 M spheres (6.8 a chunk).
- **As drawn** (occlusion on): the yeast cell 6.14 in all, 5.17 to the pose, size and level, 5.46 to the bounds; 1.17
  M chunk workgroups and 8.0k of singles, 8.47 M spheres (7.3 a chunk, 11 % of the lanes), 3.51 M of them on a pixel,
  each on one, 2.75 M atomics. MG_1184 whole 11.1k chunks and 4.2k of singles, 676k spheres, 14 pixels a sphere
  tested; zoomed to 0.3 174k chunks, 6.7 M spheres (38 a chunk), 19 pixels a sphere; PetWorld whole 49k chunks, 742k
  spheres (15 a chunk).
- So: the yeast cell's cost is its workgroups, two thirds before any read and 1.4 % the atomics; zoomed in the pixel
  loop (MG zoomed to 0.3 at 2560×1440: ray tests 9.5, atomics 5.2 of 18.0), no lane sharing a pixel there.
