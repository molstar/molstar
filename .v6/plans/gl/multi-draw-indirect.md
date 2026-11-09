# Plan: multi-draw indirect

## Status

- **Current**: done: steps 1 to 3 done, step 4 dropped after measuring, step 5 nothing to do. Multi-draw indirect (MDI) stays where it is (GPU-culled records, behind the flag); what it
  would add needs a layout whose gain comes without it, so the steps below build that and measure MDI on top only at
  the end.
  - **Step 2 done, replaying unculled instead of culling on the GPU**: indirect draws in a bundle cost the GPU about
    0.2-0.3 µs each more than direct ones, drawn or empty (10k draws, `draw-overhead/*`, Chrome 156: dynamic offsets
    1.12 -> 4.24 ms indirect, 4.08 with half empty; by `firstInstance` 0.96 -> 2.77, 2.84 all empty), as the
    draw-table plan measured before ([draw-table.md](draw-table.md) step 7): per-object GPU culling adds GPU time for
    every object every frame, and only MDI's compaction (not in bundles) would take it back. Instead a pass replayed
    from a draw list draws its rows unculled while the visible set changed within 8 frames and at least 15% of the
    present rows are visible (`getRowVisibility`: `present` is `visible` without the frustum and occlusion tests; state,
    whole clipping and level ranges kept): its recording holds as objects enter and leave the view, the GPU clipping
    what is off-screen (under 1 µs a draw); at rest culled again, recorded once. WebGPU, rotating:
    `boxes-1000-zoom0.3-rotating` 1.70 -> 0.67 ms (submit 1.55 -> 0.25, nothing recorded again), `-zoom0.1-rotating`
    1.70 -> 0.70 ms; at rest unchanged. `webgpu-draw-lists`: a view moved back and forth so that a surface leaves and
    enters it, frames as drawn directly, at most 2 chunks not replayed over the moves (12 without). Not covered:
    objects culled per instance (instance grids, `instances-1000-spacefill-unmerged-rotating` 1.62 ms,
    `assembly-60-spacefill-rotating`): their draws change with the view.
  - **Found after**: spheres without levels of detail (their default, an empty `lodLevels`) counted as having them
    (`getCullFacts` tested the list, not its length): view dependent and culled every frame, a Viewer-like scene's
    passes not replayed while rotating. `render/viewer-7qpd-rotating` (cartoon, ball-and-stick, two surfaces): 58
    chunks drawn directly -> none, every pass replayed unplanned, submit 0.33 -> 0.19 ms (the frame GPU bound, 0.86 ->
    0.83 ms).
  - **Step 1 done**: a chunk of objects counted as changing with the view where any of its objects has a solid interior
    (`isViewDependent`), for the caps of the near cut; but a pass draws the caps with the objects only where they are
    inline (blended transparency), the others in their own chunk, which compares all it read. So a solid interior
    counts only for inline caps. `units-1356-gaussian-rotating` (solid interiors on by default for Gaussian surfaces):
    WebGPU 6.69 -> 0.56 ms (submit 6.52 -> 0.45, every pass replayed unplanned); WebGL2 unchanged (no recordings).
    `webgpu-draw-lists`: the camera moved sideways and to a near cut through two surfaces, the frames as drawn
    directly, no chunk drawn directly (4 without the change).

## Where things stand

- **Availability** (Chrome 156, probed): `chromium-experimental-multi-draw-indirect` only with `--enable-unsafe-webgpu`,
  on render pass encoders, not on bundle encoders; the spec proposal (gpuweb PR 2315) unmerged and not moving. No user
  has it. Shipped by default: `indirect-first-instance` and immediates (`setImmediates` on passes and bundles,
  `maxImmediateSize` 64, WGSL `immediate_address_space`).
- **Used**: objects culled on the GPU draw their indirect records by one call per run of records
  (`render-item.ts`, capability `multiDrawIndirect`), on render passes only, so their chunks are drawn directly.
  Measured with the empty records left out: it adds about nothing (MG_1184 whole rotating 1.15 -> 1.28 ms, zoomed
  2.58 -> 2.79; [sphere-compute-raster.md](sphere-compute-raster.md) step 8).
- **Where draws cost the CPU now** (WebGPU, rotating, the benchmark's Chrome with the flag; frame, submit, draw calls):

  | case                                                 | frame   | submit  | draws |
  | ---------------------------------------------------- | ------- | ------- | ----- |
  | `render-bundles/units-1356-gaussian-rotating`        | 6.69 ms | 6.52 ms | 1363  |
  | (the same at rest, `-drawlists`: replayed)           | 0.55 ms | 0.20 ms | 7     |
  | (rotating after step 1)                              | 0.56 ms | 0.45 ms | 7     |
  | `render-bundles/boxes-1000-zoom0.3-rotating`         | 1.70 ms | 1.55 ms | 310   |
  | (the same at rest, culled)                           | 1.08 ms | 0.21 ms | 7     |
  | `merging/instances-1000-spacefill-unmerged-rotating` | 1.55 ms | 1.43 ms | 107   |
  | `merging/instances-1000-spacefill-merged-rotating`   | 1.25 ms | 1.14 ms | 107   |
  | `merging/instances-1000-spacefill-merged-gpu-rotating` | 0.78 ms | 0.53 ms | 31  |
  | `render/assembly-60-spacefill-rotating`              | 0.78 ms | 0.52 ms | 33    |
  | `render-bundles/boxes-1000-drawlists-rotating`       | 0.65 ms | 0.21 ms | 7     |

  Many objects while the camera moves cost; at rest the recorded passes replay. The 1356 units' opaque pass was drawn
  directly every rotating frame (21 chunks a frame, all of them): step 1. Zoomed in, objects entering and leaving the
  view record the pass again (the boxes: 24 of 60 frames). Instances culled on the CPU cost about twice those culled
  on the GPU.

## What MDI could do

1. **A level's culled ranges by one call**: the CPU culling's lists (merged objects, plain objects with a grid) written
   as indirect records and drawn by one call per level. Saves the encoder calls of the ranges (107 draws, about 0.6 ms
   above GPU culling's 31), for flagged users only; GPU culling, which every user has, already gets there.
2. **Many objects by one call**: the objects of a pass sharing a pipeline drawn from one bind group, their per-draw
   data (object block, level uniforms, layer offsets) in a table read by the instance index (`firstInstance`, with
   `indirect-first-instance`), geometry in the shared buffers (done, [shared-draw-resources.md](shared-draw-resources.md)).
   With MDI one call per pipeline; without it one `drawIndexedIndirect` per object, no bind group in between, and
   recordable.
3. **Fewer GPU commands**: measured as nothing on NVIDIA where it is used.

So MDI adds call counts only on top of a layout (2) and a culling (below) that give their gain without it, and it
cannot be recorded.

## Steps

1. Done: **why the 1356 units' pass is not replayed while rotating**: solid interiors counted as changing with the
   view in passes that draw their caps apart.
2. Done: **objects entering and leaving the view replayed**: unculled while the view moves, not culled per object on
   the GPU (indirect draws cost every object every frame).
3. Done: **GPU culling's `auto` threshold** 20k -> 50k instances: CPU and GPU culling even at 50k for spacefill
   rotating, the CPU ahead below by 0.1-0.3 ms and for cartoon up to 100k (numbers in
   [merged-render-objects.md](merged-render-objects.md), cases `merging/instances-{100,300,500,1000}-{spacefill,cartoon}-*`).
   One machine (RTX 4070 Ti, a fast CPU): a slower CPU moves the crossover down, a slower GPU up, where the GPU bounds
   the frame either way (measured on an integrated GPU: nothing lost).
4. Dropped after measuring: **immediates for per-draw data**. 10k draws (`draw-overhead/10k-*`, Chrome 156),
   recorded in a bundle / encoded per frame: a dynamic offset per draw 1.10 / 3.91 ms (submit 1.56), the object's 16
   bytes as immediates 1.12 / 2.19 (0.92), a 4-byte index as immediates 1.10 / 2.40 (0.93), objects by `firstInstance`
   from a storage array 1.00 / 1.52 (0.24). No gain in bundles, where Mol*'s passes replay; drawn directly
   `firstInstance` beats them and is core WebGPU, while immediates reach about 76% of browsers (none of Firefox's and
   Safari's). Per-draw data by `firstInstance` is the draw-submission plan's layer B, if directly drawn passes matter.
5. **MDI**: with step 2 there are no per-object indirect records to draw together; it stays for the GPU-culled
   records, re-measured if the proposal ships and allows bundle encoders.

## Open questions

- Whether the MDI proposal moves (gpuweb PR 2315) and allows bundle encoders; until then nothing beyond step 5.
- Unculled while moving draws the off-screen geometry too: cheap for the boxes (under 1 µs a draw), not measured for
  large off-screen meshes (an unmerged structure's surfaces zoomed into); the 15% share bounds it to about 7 times the
  visible draws.
- Objects culled per instance on the CPU (unmerged assemblies, instanced spacefill) still change their draws with the
  view; merging with GPU culling covers them where it applies.
