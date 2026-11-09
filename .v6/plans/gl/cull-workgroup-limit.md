# Plan: the merged culler past the workgroup limit

## Status

- **Built**: the culler's three dispatches through `getDispatchSize` (`prefix-sum.ts`, which the prefix sum, the grid
  kernels, the density and color smoothing already split with: no new helper), the kernels' index linear over the rows,
  the instance kernel's workgroup segment by the linear workgroup; the sphere raster's chunk dispatch through it too.
  `MergedCullerLimits.maxWorkgroupsPerDimension` lowers the limit for tests, `MergedCullerCounts.rows` counts the split
  dispatches. Test: `webgpu-gpu-culling` (a workgroup a row: the records and the image as unsplit and as the CPU's;
  the old indexing put the second row's instances in the wrong records).
- **Next limit**: the region buffers (`attributes`, `(regions + 1) × instances × bytes`: every region sized for every
  instance) pass WebGPU's 2 GiB a buffer at about 6 M instances with the raster's and the ghost's regions (the isg
  model 3×3×2: 2.98 GB). Open, see the README.

## Question

The merged culler (`merged-culler.ts`) dispatches its three kernels (cells, instances, records; `cell-cull.wgsl.ts`)
over one dimension: `ceil(n / 64)` workgroups in x. WebGPU guarantees 65,535 workgroups per dimension
(`maxComputeWorkgroupsPerDimension`), so a kernel covers at most about 4.2 M items. A mesoscale model past it fails at
its first GPU-culled frame (`kernel 'cell-cull': 83305 × 1 × 1 workgroups exceed 65535 per dimension`): the
isg_mature model repeated 3×3×2 (8.28 M instances, 5.3 M cells). The instance kernel fails next (8.28 M / 64 = 129k
workgroups).

## Design

Split the dispatches over two dimensions, as the sphere raster's chunk kernel does (`sphere-raster.ts`).

- **A shared split.** A helper next to the compute utilities: `n` workgroups as `[min(n, max), ceil(n / max)]`, `max`
  the guaranteed 65,535 (no device query: every device has it). The sphere raster's chunk dispatch takes it instead of
  its literal. A test can lower `max` to exercise the split at small sizes (a module setting, restored by the test).
- **The kernels' index.** The three entry points read `num_workgroups` and take `id.x + id.y * groups.x *
  WorkgroupSizeX` as their item, as the raster's kernels do; their bounds tests (`uCellCount`, `uInstanceCount`,
  `uRecordCount`) already end the last row's surplus. The instance kernel's segment of its workgroup's first instance
  (`instanceCells[wid.x * WorkgroupSizeX]`) takes the linear workgroup, `wid.x + wid.y * groups.x`: with `wid.x` every
  row past the first reads another workgroup's segment and puts instances in the wrong regions. Barriers stay where
  they are (every invocation reaches them).
- **The culler** dispatches through the split at both phases.
- **Not affected**: the counts are integer uniforms (exact far past these sizes); WebGL runs no culling kernels; the
  other kernels' dispatches are two-dimensional already or bounded by a texture's or grid's size.

## Tests

- The split: counts below, at and past `max`.
- Native WebGPU (`webgpu-gpu-culling`): the culled lists with `max` lowered to a few workgroups equal those without
  (levels, the raster's regions and the ghost's, which depend on the workgroup's segment). The real limit needs over
  4.2 M instances: checked by hand on the grid models, not in a lane.
