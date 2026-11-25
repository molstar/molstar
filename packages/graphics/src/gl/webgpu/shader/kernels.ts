/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { color_accumulate_wgsl, color_normalize_wgsl } from './color-smoothing.wgsl.js';
import { gd_accumulate_wgsl, gd_group_wgsl, gd_write_wgsl } from './gaussian-density.wgsl.js';
import {
  mc_chunk_allocate_wgsl,
  mc_chunk_count_wgsl,
  mc_chunk_emit_wgsl,
  mc_chunk_free_wgsl,
  mc_chunk_push_wgsl,
} from './marching-cubes-chunks.wgsl.js';
import { scan_add_wgsl, scan_block_wgsl } from './prefix-sum.wgsl.js';
import { cell_cull_wgsl, cull_args_wgsl, instance_copy_wgsl } from './cell-cull.wgsl.js';
import { scene_cell_bin_wgsl, scene_cell_spheres_wgsl } from './scene-cells.wgsl.js';
import { culled_poses_wgsl } from './culled-poses.wgsl.js';
import {
  paged_cull_args_wgsl,
  paged_cull_compose_wgsl,
  paged_cull_select_args_wgsl,
  paged_cull_groups_wgsl,
  paged_cull_select_wgsl,
  paged_cull_wants_wgsl,
  paged_cull_wgsl,
} from './paged-cull.wgsl.js';
import { paged_pick_wgsl } from './paged-pick.wgsl.js';
import { u32_runs_wgsl } from './u32-runs.wgsl.js';
import { denoise_wgsl } from './denoise.wgsl.js';
import {
  sphere_raster_args_wgsl,
  sphere_raster_blocks_args_wgsl,
  sphere_raster_blocks_wgsl,
  sphere_raster_cull_wgsl,
  sphere_raster_regions_wgsl,
  sphere_raster_wgsl,
} from './sphere-raster.wgsl.js';
import { SphereRasterChunksPerWorkgroup, SphereRasterWideChunksPerWorkgroup } from '../../compute/sphere-raster.js';
import { splat_centers_wgsl, splat_keys_wgsl, splat_scatter_wgsl } from './splat-sort.wgsl.js';
import { registerWgslKernel } from './kernel-registry.js';

/** WGSL compute kernels by the name of their `KernelCode`, like `programs.ts` for render programs */
const BuiltInKernels: [name: string, source: string][] = [
  ['scan-block', scan_block_wgsl],
  ['scan-add', scan_add_wgsl],
  ['splat-keys', splat_keys_wgsl],
  ['splat-scatter', splat_scatter_wgsl],
  ['splat-centers', splat_centers_wgsl],
  ['mc-chunk-count', mc_chunk_count_wgsl],
  ['mc-chunk-allocate', mc_chunk_allocate_wgsl],
  ['mc-chunk-emit', mc_chunk_emit_wgsl],
  ['mc-chunk-free', mc_chunk_free_wgsl],
  ['mc-chunk-push', mc_chunk_push_wgsl],
  ['color-smoothing-accumulate', color_accumulate_wgsl],
  ['color-smoothing-normalize', color_normalize_wgsl],
  ['gd-accumulate', gd_accumulate_wgsl],
  ['gd-group', gd_group_wgsl],
  ['gd-write', gd_write_wgsl],
  ['cell-cull', cell_cull_wgsl],
  ['instance-copy', instance_copy_wgsl],
  ['cull-args', cull_args_wgsl],
  ['scene-cell-bin', scene_cell_bin_wgsl],
  ['scene-cell-spheres', scene_cell_spheres_wgsl],
  ['culled-poses', culled_poses_wgsl],
  ['paged-cull-compose', paged_cull_compose_wgsl],
  ['paged-cull-select', paged_cull_select_wgsl],
  ['paged-cull-groups', paged_cull_groups_wgsl],
  ['paged-cull-select-args', paged_cull_select_args_wgsl],
  ['paged-cull-wants', paged_cull_wants_wgsl],
  ['paged-cull', paged_cull_wgsl],
  ['paged-cull-args', paged_cull_args_wgsl],
  ['paged-pick', paged_pick_wgsl],
  ['u32-runs', u32_runs_wgsl],
  ['illumination-denoise', denoise_wgsl],
  ['sphere-raster-cull', sphere_raster_cull_wgsl],
  ['sphere-raster-regions', sphere_raster_regions_wgsl],
  ['sphere-raster-blocks', sphere_raster_blocks_wgsl],
  ['sphere-raster-blocks-args', sphere_raster_blocks_args_wgsl],
  ['sphere-raster-args', sphere_raster_args_wgsl],
  ['sphere-raster', sphere_raster_wgsl(SphereRasterChunksPerWorkgroup)],
  ['sphere-raster-wide', sphere_raster_wgsl(SphereRasterWideChunksPerWorkgroup)],
];
for (const [name, source] of BuiltInKernels) registerWgslKernel(name, source);

export { getWgslKernel, getWgslKernelNames, registerWgslKernel } from './kernel-registry.js';
