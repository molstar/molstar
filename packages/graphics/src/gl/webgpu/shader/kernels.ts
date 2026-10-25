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

/** WGSL compute kernels by the name of their `KernelCode`, like `programs.ts` for render programs */
const Kernels = new Map<string, string>([
  ['scan-block', scan_block_wgsl],
  ['scan-add', scan_add_wgsl],
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
]);

/** Registers the WGSL source of a kernel; the bindings are declared from its schema, see `getWgslKernelDeclarations` */
export function registerWgslKernel(name: string, source: string) {
  if (Kernels.has(name) && Kernels.get(name) !== source) throw new Error(`kernel '${name}' is already registered`);
  Kernels.set(name, source);
}

export function getWgslKernel(name: string): string | undefined {
  return Kernels.get(name);
}

/** The names of the registered kernels (e.g. for a check that every one compiles) */
export function getWgslKernelNames(): string[] {
  return [...Kernels.keys()];
}
