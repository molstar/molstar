/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { GraphicsContext } from './context.js';

/**
 * What an object's content hides, as optical densities (1/Å) on a grid over its box (`dimensions` cells, x fastest):
 * the depth pyramid primed from it (`HiZPass.setOccluder`), a ray's depth where its optical depth reaches the
 * occlusion's threshold, so the GPU's occlusion holds through the gaps between molecules (the pyramid's farthest depth
 * otherwise their background)
 */
export interface OccluderGrid {
  readonly min: Vec3;
  readonly max: Vec3;
  readonly dimensions: readonly [number, number, number];
  readonly density: Float32Array;
  /** its state of `density` where the array is written in place (a new version each time) */
  readonly version?: number;
  /**
   * what changed from version `from` (the same `density`): the boxes of cells written, `[x0, y0, z0, x1, y1, z1)` six a
   * box: its upload those cells' alone where the last uploaded is that version
   */
  readonly update?: { readonly from: number; readonly boxes: Int32Array };
}

/** per object's state (its render object's and renderable's, the same object: a merged object's values copied) its grid */
const grids = new WeakMap<object, OccluderGrid>();

/** `grid` the occluder of the object of `state` (`GraphicsRenderObject.state`; none: none) */
export function setOccluderGrid(state: object, grid: OccluderGrid | undefined) {
  if (grid) grids.set(state, grid);
  else grids.delete(state);
}

export function getOccluderGrid(state: object): OccluderGrid | undefined {
  return grids.get(state);
}

/** whether `graphics` primes its occlusion from a grid (`HiZPass.setOccluder`: WebGPU's): one worth computing */
export function usesOccluderGrid(graphics: GraphicsContext | undefined) {
  return graphics?.capabilities.backend === 'webgpu';
}
