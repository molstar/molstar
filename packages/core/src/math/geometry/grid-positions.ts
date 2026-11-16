/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../linear-algebra/3d/vec3.js';

/**
 * Positions on a grid, 8 B an element where float32 is 12 B: `origin + step * u`, `u` 21 bits an axis (0 to
 * `GridPositions.Steps`) in two u32s, `x | y << 21` and `y >> 11 | z << 10`. Over a box of extent `e` within
 * `e / 2^22` (the yeast cell's 47,000 Å: 0.011 Å; float32 there 0.002 Å). Readers that go over all elements branch on
 * the kind once and decode where they read (`GridPositions.x`, `y`, `z`), no decoded copy made
 */
export interface GridPositions {
  readonly kind: 'grid';
  /** two u32s an element */
  readonly data: Uint32Array;
  readonly origin: Vec3;
  readonly step: number;
}

/** Positions as float32 (x, y, z an element) or on a grid */
export type Positions = Float32Array | GridPositions;

export namespace GridPositions {
  /** The largest value of an axis */
  export const Steps = (1 << 21) - 1;

  export function is(p: Positions | undefined): p is GridPositions {
    return !!p && !(p instanceof Float32Array) && p.kind === 'grid';
  }

  /** Element `i`'s x on the grid of `data` (the loops over all inline these) */
  export function x(data: Uint32Array, i: number) {
    return data[i * 2] & 0x1fffff;
  }
  export function y(data: Uint32Array, i: number) {
    return (data[i * 2] >>> 21) | ((data[i * 2 + 1] & 0x3ff) << 11);
  }
  export function z(data: Uint32Array, i: number) {
    return (data[i * 2 + 1] >>> 10) & 0x1fffff;
  }

  /** Element `i`'s position of `p` (either kind) into `out` (x, y, z) */
  export function get<T extends Vec3 | number[]>(p: Positions, i: number, out: T): T {
    if (p instanceof Float32Array) {
      out[0] = p[i * 3];
      out[1] = p[i * 3 + 1];
      out[2] = p[i * 3 + 2];
      return out;
    }
    const { data, origin, step } = p;
    out[0] = origin[0] + step * x(data, i);
    out[1] = origin[1] + step * y(data, i);
    out[2] = origin[2] + step * z(data, i);
    return out;
  }

  /** The number of elements of `p` */
  export function count(p: Positions) {
    return p instanceof Float32Array ? p.length / 3 : p.data.length / 2;
  }

  /**
   * The grid of a box: `origin` its min, `step` its larger extent over `Steps` (a little over: the float32 rounding of
   * the readers within it); 1 for a point
   */
  export function gridOf(min: ArrayLike<number>, max: ArrayLike<number>) {
    const extent = Math.max(max[0] - min[0], max[1] - min[1], max[2] - min[2]);
    const step = extent > 0 ? (extent * (1 + 1e-6)) / Steps : 1;
    return { origin: Vec3.create(min[0], min[1], min[2]), step };
  }

  /** Writes element `i` of `data` at `x, y, z` on the grid (`origin`, `1 / step`), clamped to it */
  export function set(data: Uint32Array, i: number, origin: Vec3, inv: number, x: number, y: number, z: number) {
    const ux = Math.max(0, Math.min(Steps, Math.round((x - origin[0]) * inv)));
    const uy = Math.max(0, Math.min(Steps, Math.round((y - origin[1]) * inv)));
    const uz = Math.max(0, Math.min(Steps, Math.round((z - origin[2]) * inv)));
    data[i * 2] = (ux | (uy << 21)) >>> 0;
    data[i * 2 + 1] = ((uy >>> 11) | (uz << 10)) >>> 0;
  }

  /** `positions` (float32, x, y, z an element) on a grid of their box (`gridOf`) */
  export function fromFloat(positions: Float32Array): GridPositions {
    const n = positions.length / 3;
    if (n === 0) return { kind: 'grid', data: new Uint32Array(0), origin: Vec3(), step: 1 };
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n * 3; i += 3) {
      for (let k = 0; k < 3; ++k) {
        const v = positions[i + k];
        if (v < min[k]) min[k] = v;
        if (v > max[k]) max[k] = v;
      }
    }
    const { origin, step } = gridOf(min, max);
    const data = new Uint32Array(n * 2);
    const inv = 1 / step;
    for (let i = 0; i < n; ++i) set(data, i, origin, inv, positions[i * 3], positions[i * 3 + 1], positions[i * 3 + 2]);
    return { kind: 'grid', data, origin, step };
  }

  /** `p` as float32 (x, y, z an element): `p` itself where float32 */
  export function toFloat(p: Positions): Float32Array {
    if (p instanceof Float32Array) return p;
    const n = count(p);
    const out = new Float32Array(n * 3);
    const tmp = [0, 0, 0];
    for (let i = 0; i < n; ++i) {
      get(p, i, tmp);
      out[i * 3] = tmp[0];
      out[i * 3 + 1] = tmp[1];
      out[i * 3 + 2] = tmp[2];
    }
    return out;
  }
}
