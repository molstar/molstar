/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat3, Vec3 } from '../linear-algebra.js';
import type { Cell } from './spacegroup/cell.js';

export { PeriodicBox };

/**
 * A simulation's periodic box: the lattice its coordinates repeat on. Unlike a crystal's `Cell`, it keeps the box's
 * own basis (not rebuilt in the standard frame) and its origin.
 */
interface PeriodicBox {
  /** The box vectors a, b, c as columns: cartesian = basis * fractional */
  readonly basis: Mat3;
  /** fractional = inverse * cartesian */
  readonly inverse: Mat3;
  /** The box's corner; the lattice does not depend on it, only wrapping into the box does */
  readonly origin: Vec3;
}

namespace PeriodicBox {
  export function create(a: Vec3, b: Vec3, c: Vec3, origin: Vec3 = Vec3()): PeriodicBox {
    const basis = Mat3.fromColumns(Mat3(), a, b, c);
    return { basis, inverse: Mat3.invert(Mat3(), basis), origin };
  }

  /** The box of a cell, in the cell's standard frame (a along x, b in the xy plane) and at the origin. */
  export function fromCell(cell: Cell): PeriodicBox {
    const basis = Mat3.fromMat4(Mat3(), cell.fromFractional);
    return { basis, inverse: Mat3.invert(Mat3(), basis), origin: Vec3() };
  }

  /** A box whose vectors are all longer than zero and not in one plane. */
  export function isValid(box: PeriodicBox): boolean {
    const d = Mat3.determinant(box.basis);
    return Number.isFinite(d) && Math.abs(d) > 1e-6;
  }

  /**
   * The displacement `d` moved by a lattice vector to the shortest image. Rounding the fractional coordinates gives
   * the image in the box-shaped cell around zero, which is the shortest only in a rectangular box. In a skewed box the
   * shorter images are searched exactly: an image no longer than the rounded one, |m|, has its fractional coordinate i
   * within |m| · |row i of the inverse| of zero, which bounds the lattice shifts to try (the 27 around the rounded
   * image for GROMACS's dodecahedron and octahedron; more for thin, tilted boxes) (perf-audit-3 02/B7).
   */
  export function minimumImage(out: Vec3, box: PeriodicBox, d: Vec3): Vec3 {
    const f = Vec3.transformMat3(out, d, box.inverse);
    const fx = f[0] - Math.round(f[0]),
      fy = f[1] - Math.round(f[1]),
      fz = f[2] - Math.round(f[2]);
    const m = box.basis;
    const x = fx * m[0] + fy * m[3] + fz * m[6],
      y = fx * m[1] + fy * m[4] + fz * m[7],
      z = fx * m[2] + fy * m[5] + fz * m[8];
    if (m[1] === 0 && m[2] === 0 && m[3] === 0 && m[5] === 0 && m[6] === 0 && m[7] === 0) {
      return Vec3.set(out, x, y, z);
    }

    let best = x * x + y * y + z * z;
    const r = Math.sqrt(best);
    const v = box.inverse;
    const hx = r * Math.hypot(v[0], v[3], v[6]),
      hy = r * Math.hypot(v[1], v[4], v[7]),
      hz = r * Math.hypot(v[2], v[5], v[8]);
    let bx = x,
      by = y,
      bz = z;
    for (let i = Math.ceil(-fx - hx), il = Math.floor(-fx + hx); i <= il; i++) {
      for (let j = Math.ceil(-fy - hy), jl = Math.floor(-fy + hy); j <= jl; j++) {
        for (let k = Math.ceil(-fz - hz), kl = Math.floor(-fz + hz); k <= kl; k++) {
          const sx = x + i * m[0] + j * m[3] + k * m[6],
            sy = y + i * m[1] + j * m[4] + k * m[7],
            sz = z + i * m[2] + j * m[5] + k * m[8];
          const l = sx * sx + sy * sy + sz * sz;
          if (l < best) {
            best = l;
            bx = sx;
            by = sy;
            bz = sz;
          }
        }
      }
    }
    return Vec3.set(out, bx, by, bz);
  }
}
