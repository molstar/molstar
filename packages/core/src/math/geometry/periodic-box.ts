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
   * The displacement `d` moved by a lattice vector to the image nearest to zero in fractional coordinates: exact for
   * a rectangular box, and for any box GROMACS accepts while |d| is under half the shortest box vector.
   */
  export function minimumImage(out: Vec3, box: PeriodicBox, d: Vec3): Vec3 {
    const f = Vec3.transformMat3(out, d, box.inverse);
    f[0] -= Math.round(f[0]);
    f[1] -= Math.round(f[1]);
    f[2] -= Math.round(f[2]);
    return Vec3.transformMat3(out, f, box.basis);
  }
}
