/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { LammpsBox } from '@molstar/io/reader/lammps/schema';

/**
 * A LAMMPS box as a periodic box (its corner as the origin) and its cell, in Å: a = (lx, 0, 0), b = (xy, ly, 0),
 * c = (xz, yz, lz), already in the cell's standard frame. Only a box periodic in every dimension: a fixed or
 * shrink-wrapped dimension is not a lattice, and its molecules must not be imaged across it.
 */
export function lammpsPeriodicBox(box: LammpsBox, scale: number): { box: PeriodicBox; cell: Cell } | undefined {
  if (!box.periodicity.every((p) => p.startsWith('p'))) return;
  const [lx, ly, lz] = box.length;
  const [xy, xz, yz] = box.tilt ?? [0, 0, 0];
  const a = Vec3.create(lx * scale, 0, 0);
  const b = Vec3.create(xy * scale, ly * scale, 0);
  const c = Vec3.create(xz * scale, yz * scale, lz * scale);
  const origin = Vec3.create(box.lower[0] * scale, box.lower[1] * scale, box.lower[2] * scale);
  const periodic = PeriodicBox.create(a, b, c, origin);
  if (!PeriodicBox.isValid(periodic)) return;
  return { box: periodic, cell: Cell.fromBasis(a, b, c) };
}
