/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Column } from '@molstar/core/data/db';

export interface GroHeader {
  title: string;
  timeInPs: number;
  /** number of decimal places */
  precision: { position: number; velocity: number };
  hasVelocities: boolean;
  /** The box's diagonal, nm: v1(x) v2(y) v3(z) */
  box: [number, number, number];
  /** The box vectors, nm: v1, v2, v3 as x y z each (a triclinic box's off-diagonal values; zero when not given) */
  boxVectors: [number, number, number, number, number, number, number, number, number];
}

export interface GroAtoms {
  count: number;
  residueNumber: Column<number>;
  residueName: Column<string>;
  atomName: Column<string>;
  atomNumber: Column<number>;
  x: Column<number>;
  y: Column<number>;
  z: Column<number>;
  vx: Column<number>;
  vy: Column<number>;
  vz: Column<number>;
}

export interface GroStructure {
  header: Readonly<GroHeader>;
  atoms: Readonly<GroAtoms>;
}

export interface GroFile {
  structures: GroStructure[];
}
