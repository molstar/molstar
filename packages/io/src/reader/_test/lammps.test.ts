/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parseLammpsData } from '../lammps/data/parser.js';
import { parseLammpsTrajectory } from '../lammps/traj/parser.js';
import { plainResult } from './plain-result.js';

const dataString = `LAMMPS data file, a water 中

3 atoms
2 atom types
2 bonds

-5.0 5.0 xlo xhi
-5.0 5.0 ylo yhi
-5.0 5.0 zlo zhi
1.0 0.0 0.0 xy xz yz

Masses

# type mass
1 15.9994 # Ö
2 1.008

Atoms # full

1 1 1 -0.834 0.0 0.0 0.0
2 1 2 0.417 0.9572 0.0 0.0
3 1 2 0.417 -0.24 0.927 0.0

Velocities

1 0.0 0.0 0.0
2 0.0 0.0 0.0
3 0.0 0.0 0.0

Bonds

1 1 1 2
2 1 1 3
`;

const dump = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../data/examples/trajectory/periodic/long-chain.lammpstrj'),
  'utf8',
);

describe('lammps data reader', () => {
  it('parses bytes as it parses the text, non-ASCII symbols included', async () => {
    const fromString = await parseLammpsData(dataString).run();
    const fromBytes = await parseLammpsData(utf8Encode(dataString)).run();
    if (fromBytes.isError) throw new Error(fromBytes.message);
    const { atoms, bonds, masses, box } = fromBytes.result;
    expect(atoms.count).toBe(3);
    expect(atoms.x.value(1)).toBeCloseTo(0.9572, 6);
    expect(atoms.charge.value(0)).toBeCloseTo(-0.834, 6);
    expect(bonds.atomIdB.toArray()).toEqual([2, 3]);
    expect(masses?.symbol.toArray()).toEqual(['Ö', '']);
    expect(box?.tilt).toEqual([1, 0, 0]);
    expect(plainResult(fromBytes)).toEqual(plainResult(fromString));
  });
});

describe('lammps trajectory reader', () => {
  it('parses bytes as it parses the text', async () => {
    const fromString = await parseLammpsTrajectory(dump).run();
    const fromBytes = await parseLammpsTrajectory(utf8Encode(dump)).run();
    if (fromBytes.isError) throw new Error(fromBytes.message);
    expect(fromBytes.result.frames.length).toBeGreaterThan(1);
    expect(fromBytes.result.frames[0].image).toBeDefined();
    expect(plainResult(fromBytes)).toEqual(plainResult(fromString));
  });
});
