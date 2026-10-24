/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parseCube } from '../cube/parser.js';

const cube = `Molécule test, Ångström
 Electron density
    2   -1.000000   -2.000000    0.500000
    2    0.500000    0.000000    0.000000
    2    0.000000    0.500000    0.000000
    3    0.000000    0.000000    0.500000
    8    8.000000    0.000000    0.000000    0.221665
    1    1.000000    0.000000    1.430901   -0.886658
 1.00000E-01 -2.0E+00 3.5
  4.0 5.0 6.0
 7.0 8.0
 9.0 10.0 11.0  1.2e-3
`;

const orbitals = `orbitals
 MO
   -1    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    1    0.000000    1.000000    0.000000
    1    0.000000    0.000000    1.000000
    1    1.000000    0.000000    0.000000    0.000000
    2    5    6
 0.1 0.2\r
 0.3 0.4\r
`;

describe('cube reader', () => {
  it('parses the header, atoms and values', async () => {
    const parsed = await parseCube(cube, 'test').run();
    if (parsed.isError) throw new Error(parsed.message);
    const { header, atoms, values } = parsed.result;
    expect(header.comment1).toBe('Molécule test, Ångström');
    expect(header.comment2).toBe('Electron density');
    expect(header.atomCount).toBe(2);
    expect(Array.from(header.dim)).toEqual([2, 2, 3]);
    expect(header.orbitals).toBe(false);
    expect(header.dataSetIds).toEqual([0]);
    expect(atoms.count).toBe(2);
    expect(atoms.number.value(0)).toBe(8);
    expect(atoms.nuclearCharge.value(1)).toBe(1);
    expect(Array.from(values)).toEqual([0.1, -2, 3.5, 4, 5, 6, 7, 8, 9, 10, 11, 0.0012]);
  });

  it('parses orbitals with their data set ids', async () => {
    const parsed = await parseCube(orbitals, 'test').run();
    if (parsed.isError) throw new Error(parsed.message);
    const { header, values } = parsed.result;
    expect(header.orbitals).toBe(true);
    expect(header.dataSetIds).toEqual([5, 6]);
    expect(Array.from(values)).toEqual([0.1, 0.2, 0.3, 0.4]);
  });

  it('parses the same from a string and from its UTF-8 bytes', async () => {
    for (const text of [cube, orbitals]) {
      const fromString = await parseCube(text, 'test').run();
      const fromBytes = await parseCube(utf8Encode(text), 'test').run();
      if (fromString.isError || fromBytes.isError) throw new Error('parse failed');
      expect(fromBytes.result.header).toEqual(fromString.result.header);
      expect(fromBytes.result.values).toEqual(fromString.result.values);
      expect(fromBytes.result.atoms.x.toArray()).toEqual(fromString.result.atoms.x.toArray());
    }
  });
});
