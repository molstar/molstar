/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { chainsOf, readStructure, structureOf } from '../index.js';
import { bytesOf } from './fixtures.js';

// The chain index the complex, the structure providers and the link read (chains.ts).

describe('chainsOf', () => {
  it('lists water chains with no elements; every other chain indexes its elements', async () => {
    const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
    const { chains, chainAt, elementAt } = chainsOf(structure);
    const waters = chains.filter((c) => c.kind === 'water');
    const others = chains.filter((c) => c.kind !== 'water');
    expect(waters.length).toBeGreaterThan(0);
    expect(others.some((c) => c.kind === 'polymer')).toBe(true);
    expect(others.some((c) => c.kind === 'ligand')).toBe(true);
    for (const c of waters) expect(c.count).toBe(0);
    for (const c of others) expect(c.count).toBeGreaterThan(0);
    // A water atom is no chain's element; every element of the others is mapped back.
    for (const c of waters) {
      const at = chainAt.get(c.unit.id)!;
      const k = chains.indexOf(c);
      expect(at.includes(k)).toBe(false);
    }
    for (const c of others) {
      const k = chains.indexOf(c);
      const cAt = chainAt.get(c.unit.id)!;
      const eAt = elementAt.get(c.unit.id)!;
      for (let e = 0; e < c.count; e++) {
        for (let i = c.start[e]; i < c.end[e]; i++) {
          expect(cAt[i]).toBe(k);
          expect(eAt[i]).toBe(e);
        }
      }
    }
  });
});
