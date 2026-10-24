/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { StructConn } from '@molstar/model/formats/structure/property/bonds/struct_conn';
import { Structure, Unit } from '@molstar/model/model/structure';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { MoleculeType } from '../../../../model/types.js';
import { computeInterUnitBonds } from '../inter-compute.js';
import { getElementThreshold, getPairingThreshold, MaxPairingThreshold } from '../common.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../../data/examples');

async function load1tup() {
  const parsed = await CIF.parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, '1tup.bcif')))).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

describe('inter-unit distance bonds', () => {
  it('MaxPairingThreshold is the largest pairing threshold of any element pair (perf-audit-2 05/P16)', () => {
    let max = 0;
    for (let i = -1; i < 130; i++) {
      for (let j = -1; j < 130; j++) {
        max = Math.max(max, getPairingThreshold(i, j, getElementThreshold(i), getElementThreshold(j)));
      }
    }
    expect(MaxPairingThreshold).toBe(max);
    expect(MaxPairingThreshold).toBeLessThan(4);
  });

  it('are searched at the largest pairing threshold, not the 4 Å candidate radius', async () => {
    const model = await load1tup();
    const structure = Structure.ofModel(model);
    const radii: number[] = [];
    const spies = structure.units.map((u) => {
      const find = u.lookup3d.find.bind(u.lookup3d);
      return vi.spyOn(u.lookup3d, 'find').mockImplementation((x, y, z, r) => (radii.push(r), find(x, y, z, r)));
    });
    const bonds = computeInterUnitBonds(structure);
    for (const s of spies) s.mockRestore();
    expect(bonds.edgeCount).toBeGreaterThan(0);
    expect(radii.filter((r) => r === 4).length).toBe(0);
    expect(radii.filter((r) => r === MaxPairingThreshold).length).toBeGreaterThan(0);
  });
});

// Water and ion units pair with another unit only through a common struct_conn record (1tup: Zn, metalc); a unit's
// partners are read once, not per pair and computation (perf-audit-2 05/P17).
describe('inter-unit bonds of ion units', () => {
  it('read the struct_conn rows of a unit once; the bonds are the same', async () => {
    const model = await load1tup();
    const structure = Structure.ofModel(model);
    const isIon = (u: Unit) =>
      Unit.isAtomic(u) &&
      model.atomicHierarchy.derived.residue.moleculeType[u.residueIndex[u.elements[0]]] === MoleculeType.Ion;
    const ions = structure.units.filter(isIon);
    expect(ions.length).toBeGreaterThan(0);

    const structConn = StructConn.Provider.get(model)!;
    const get = vi.spyOn(structConn.byAtomIndex, 'get');
    const a = computeInterUnitBonds(structure);
    const firstReads = get.mock.calls.length;
    expect(firstReads).toBeGreaterThan(0);

    // the second computation reads only the rows of the atoms that get bonds
    const b = computeInterUnitBonds(structure);
    const secondReads = get.mock.calls.length - firstReads;
    get.mockRestore();
    expect(secondReads).toBeLessThan(firstReads);

    expect(b.edgeCount).toBe(a.edgeCount);
    // every Zn is bonded through its metalc records
    for (const u of ions) expect(a.getConnectedUnits(u.id).length).toBeGreaterThan(0);
  });
});

// A child made of whole parent units (a component of chain units) shares the parent's unit-pair lists between its units
// instead of rebuilding every edge (perf-audit-2 05/P20).
describe('inter-unit bonds of a child of whole units', () => {
  it('are the parent bonds between its units, without a rebuild', async () => {
    const parent = Structure.ofModel(await load1tup());
    const bonds = parent.interUnitBonds;
    // the Zn ions and the units bonded to them, and one unit without bonds
    const zn = parent.units.filter((u) => bonds.getConnectedUnits(u.id).length > 0);
    const lone = parent.units.find((u) => bonds.getConnectedUnits(u.id).length === 0)!;
    const units = [...zn.slice(0, zn.length - 1), lone].sort((a, b) => a.id - b.id);
    expect(units.length).toBeGreaterThan(2);
    const child = Structure.create(units, { parent });

    const add = vi.spyOn(InterUnitGraph.Builder.prototype, 'add');
    const childBonds = child.interUnitBonds;
    const adds = add.mock.calls.length;
    add.mockRestore();
    expect(adds).toBe(0);

    const ids = new Set(units.map((u) => u.id));
    const key = (e: { unitA: number; unitB: number; indexA: number; indexB: number }) =>
      `${e.unitA}:${e.indexA}-${e.unitB}:${e.indexB}`;
    const expected = bonds.edges.filter((e) => ids.has(e.unitA) && ids.has(e.unitB)).map(key);
    expect(expected.length).toBeGreaterThan(0);
    expect(expected.length).toBeLessThan(bonds.edgeCount);
    expect(childBonds.edges.map(key).sort()).toEqual(expected.sort());
    expect(childBonds.edgeCount).toBe(expected.length);
    for (const e of childBonds.edges) {
      expect(childBonds.getEdgeIndex(e.indexA, e.unitA, e.indexB, e.unitB)).toBeGreaterThanOrEqual(0);
    }
  });
});
