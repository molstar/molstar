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
import { MoleculeType } from '../../../../model/types.js';
import { computeInterUnitBonds } from '../inter-compute.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../../data/examples');

// Water and ion units pair with another unit only through a common struct_conn record (1tup: Zn, metalc); a unit's
// partners are read once, not per pair and computation (perf-audit-2 05/P17).
describe('inter-unit bonds of ion units', () => {
  it('read the struct_conn rows of a unit once; the bonds are the same', async () => {
    const parsed = await CIF.parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, '1tup.bcif')))).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
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
