/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import { calcUnitBackboneHbonds } from '../secondary-structure/dssp/backbone-hbonds.js';
import { getUnitProteinInfo } from '../secondary-structure/dssp/protein-info.js';
import { calcUnitProteinTraceLookup3D } from '../secondary-structure/dssp/trace-lookup.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../data/examples');

describe('DSSP backbone H-bonds', () => {
  it('store one energy per edge slot (perf-audit-2 06/B8)', async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const unit = Structure.ofModel(model).units.find(Unit.isAtomic)!;
    const info = getUnitProteinInfo(unit);
    const hbonds = calcUnitBackboneHbonds(unit, info, calcUnitProteinTraceLookup3D(unit, info.residueIndices));
    const { energies } = hbonds.edgeProps;
    expect(hbonds.edgeCount).toBeGreaterThan(50);
    // indexed like `a` and `b` (the slots), as every edge property of the graph
    expect(energies).toBeInstanceOf(Float32Array);
    expect(energies.length).toBe(hbonds.a.length);
    for (let i = 0; i < energies.length; i++) {
      expect(energies[i]).toBeLessThanOrEqual(-0.5);
      expect(energies[i]).toBeGreaterThan(-10);
    }
  });
});
