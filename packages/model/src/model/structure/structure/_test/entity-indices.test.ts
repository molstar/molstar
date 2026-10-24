/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../../../../formats/structure/mmcif.js';
import { StructureElement, StructureProperties } from '../../../structure.js';
import { Structure } from '../structure.js';
import { Unit } from '../unit.js';

async function load8zzh() {
  const file = path.resolve(import.meta.dirname, '../../../../../../../data/examples/ihm-coarse-8zzh.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

describe('Structure.entityIndices', () => {
  it('reads the entities of coarse units by their own chains (perf-audit-2 05/B2)', async () => {
    const structure = await load8zzh();
    const coarse = structure.units.filter((u) => !Unit.isAtomic(u));
    expect(coarse.length).toBeGreaterThan(0);

    // every element's entity, read one by one
    const expected = new Set<number>();
    const l = StructureElement.Location.create(structure);
    for (const unit of structure.units) {
      const prop = Unit.isAtomic(unit) ? StructureProperties.entity.key : StructureProperties.coarse.entityKey;
      l.unit = unit;
      for (let i = 0; i < unit.elements.length; i++) {
        l.element = unit.elements[i];
        expected.add(prop(l));
      }
    }
    expect(expected.size).toBeGreaterThan(0);
    expect([...structure.entityIndices]).toEqual([...expected].sort((a, b) => a - b));
  });
});
