/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { RuntimeContext } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Script } from '@molstar/model/script/script';
import { StructureSelection } from '../../../query/selection.js';
import { Structure } from '../../structure.js';
import type { UnitIndex } from '../element.js';
import { Loci } from '../loci.js';
import '@molstar/model/script/runtime/query/table'; // (registers the symbols of the queries)

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

describe('Loci.toExpression', () => {
  it('keeps the elements of each model apart (perf-audit-2 05/B11)', async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1grm_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
    const structure = await Structure.ofTrajectory(trajectory, RuntimeContext.Synchronous);
    expect(structure.models.length).toBeGreaterThan(1);

    // the first unit of the first two models: both have operator 1_555
    const [m0, m1] = structure.models;
    const u0 = structure.units.find((u) => u.model === m0)!;
    const u1 = structure.units.find((u) => u.model === m1)!;
    expect(u0.conformation.operator.name).toBe(u1.conformation.operator.name);
    const loci = Loci(structure, [
      { unit: u0, indices: OrderedSet.ofSortedArray([0, 1, 2] as UnitIndex[]) },
      { unit: u1, indices: OrderedSet.ofSortedArray([5, 6, 7] as UnitIndex[]) },
    ]);

    const selection = Script.getStructureSelection(Loci.toExpression(loci), structure);
    const back = StructureSelection.toLociWithSourceUnits(selection);
    const picked = back.elements.map((e) => [
      e.unit.model.modelNum,
      ...OrderedSet.toArray(e.indices).map((i) => e.unit.elements[i]),
    ]);
    expect(picked).toEqual([
      [m0.modelNum, ...[0, 1, 2].map((i) => u0.elements[i])],
      [m1.modelNum, ...[5, 6, 7].map((i) => u1.elements[i])],
    ]);
  });
});
