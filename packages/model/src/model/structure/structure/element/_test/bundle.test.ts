/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { ElementIndex } from '../../../model/indexing.js';
import { Structure } from '../../structure.js';
import { StructureSymmetry } from '../../symmetry.js';
import type { Unit } from '../../unit.js';
import { Bundle } from '../bundle.js';
import type { UnitIndex } from '../element.js';
import { Loci } from '../loci.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

async function crystalCell() {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1crn.cif'), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  const s = Structure.ofModel(model);
  return StructureSymmetry.buildSymmetryRange(s, Vec3.create(0, 0, 0), Vec3.create(0, 0, 0)).run();
}

function firstElements(unit: Unit, n: number) {
  return SortedArray.ofSortedArray<ElementIndex>(unit.elements.slice(0, n) as any);
}

describe('Bundle', () => {
  it('toStructure maps each copy of a unit through its own elements (perf-audit-2 05/B10)', async () => {
    const cell = await crystalCell();
    const [u1, u2] = cell.units;
    expect(u1.invariantId).toBe(u2.invariantId);

    // copies of one invariant unit holding different elements: 10 and 20 atoms
    const child = Structure.create([u1.getChild(firstElements(u1, 10)), u2.getChild(firstElements(u2, 20))], {
      parent: cell,
    });
    const indices = OrderedSet.ofBounds(0, 10) as OrderedSet<UnitIndex>;
    const loci = Loci(child, [
      { unit: child.units[0], indices },
      { unit: child.units[1], indices },
    ]);

    const bundle = Bundle.fromLoci(loci);
    expect(bundle.elements.length).toBe(1);
    expect(bundle.elements[0].groupedUnits[0].length).toBe(2);

    const s = Bundle.toStructure(bundle, child);
    expect(s.elementCount).toBe(20);
    expect(Array.from(s.units[1].elements)).toEqual(Array.from(u2.elements.slice(0, 10)));
    expect(Loci.size(Bundle.toLoci(bundle, child))).toBe(20);
  });
});
