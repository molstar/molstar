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
import {
  Queries,
  Structure,
  StructureProperties as P,
  StructureQuery,
  StructureSelection,
} from '@molstar/model/model/structure';
import { StructureSubsetBuilder } from '../../../structure/util/subset-builder.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

async function load1cbs() {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

// Chain-level atom groups (polymer, protein, ligand components) set a unit whose every chain passes as a whole, without
// copying its elements; a unit with some chains passing gets their ranges (perf-audit-2 05/P19).
describe('atomGroups at chain level', () => {
  it('takes whole units without copying their elements', async () => {
    const structure = await load1cbs();
    const addRange = vi.spyOn(StructureSubsetBuilder.prototype, 'addElementRange');
    const polymer = Queries.generators.atoms({ entityTest: (ctx) => P.entity.type(ctx.element) === 'polymer' });
    const sel = StructureQuery.run(polymer, structure);
    const calls = addRange.mock.calls.length;
    addRange.mockRestore();

    const s = StructureSelection.unionStructure(sel);
    const polymerUnits = structure.units.filter(
      (u) => P.entity.type({ structure, unit: u, element: u.elements[0] } as any) === 'polymer',
    );
    expect(polymerUnits.length).toBeGreaterThan(0);
    expect(s.units.length).toBe(polymerUnits.length);
    for (const u of s.units) expect(u).toBe(structure.unitMap.get(u.id));
    expect(calls).toBe(0);
  });

  it('selects by chain as a scan of the elements does', async () => {
    const structure = await load1cbs();
    // every unit but the first chain's passes
    const asyms = new Set(
      structure.units.map((u) => P.chain.label_asym_id({ structure, unit: u, element: u.elements[0] } as any)),
    );
    const [first] = asyms;
    const q = Queries.generators.atoms({ chainTest: (ctx) => P.chain.label_asym_id(ctx.element) !== first });
    const s = StructureSelection.unionStructure(StructureQuery.run(q, structure));
    let expected = 0;
    for (const u of structure.units) {
      for (const e of u.elements) {
        if (P.chain.label_asym_id({ structure, unit: u, element: e } as any) !== first) expected++;
      }
    }
    expect(s.elementCount).toBe(expected);
  });
});
