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
import {
  Queries,
  type QueryContext,
  Structure,
  StructureProperties,
  StructureQuery,
  StructureSelection,
} from '@molstar/model/model/structure';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

async function load1crn() {
  const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1crn.cif'), 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

describe('pick', () => {
  it('restores the current structure after a nested pick (perf-audit-2 05/B6)', async () => {
    const structure = await load1crn();
    const inner = Queries.filters.pick(Queries.generators.atoms(), () => true);
    // per residue: runs the inner pick, then reads its own residue
    const outer = Queries.filters.pick(Queries.generators.residues(), (ctx) => {
      const count = StructureSelection.structureCount(inner(ctx));
      return count > 0 && Queries.atomset.atomCount(ctx) >= 6;
    });
    const sel = StructureQuery.run(outer, structure);
    const residues = StructureSelection.structureCount(StructureQuery.run(Queries.generators.residues(), structure));
    const big = StructureSelection.structureCount(
      StructureQuery.run(
        Queries.filters.pick(Queries.generators.residues(), (ctx) => Queries.atomset.atomCount(ctx) >= 6),
        structure,
      ),
    );
    expect(big).toBeGreaterThan(0);
    expect(big).toBeLessThan(residues);
    expect(StructureSelection.structureCount(sel)).toBe(big);
  });

  it('restores the current structure after a nested withSameAtomProperties (perf-audit-2 05/B6)', async () => {
    const structure = await load1crn();
    const props = (ctx: QueryContext) => StructureProperties.residue.label_comp_id(ctx.element);
    const inner = Queries.filters.withSameAtomProperties(
      Queries.generators.residues(),
      Queries.generators.residues({ residueTest: (ctx) => props(ctx) === 'ALA' }),
      props,
    );
    // per residue: runs the inner filter, then reads its own residue
    const outer = Queries.filters.pick(Queries.generators.residues(), (ctx) => {
      StructureSelection.structureCount(inner(ctx));
      return Queries.atomset.atomCount(ctx) >= 6;
    });
    const big = Queries.filters.pick(Queries.generators.residues(), (ctx) => Queries.atomset.atomCount(ctx) >= 6);
    expect(StructureSelection.structureCount(StructureQuery.run(outer, structure))).toBe(
      StructureSelection.structureCount(StructureQuery.run(big, structure)),
    );
  });
});
