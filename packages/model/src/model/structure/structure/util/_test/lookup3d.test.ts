/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Queries, StructureProperties, StructureQuery, StructureSelection } from '@molstar/model/model/structure';
import { Structure } from '../../structure.js';

/** Two chains, two units: A has one atom at x = `a`, B one at x = `b` */
async function twoUnits(a: number, b: number) {
  const cif = `data_T
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 1 ${a} 0 0 1 A 1
ATOM 2 C CA ALA B 2 1 ${b} 0 0 1 B 1
`;
  const parsed = await CIF.parse(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

describe('StructureLookup3D.approxNearest', () => {
  it('gives the nearest over the units, once (perf-audit-2 02/B1)', async () => {
    // the nearer atom in either unit, on either side: some of these search the farther unit first
    for (const [a, b] of [
      [0.3, 2.9],
      [2.9, 0.3],
      [0.3, -2.9],
      [-2.9, 0.3],
    ]) {
      const structure = await twoUnits(a, b);
      expect(structure.units.length).toBe(2);

      const r = structure.lookup3d.approxNearest(0, 0, 0, 5);
      expect(r.count).toBe(1);
      const unit = r.units[0];
      expect(structure.model.atomicConformation.x[unit.elements[r.indices[0]]]).toBeCloseTo(0.3, 5);
    }
  });
});

describe('include surroundings with an atom radius', () => {
  it('reaches units beyond the radius but within the atoms radii (perf-audit-2 05/B12)', async () => {
    // vdW 1.7 each: the surfaces are 0.1 apart, inside the 0.5 radius; B's unit sphere is 3.5 from A
    const structure = await twoUnits(0, 3.5);
    const query = Queries.modifiers.includeSurroundings(
      Queries.generators.atoms({ chainTest: (ctx) => StructureProperties.chain.label_asym_id(ctx.element) === 'A' }),
      { radius: 0.5, elementRadius: (ctx) => StructureProperties.atom.vdw_radius(ctx.element) },
    );
    const s = StructureSelection.unionStructure(StructureQuery.run(query, structure));
    expect(s.elementCount).toBe(2);
  });
});
