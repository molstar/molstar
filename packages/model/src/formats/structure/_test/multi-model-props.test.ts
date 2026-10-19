/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../mmcif.js';
import { ModelSymmetry } from '../property/symmetry.js';
import { IndexPairBonds } from '../property/bonds/index-pair.js';
import type { Model } from '@molstar/model/model/structure';

/** Two models of one residue (same topology, so model 2 shares model 1's hierarchy), a cell, an assembly and two bonds each. */
function mmcif() {
  const atoms: string[] = [];
  for (const model of [1, 2]) {
    ['N', 'CA', 'C'].forEach((name, i) => {
      const id = 3 * (model - 1) + i + 1;
      atoms.push(`ATOM ${id} ${name[0]} ${name} ALA A 1 1 ${i + model} 0 0 1 10 1 ALA A ${name} ${model}`);
    });
  }
  return `data_T
_cell.entry_id T
_cell.length_a 10
_cell.length_b 20
_cell.length_c 30
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.entry_id T
_symmetry.space_group_name_H-M 'P 1 21 1'
_pdbx_struct_assembly.id 1
_pdbx_struct_assembly.details author_defined_assembly
_pdbx_struct_assembly.oligomeric_details monomeric
_pdbx_struct_assembly.oligomeric_count 1
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression 1
_pdbx_struct_assembly_gen.asym_id_list A
_pdbx_struct_oper_list.id 1
_pdbx_struct_oper_list.type 'identity operation'
_pdbx_struct_oper_list.matrix[1][1] 1
_pdbx_struct_oper_list.matrix[1][2] 0
_pdbx_struct_oper_list.matrix[1][3] 0
_pdbx_struct_oper_list.vector[1] 0
_pdbx_struct_oper_list.matrix[2][1] 0
_pdbx_struct_oper_list.matrix[2][2] 1
_pdbx_struct_oper_list.matrix[2][3] 0
_pdbx_struct_oper_list.vector[2] 0
_pdbx_struct_oper_list.matrix[3][1] 0
_pdbx_struct_oper_list.matrix[3][2] 0
_pdbx_struct_oper_list.matrix[3][3] 1
_pdbx_struct_oper_list.vector[3] 0
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
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
${atoms.join('\n')}
loop_
_molstar_bond_site.atom_id_1
_molstar_bond_site.atom_id_2
_molstar_bond_site.value_order
_molstar_bond_site.type_id
1 2 sing covale
2 3 sing covale
4 5 sing covale
5 6 sing covale
`;
}

async function models(): Promise<Model[]> {
  const parsed = await CIF.parseText(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  const ret: Model[] = [];
  for (let i = 0; i < trajectory.frameCount; i++) ret.push(trajectory.getFrameAtIndex(i) as Model);
  return ret;
}

// Models with the topology of the previous one share its `customProperties`, not its dynamic property data
describe('dynamic format properties of a multi-model file', () => {
  it('the models share customProperties', async () => {
    const [m1, m2] = await models();
    expect(m2.atomicHierarchy).toBe(m1.atomicHierarchy);
    expect(m2.customProperties).toBe(m1.customProperties);
    expect(m2._dynamicPropertyData).not.toBe(m1._dynamicPropertyData);
  });

  it('every model has its symmetry, in either order', async () => {
    const [m1, m2] = await models();
    expect(ModelSymmetry.Provider.get(m1)?.assemblies.length).toBe(1);
    expect(ModelSymmetry.Provider.get(m2)?.assemblies.length).toBe(1);
    expect(ModelSymmetry.Provider.get(m2)?.spacegroup.name).toBe('P 1 21 1');

    const [n1, n2] = await models();
    expect(ModelSymmetry.Provider.get(n2)?.assemblies.length).toBe(1);
    expect(ModelSymmetry.Provider.get(n1)?.assemblies.length).toBe(1);
  });

  it('every model has its index pair bonds', async () => {
    const [m1, m2] = await models();
    // each model reads its own rows of the bond site: 2 bonds, stored in both directions
    for (const m of [m1, m2]) {
      const { bonds } = IndexPairBonds.Provider.get(m)!;
      expect(bonds.edgeCount).toBe(2);
      expect(Array.from(bonds.a)).toEqual([0, 1, 1, 2]);
      expect(Array.from(bonds.b)).toEqual([1, 0, 2, 1]);
    }
  });

  it('a value that was set stays, an unset one is obtained again per model', async () => {
    const [m1, m2] = await models();
    const set = ModelSymmetry.Provider.get(m1)!;
    const other = { ...set, assemblies: [] };
    ModelSymmetry.Provider.set(m2, other);
    expect(ModelSymmetry.Provider.get(m2)).toBe(other);
    expect(ModelSymmetry.Provider.get(m1)).toBe(set);
  });
});
