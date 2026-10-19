/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { Model } from '@molstar/model/model/structure/model';
import { BondType } from '@molstar/model/model/structure/model/types';
import { Structure, Unit } from '@molstar/model/model/structure';
import { trajectoryFromMmCIF } from '../mmcif.js';

// Entity 1, a polymer, in chains A (ALA GLY) and B (ALA GLY SER); entity 2, a zinc ligand (ZN bonded to N1 in the
// component dictionary), in chain C. No entity_poly_seq, so sequences come from the hierarchy.
const atoms = [
  ['ATOM', 'N', 'N', 'ALA', 'A', 1, 1],
  ['ATOM', 'C', 'CA', 'ALA', 'A', 1, 1],
  ['ATOM', 'N', 'N', 'GLY', 'A', 1, 2],
  ['ATOM', 'C', 'CA', 'GLY', 'A', 1, 2],
  ['ATOM', 'N', 'N', 'ALA', 'B', 1, 1],
  ['ATOM', 'C', 'CA', 'ALA', 'B', 1, 1],
  ['ATOM', 'N', 'N', 'GLY', 'B', 1, 2],
  ['ATOM', 'C', 'CA', 'GLY', 'B', 1, 2],
  ['ATOM', 'N', 'N', 'SER', 'B', 1, 3],
  ['ATOM', 'C', 'CA', 'SER', 'B', 1, 3],
  ['HETATM', 'ZN', 'ZN', 'ZNL', 'C', 2, '.'],
  ['HETATM', 'N', 'N1', 'ZNL', 'C', 2, '.'],
] as const;

function mmcif() {
  const rows = atoms.map(([group, type, name, comp, asym, entity, seq], i) => {
    const x = 1.5 * i;
    const authSeq = seq === '.' ? 1 : seq;
    return `${group} ${i + 1} ${type} ${name} ${comp} ${asym} ${entity} ${seq} ${x} 0 0 1 10 ${authSeq} ${comp} ${asym} ${name} 1`;
  });
  return `data_T
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
ZNL ZN N1 SING N
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
${rows.join('\n')}
`;
}

async function model() {
  const parsed = await CIF.parse(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

describe('mmCIF model', () => {
  it("a polymer entity's sequence is its first chain's", async () => {
    const m = await model();
    const seq = m.sequence.byEntityKey[m.entities.getEntityIndex('1')];
    expect(seq.sequence.code.toArray().join('')).toBe('AG');
    expect(m.sequence.sequences.length).toBe(1);
  });

  it('the inverted atom source index of a file-ordered model is the identity', async () => {
    const m = await model();
    const { isIdentity, invertedIndex } = Model.getInvertedAtomSourceIndex(m);
    expect(Array.from(invertedIndex)).toEqual(atoms.map((_, i) => i));
    expect(isIdentity).toBe(true);
  });

  it('a metal bond from the component dictionary is metallic, not covalent', async () => {
    const m = await model();
    const ligand = Structure.ofModel(m).units.find((u) => u.elements.length === 2) as Unit.Atomic;
    const { edgeCount, edgeProps } = ligand.bonds;
    expect(edgeCount).toBe(1);
    expect(edgeProps.flags[0]).toBe(BondType.Flag.MetallicCoordination);
    expect(edgeProps.key[0]).toBeGreaterThanOrEqual(0);
  });
});
