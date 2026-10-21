/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { StructConn } from '../property/bonds/struct_conn.js';
import { trajectoryFromMmCIF } from '../mmcif.js';

// Entities 1 (chains A, B) and 2 (chain C, and a stray chain A row: one label_asym_id in two entities, as in a
// malformed file). struct_conn rows: within A, A to C, a partner on an atom that does not exist, a partner in an
// asym id that does not exist, the stray A's atom (only entity 2's A has it), an atom both A chains have (entity 1's
// comes first).
const atoms = [
  ['A', 1, 1, 'CYS', 'SG'],
  ['A', 1, 2, 'CYS', 'SG'],
  ['B', 1, 1, 'CYS', 'SG'],
  ['C', 2, '.', 'ZN', 'ZN'],
  ['A', 2, '.', 'ZN', 'ZN'],
  ['A', 2, '.', 'ZN', 'SG'],
] as const;

const conn = [
  ['disulf', 'A', 1, 'SG', 'A', 2, 'SG'],
  ['metalc', 'A', 1, 'SG', 'C', '.', 'ZN'],
  ['metalc', 'A', 1, 'XX', 'C', '.', 'ZN'],
  ['metalc', 'Q', 1, 'SG', 'C', '.', 'ZN'],
  ['metalc', 'A', '.', 'ZN', 'B', 1, 'SG'],
  ['metalc', 'A', 1, 'SG', 'A', '.', 'SG'],
] as const;

function mmcif() {
  const rows = atoms.map(([asym, entity, seq, comp, atom], i) => {
    const authSeq = seq === '.' ? 1 : seq;
    return `ATOM ${i + 1} C ${atom} ${comp} ${asym} ${entity} ${seq} ${i} 0 0 ${authSeq} ${comp} ${asym} ${atom} 1`;
  });
  const conns = conn.map(
    ([type, a1, s1, n1, a2, s2, n2], i) =>
      `c${i} ${type} ${a1} ${s1} ${s1 === '.' ? 1 : s1} ${n1} ${a2} ${s2} ${s2 === '.' ? 1 : s2} ${n2}`,
  );
  return `data_T
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
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
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
${rows.join('\n')}
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_label_atom_id
${conns.join('\n')}
`;
}

describe('struct_conn', () => {
  it("finds each partner's atom in the entity its asym id belongs to", async () => {
    const parsed = await CIF.parse(mmcif()).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const { atomSourceIndex } = model.atomicHierarchy;
    const entries = StructConn.Provider.get(model)!.entries.map((e) => [
      e.rowIndex,
      atomSourceIndex.value(e.partnerA.atomIndex),
      atomSourceIndex.value(e.partnerB.atomIndex),
    ]);
    // [row, source atom of partner A, of partner B]: rows 2 and 3 have no partner atom
    expect(entries).toEqual([
      [0, 0, 1],
      [1, 0, 3],
      [4, 4, 2],
      [5, 0, 0],
    ]);
  });
});
