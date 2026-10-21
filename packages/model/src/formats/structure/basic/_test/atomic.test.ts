/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../../mmcif.js';

/** A one-model mmCIF of these rows: [asym, entity, seq, comp, atom]. */
function mmcif(rows: [string, number, number | '.', string, string][]) {
  const lines = rows.map(([asym, entity, seq, comp, atom], i) => {
    const authSeq = seq === '.' ? 1 : seq;
    return `ATOM ${i + 1} C ${atom} ${comp} ${asym} ${entity} ${seq} ${i} 0 0 ${authSeq} ${comp} ${asym} ${atom} 1`;
  });
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
${lines.join('\n')}
`;
}

async function residueSourceIndex(rows: Parameters<typeof mmcif>[0]) {
  const parsed = await CIF.parse(mmcif(rows)).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  return Array.from(model.atomicHierarchy.residueSourceIndex.toArray());
}

describe('residue source index', () => {
  it('is the identity for residues in file order', async () => {
    const index = await residueSourceIndex([
      ['A', 1, 1, 'ALA', 'N'],
      ['A', 1, 1, 'ALA', 'CA'],
      ['A', 1, 2, 'GLY', 'N'],
      ['B', 1, 1, 'SER', 'N'],
      ['C', 2, '.', 'ZN', 'ZN'],
    ]);
    expect(index).toEqual([0, 1, 2, 3]);
  });

  it("is each residue's rank by its first row in the file when the model sorts them", async () => {
    // Chain A's seq ids reversed: sorted (entities in their first rows' order), residues are C, A1, A2, B1.
    const index = await residueSourceIndex([
      ['C', 2, '.', 'ZN', 'ZN'],
      ['A', 1, 2, 'GLY', 'N'],
      ['A', 1, 1, 'ALA', 'N'],
      ['A', 1, 1, 'ALA', 'CA'],
      ['B', 1, 1, 'SER', 'N'],
    ]);
    expect(index).toEqual([0, 2, 1, 3]);
  });
});
