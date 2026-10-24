/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CIF } from '@molstar/io/reader/cif';
import { MolScriptBuilder as B } from '@molstar/query-language/language/builder';
import { describe, expect, it } from 'vitest';
import { trajectoryFromMmCIF } from '../../../../formats/structure/mmcif.js';
import { Script } from '../../../../script/script.js';
import { Structure, StructureSelection } from '../../../structure.js';

// `queryInSelection` in the complement (MolQL's `not`, as the PyMOL and VMD transpilers write it): the complement of
// nothing is everything.

const mmcif = `data_T
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
ATOM 1 N N ALA A 1 1 0 0 0 1 A 1
ATOM 2 C CA ALA A 1 1 1 0 0 1 A 1
ATOM 3 O O ALA A 1 1 2 0 0 1 A 1
`;

async function structure() {
  const parsed = await CIF.parse(mmcif).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  return Structure.ofModel(model);
}

const count = (s: Structure, test: ReturnType<typeof B.core.rel.eq>) =>
  StructureSelection.unionStructure(
    Script.getStructureSelection(
      B.struct.generator.queryInSelection({
        0: B.struct.generator.atomGroups({ 'atom-test': test }),
        query: B.struct.generator.all(),
        'in-complement': true,
      }),
      s,
    ),
  ).elementCount;

describe('queryInSelection in the complement', () => {
  it('is everything when the selection is empty', async () => {
    const s = await structure();
    expect(count(s, B.core.rel.eq([B.acp('elementSymbol'), B.es('H')]))).toBe(3);
  });

  it('is the rest otherwise', async () => {
    const s = await structure();
    expect(count(s, B.core.rel.eq([B.acp('elementSymbol'), B.es('O')]))).toBe(2);
  });
});
