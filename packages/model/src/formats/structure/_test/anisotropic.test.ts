/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../mmcif.js';
import { AtomSiteAnisotrop } from '../property/anisotropic.js';
import type { Model } from '@molstar/model/model/structure';

/**
 * Two models of one residue (model 2 shares model 1's hierarchy), atom ids from `firstId`, and an ADP for every atom
 * with U11 = its atom id.
 */
function mmcif(firstId: number) {
  const atoms: string[] = [];
  const aniso: string[] = [];
  for (const model of [1, 2]) {
    ['N', 'CA', 'C'].forEach((name, i) => {
      const id = firstId + 3 * (model - 1) + i;
      atoms.push(`ATOM ${id} ${name[0]} ${name} ALA A 1 1 ${i + model} 0 0 1 10 1 ALA A ${name} ${model}`);
      aniso.push(`${id} ${id} 0.1 0.1 0 0 0`);
    });
  }
  return `data_T
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
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
${aniso.join('\n')}
`;
}

async function models(firstId: number): Promise<Model[]> {
  const parsed = await CIF.parse(mmcif(firstId)).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  const out: Model[] = [];
  for (let i = 0; i < trajectory.frameCount; i++) out.push(await trajectory.getFrameAtIndex(i));
  return out;
}

/** U11 of each element's ADP (its atom id), -1 without one */
function u11(model: Model) {
  const { data, elementToAnsiotrop } = AtomSiteAnisotrop.Provider.get(model)!;
  return Array.from(elementToAnsiotrop, (ai) => (ai < 0 ? -1 : data.U.value(ai)[0]));
}

describe('atom_site_anisotrop mapping (perf-audit-2 04/B3)', () => {
  it('maps every atom of each model, the last one included', async () => {
    const [m1, m2] = await models(1);
    expect(m2.atomicHierarchy).toBe(m1.atomicHierarchy);
    expect(u11(m1)).toEqual([1, 2, 3]);
    expect(u11(m2)).toEqual([4, 5, 6]);
  });

  it('maps atom ids that do not start at 1', async () => {
    const [m1, m2] = await models(101);
    expect(u11(m1)).toEqual([101, 102, 103]);
    expect(u11(m2)).toEqual([104, 105, 106]);
  });
});
