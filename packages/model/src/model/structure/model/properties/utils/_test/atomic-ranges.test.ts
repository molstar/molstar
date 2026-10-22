/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../../../../../../formats/structure/mmcif.js';
import type { AtomicConformation } from '../../atomic/conformation.js';
import { getAtomicRanges } from '../atomic-ranges.js';

/** A one-model mmCIF of DNA chain A: [seq, comp, atom, x]. */
function mmcif(rows: [number, string, string, number][]) {
  const lines = rows.map(
    ([seq, comp, atom, x], i) => `ATOM ${i + 1} ${atom[0]} "${atom}" ${comp} A 1 ${seq} ${x} 0 0 ${seq} ${comp} A 1`,
  );
  return `data_T
loop_
_entity.id
_entity.type
1 polymer
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
_atom_site.pdbx_PDB_model_num
${lines.join('\n')}
`;
}

async function model(rows: Parameters<typeof mmcif>[0]) {
  const parsed = await CIF.parse(mmcif(rows)).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

/** coordinate arrays that throw when read outside their elements */
function strict(conformation: AtomicConformation): AtomicConformation {
  const guard = (a: ArrayLike<number>) =>
    new Proxy(a, {
      get(target, p, receiver) {
        if (typeof p === 'string' && /^-?\d+$/.test(p) && (+p < 0 || +p >= target.length)) {
          throw new Error(`coordinate ${p} read`);
        }
        return Reflect.get(target, p, receiver);
      },
    });
  return { ...conformation, x: guard(conformation.x), y: guard(conformation.y), z: guard(conformation.z) };
}

describe('atomic ranges', () => {
  it('never read coordinates of a missing backbone atom (perf-audit-2 02/B5)', async () => {
    // the 5' nucleotide has no P (backbone start and coarse backbone): the cyclic check and the link to the next
    // residue have no atom to measure from; the others are connected (O3'–P 1.6 Å)
    const m = await model([
      [1, 'DA', "O3'", 0],
      [2, 'DT', 'P', 1.6],
      [2, 'DT', "O3'", 5],
      [3, 'DG', 'P', 6.6],
      [3, 'DG', "O3'", 10],
    ]);
    const { atomicHierarchy: h, entities, sequence } = m;
    const ranges = getAtomicRanges(h, entities, strict(m.atomicConformation), sequence);
    expect(ranges.cyclicPolymerMap.size).toBe(0);
    expect(Array.from(ranges.polymerRanges)).toEqual(Array.from(m.atomicRanges.polymerRanges));
    expect(Array.from(ranges.gapRanges)).toEqual(Array.from(m.atomicRanges.gapRanges));
  });
});
