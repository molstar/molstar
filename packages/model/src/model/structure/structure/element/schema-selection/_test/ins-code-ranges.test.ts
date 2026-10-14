/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '../../../../structure.js';
import type { SchemaItem } from '../../schema.js';
import { IndicesAndSortings } from '../indexing.js';
import { substructureFromRows } from '../row-structure.js';
import { getAtomRangesForRow } from '../selections.js';

/** The residues are given in the order of the file, as `auth_asym_id:auth_seq_id[ins code]`; every residue has 2 atoms (N, CA). */
const Chains: Record<string, string[]> = {
  // an insertion series after 11 (11A, 11B), and one after 13
  A: ['10', '11', '11A', '11B', '12', '13', '13A', '14'],
  // descending numbers
  B: ['20', '19', '19A', '18'],
  // only the residue with an insertion code is there (11 is missing)
  C: ['10', '11A', '12'],
  // a ligand chain, which belongs to the auth chain A (as the ligands of a polymer chain do)
  L: ['500', '501'],
  // the waters of the auth chain A, numbered again from the start (as in many files), and one which only they have
  W: ['10', '11', '12', '13', '14', '700'],
};
/** The `auth_asym_id` of the chains which are not their own auth chain */
const AuthOf: Record<string, string> = { L: 'A', W: 'A' };

function mmcif(): string {
  const lines = [
    'data_T',
    'loop_',
    ...'group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'
      .split(' ')
      .map((c) => `_atom_site.${c}`),
  ];
  let id = 1;
  for (const [asym, residues] of Object.entries(Chains)) {
    residues.forEach((residue, i) => {
      const [, seq, ins] = /^(\d+)([A-Z]?)$/.exec(residue)!;
      const comp = asym === 'W' ? 'HOH' : 'ALA';
      ['N', 'CA'].forEach((name, a) => {
        lines.push(
          `ATOM ${id++} ${name[0]} ${name} . ${comp} ${asym} ${asym === 'W' ? 2 : 1} ${i + 1} ${ins || '?'} ${a} ${i} 0 1 10 ${seq} ${comp} ${AuthOf[asym] ?? asym} ${name} 1`,
        );
      });
    });
  }
  return `${lines.join('\n')}\n`;
}

describe('ranges of auth_seq_id with insertion codes', async () => {
  const parsed = await CIF.parseText(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
  const structure = Structure.ofModel(model);
  const indices = IndicesAndSortings.get(model);
  const { residueAtomSegments } = model.atomicHierarchy;

  /** The selected residues as `chain:residue`, in the order of the structure, and the atoms they have (2 each). */
  function select(row: SchemaItem, atomsPerResidue = 2) {
    const result = substructureFromRows(structure, [row]);
    const residues: string[] = [];
    let nAtoms = 0;
    for (const unit of result.units) {
      for (let i = 0; i < unit.elements.length; i++) {
        nAtoms++;
        const iRes = residueAtomSegments.index[unit.elements[i]];
        const label = `${model.atomicHierarchy.chains.auth_asym_id.value(model.atomicHierarchy.chainAtomSegments.index[unit.elements[i]])}:${model.atomicHierarchy.residues.auth_seq_id.value(iRes)}${model.atomicHierarchy.residues.pdbx_PDB_ins_code.value(iRes)}`;
        if (residues[residues.length - 1] !== label) residues.push(label);
      }
    }
    expect(nAtoms).toBe(residues.length * atomsPerResidue);
    // (the range of the row is the same when it is resolved on its own)
    const ranges = getAtomRangesForRow(row, model, { instanceId: '1_555', name: '1_555' }, indices);
    let n = 0;
    if (ranges) for (let i = 0; i < ranges.from.length; i++) n += ranges.to[i] - ranges.from[i];
    expect(n).toBe(nAtoms);
    return residues;
  }
  const chainA = (...residues: string[]) => residues.map((r) => `A:${r}`);

  it('is a range of residues in the order of the file, between the first residue of the start and the last of the end', () => {
    const A = { label_asym_id: 'A' };
    expect(select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: 'A' })).toEqual(
      chainA('11A', '11B', '12', '13', '13A', '14'),
    );
    expect(select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: 'B', end_auth_seq_id: 13 })).toEqual(
      chainA('11B', '12', '13', '13A'),
    );
    expect(
      select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 13, end_auth_ins_code: 'A' }),
    ).toEqual(chainA('11A', '11B', '12', '13', '13A'));
    expect(select({ ...A, beg_auth_seq_id: 10, end_auth_seq_id: 11, end_auth_ins_code: 'A' })).toEqual(
      chainA('10', '11', '11A'),
    );
    expect(
      select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 11, end_auth_ins_code: 'A' }),
    ).toEqual(chainA('11A'));
  });

  it('an empty insertion code is the residue without one', () => {
    const A = { label_asym_id: 'A' };
    expect(select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: '', end_auth_seq_id: 12 })).toEqual(
      chainA('11', '11A', '11B', '12'),
    );
    // (the end without a code is the last residue of the number, the end with the empty code is the residue without a code)
    expect(select({ ...A, beg_auth_seq_id: 12, end_auth_seq_id: 13 })).toEqual(chainA('12', '13', '13A'));
    expect(select({ ...A, beg_auth_seq_id: 12, end_auth_seq_id: 13, end_auth_ins_code: '' })).toEqual(
      chainA('12', '13'),
    );
    expect(
      select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: '', end_auth_seq_id: 11, end_auth_ins_code: '' }),
    ).toEqual(chainA('11'));
    // (a start without a code is the first residue of the number)
    expect(select({ ...A, beg_auth_seq_id: 13, end_auth_ins_code: 'A', end_auth_seq_id: 13 })).toEqual(
      chainA('13', '13A'),
    );
  });

  it('follows the order of the file, also when the numbers descend', () => {
    const B = { label_asym_id: 'B' };
    expect(select({ ...B, beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 18 })).toEqual([
      'B:19',
      'B:19A',
      'B:18',
    ]);
    expect(select({ ...B, beg_auth_seq_id: 19, beg_auth_ins_code: 'A', end_auth_seq_id: 18 })).toEqual([
      'B:19A',
      'B:18',
    ]);
    expect(select({ ...B, beg_auth_seq_id: 20, end_auth_seq_id: 19, end_auth_ins_code: '' })).toEqual(['B:20', 'B:19']);
    // (without insertion codes it is a range of values, which is empty when it descends)
    expect(select({ ...B, beg_auth_seq_id: 19, end_auth_seq_id: 18 })).toEqual([]);
  });

  it('selects nothing from a chain which does not have a bound, or for bounds in the reverse order', () => {
    const A = { label_asym_id: 'A' };
    expect(select({ ...A, beg_auth_seq_id: 11, beg_auth_ins_code: 'C' })).toEqual([]);
    expect(select({ ...A, beg_auth_seq_id: 99, beg_auth_ins_code: 'A' })).toEqual([]);
    expect(select({ ...A, beg_auth_seq_id: 10, end_auth_seq_id: 99, end_auth_ins_code: 'A' })).toEqual([]);
    expect(
      select({ ...A, beg_auth_seq_id: 13, beg_auth_ins_code: 'A', end_auth_seq_id: 11, end_auth_ins_code: 'A' }),
    ).toEqual([]);
    expect(select({ label_asym_id: 'B', beg_auth_seq_id: 18, end_auth_seq_id: 20, end_auth_ins_code: '' })).toEqual([]);
    // chain C has no residue 11 without a code: a bound without a code is the number, which it does not have
    expect(select({ label_asym_id: 'C', beg_auth_seq_id: 11, beg_auth_ins_code: '' })).toEqual([]);
  });

  it('an insertion code without the number of its bound selects nothing (not everything)', () => {
    expect(select({ beg_auth_ins_code: 'A' })).toEqual([]);
    expect(select({ label_asym_id: 'A', end_auth_ins_code: 'A' })).toEqual([]);
    expect(select({ label_asym_id: 'A', beg_auth_seq_id: 10, end_auth_ins_code: 'A' })).toEqual([]);
    expect(select({ label_asym_id: 'A', beg_auth_ins_code: '', end_auth_seq_id: 12 })).toEqual([]);
  });

  it('is resolved in each chain which the row selects', () => {
    // A: 11A..13A, B: no 11, C: has 11A but no 13A
    expect(
      select({ beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 13, end_auth_ins_code: 'A' }),
    ).toEqual(chainA('11A', '11B', '12', '13', '13A'));
    // C has 11A and 12, A has them too
    expect(select({ beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 12 })).toEqual([
      ...chainA('11A', '11B', '12'),
      'C:11A',
      'C:12',
    ]);
    // the chain is selected by its auth_asym_id and by an entity, too
    expect(select({ auth_asym_id: 'C', beg_auth_seq_id: 11, beg_auth_ins_code: 'A' })).toEqual(['C:11A', 'C:12']);
    expect(select({ label_entity_id: '1', beg_auth_seq_id: 12, end_auth_seq_id: 12, end_auth_ins_code: '' })).toEqual([
      'A:12',
      'C:12',
    ]);
  });

  it('is resolved in the auth chain: all the chains with the auth_asym_id (a polymer and its ligands), in the order of the file', () => {
    // (chain L, the ligands, has the auth_asym_id A: a range which starts in the polymer and ends in the ligands is one range)
    const range = { beg_auth_seq_id: 13, beg_auth_ins_code: 'A', end_auth_seq_id: 501, end_auth_ins_code: '' };
    expect(select({ auth_asym_id: 'A', ...range })).toEqual(['A:13A', 'A:14', 'A:500', 'A:501']);
    // (the chains which the row selects get the part of the range which is theirs)
    expect(select({ label_asym_id: 'A', ...range })).toEqual(chainA('13A', '14'));
    expect(select({ label_asym_id: 'L', ...range })).toEqual(['A:500', 'A:501']);
    // an open end is the end of the auth chain, an open beginning is its beginning
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 14, end_auth_ins_code: '', end_auth_seq_id: 500 })).toEqual([
      'A:14',
      'A:500',
    ]);
    expect(select({ auth_asym_id: 'A', end_auth_seq_id: 10, end_auth_ins_code: '' })).toEqual(['A:10']);
    // bounds in the reverse order (the ligands come after the polymer in the file), or which the auth chain does not have, select nothing
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 500, end_auth_seq_id: 11, end_auth_ins_code: 'B' })).toEqual(
      [],
    );
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 14, beg_auth_ins_code: 'A' })).toEqual([]);
    expect(select({ auth_asym_id: 'B', beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 501 })).toEqual([]);
    // (the ligands of A do not make a range of chain B)
    expect(select({ auth_asym_id: 'B', beg_auth_seq_id: 19, beg_auth_ins_code: '' })).toEqual([
      'B:19',
      'B:19A',
      'B:18',
    ]);
  });

  it('looks for a bound in the chains which are not waters first (waters are numbered again from the start)', () => {
    // (the waters W:10..14 have numbers of the polymer, and come after it and the ligands in the file: without the rule the end 11 is the water 11)
    expect(select({ label_asym_id: 'A', beg_auth_seq_id: 10, end_auth_seq_id: 11, end_auth_ins_code: '' })).toEqual(
      chainA('10', '11'),
    );
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 10, end_auth_seq_id: 11, end_auth_ins_code: '' })).toEqual(
      chainA('10', '11'),
    );
    expect(
      select({
        auth_asym_id: 'A',
        beg_auth_seq_id: 12,
        beg_auth_ins_code: '',
        end_auth_seq_id: 14,
        end_auth_ins_code: '',
      }),
    ).toEqual(chainA('12', '13', '13A', '14'));
    // the waters get no part of such a range, even when the row asks for them
    expect(select({ label_asym_id: 'W', beg_auth_seq_id: 10, end_auth_seq_id: 11, end_auth_ins_code: '' })).toEqual([]);
    // a bound which only the waters have is found in them (the range then goes from the polymer to the water, through the ligands)
    const toWater = { beg_auth_seq_id: 14, beg_auth_ins_code: '', end_auth_seq_id: 700, end_auth_ins_code: '' };
    const all = ['A:14', 'A:500', 'A:501', 'A:10', 'A:11', 'A:12', 'A:13', 'A:14', 'A:700'];
    expect(select({ auth_asym_id: 'A', ...toWater })).toEqual(all);
    expect(select({ label_asym_id: 'W', ...toWater })).toEqual(all.slice(3));
    expect(select({ label_asym_id: 'L', ...toWater })).toEqual(['A:500', 'A:501']);
    // (an open end is the end of the auth chain: the waters)
    expect(select({ auth_asym_id: 'A', beg_auth_seq_id: 501, beg_auth_ins_code: '' })).toEqual([
      'A:501',
      ...all.slice(3),
    ]);
  });

  it('the other criteria of the row filter the residues of the range', () => {
    const range = { label_asym_id: 'A', beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 13 };
    expect(select({ ...range, label_seq_id: 5 })).toEqual(chainA('12')); // (the 5th residue of A is 12)
    expect(select({ ...range, label_seq_id: 1 })).toEqual([]);
    expect(select({ ...range, pdbx_PDB_ins_code: '' })).toEqual(chainA('12', '13'));
    expect(select({ ...range, pdbx_PDB_ins_code: 'B' })).toEqual(chainA('11B'));
    expect(select({ ...range, residue_index: 2 })).toEqual(chainA('11A'));
    expect(select({ ...range, auth_seq_id: 13 })).toEqual(chainA('13', '13A'));
    expect(select({ ...range, label_atom_id: 'CA' }, 1)).toEqual(chainA('11A', '11B', '12', '13', '13A'));
    expect(select({ ...range, beg_label_seq_id: 4, end_label_seq_id: 6 })).toEqual(chainA('11B', '12', '13'));
    expect(select({ ...range, label_comp_id: 'GLY' })).toEqual([]);
  });

  it('applies to the atom given by atom_id and atom_index, too', () => {
    const range = { beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 12 };
    // atoms of A: 2 per residue, ids from 1: residue 11A is the 3rd residue (ids 5, 6), residue 13A the 7th
    expect(select({ ...range, label_asym_id: 'A', atom_id: 5 }, 1)).toEqual(chainA('11A'));
    expect(select({ ...range, atom_index: 4 }, 1)).toEqual(chainA('11A')); // (the 3rd residue, the first atom)
    expect(select({ ...range, atom_id: 3 }, 1)).toEqual([]); // (the 2nd residue, 11: before the range)
    expect(select({ ...range, atom_id: 13 }, 1)).toEqual([]); // (13A: after the range)
    expect(select({ ...range, atom_index: 12 }, 1)).toEqual([]);
  });

  it('is resolved in the chain of the residue which another criterion finds in the whole model', () => {
    // (residue_index and pdbx_PDB_ins_code are looked up in all chains at once; the range is then still the one of the chain of the residue)
    const iB19A = 8 + 2; // (the 3rd residue of B, after the 8 residues of A)
    expect(select({ beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 18, residue_index: iB19A })).toEqual([
      'B:19A',
    ]);
    expect(
      select({ beg_auth_seq_id: 11, beg_auth_ins_code: 'A', end_auth_seq_id: 12, pdbx_PDB_ins_code: 'A' }),
    ).toEqual(['A:11A', 'C:11A']);
    expect(select({ beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 18, residue_index: 3 })).toEqual([]);
  });

  it('applies to an atom of a range which descends in values', () => {
    // atoms of B start at id 17: B:20 (17, 18), B:19 (19, 20), B:19A (21, 22), B:18 (23, 24)
    expect(select({ beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 18, atom_id: 21 }, 1)).toEqual([
      'B:19A',
    ]);
    expect(select({ beg_auth_seq_id: 19, beg_auth_ins_code: '', end_auth_seq_id: 18, atom_id: 17 }, 1)).toEqual([]);
  });

  it('a range without insertion codes is still a range of values (it includes all the insertion codes of its numbers)', () => {
    expect(select({ label_asym_id: 'A', beg_auth_seq_id: 11, end_auth_seq_id: 13 })).toEqual(
      chainA('11', '11A', '11B', '12', '13', '13A'),
    );
    expect(select({ label_asym_id: 'A', end_auth_seq_id: 11 })).toEqual(chainA('10', '11', '11A', '11B'));
  });
});
