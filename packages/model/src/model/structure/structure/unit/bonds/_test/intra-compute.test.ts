/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { computeIntraUnitBonds } from '../intra-compute.js';
import type { BondComputationProps } from '../common.js';

// Pins the intra-unit bonds computed from distances and the component dictionary, as sorted (unit, a < b, order,
// flags, key) rows. Recorded before the 02/P3 and 05/P7 changes (.v6/perf-audit/), which must keep them; the edge
// order in the graph may change.

const examples = path.resolve(import.meta.dirname, '../../../../../../../../../data/examples');

async function parseCif(data: string): Promise<Model> {
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
}

const readExample = (name: string) => fs.readFileSync(path.join(examples, name), 'utf8');

/** The CA atoms of 1crn only: a coarse-grained model */
function caOnly(data: string) {
  return data
    .split('\n')
    .filter((l) => !l.startsWith('ATOM') || / CA {2}\. /.test(l))
    .join('\n');
}

/**
 * Two alanines with hydrogens (HB1, HB2 are not in the dictionary: distance), a struct_conn between their CB atoms
 * (no distance bonds between the two residues), and a residue of 300 carbons on a line 1.5 Å apart whose dictionary
 * bonds consecutive atoms, one pair at 3.0 Å and one at 4.5 Å (beyond the 4 Å radius)
 */
function syntheticCif() {
  const ala: [string, string, number, number, number][] = [
    ['N', 'N', 0, 0, 0],
    ['CA', 'C', 1.46, 0, 0],
    ['C', 'C', 2.0, 1.42, 0],
    ['O', 'O', 1.25, 2.39, 0],
    ['CB', 'C', 2.0, -0.77, 1.2],
    ['H', 'H', -0.5, -0.85, 0],
    ['HA', 'H', 1.8, -0.5, -0.9],
    ['HB1', 'H', 3.09, -0.77, 1.2],
    ['HB2', 'H', 1.64, -1.8, 1.2],
  ];
  const rows: string[] = [];
  let id = 1;
  for (const [seq, shift] of [
    [1, 0],
    [2, 3.0],
  ]) {
    for (const [name, el, ax, ay, az] of ala) {
      rows.push(
        `ATOM ${id++} ${el} ${name} . ALA A 1 ${seq} ? ${(ax + shift).toFixed(3)} ${ay.toFixed(3)} ${az.toFixed(3)} 1 0 ${seq} A 1`,
      );
    }
  }
  for (let i = 1; i <= 300; i++) {
    rows.push(`HETATM ${id++} C C${i} . BIG B 2 . ? ${(1.5 * i).toFixed(3)} 0.000 50.000 1 0 1 B 1`);
  }
  const bonds = [
    ['ALA', 'N', 'CA'],
    ['ALA', 'CA', 'C'],
    ['ALA', 'C', 'O'],
    ['ALA', 'CA', 'CB'],
    ['ALA', 'N', 'H'],
    ['ALA', 'CA', 'HA'],
    ...Array.from({ length: 299 }, (_, i) => ['BIG', `C${i + 1}`, `C${i + 2}`]),
    ['BIG', 'C10', 'C12'],
    ['BIG', 'C20', 'C23'],
  ];
  return `data_TEST
_entry.id TEST
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
loop_
_chem_comp.id
_chem_comp.type
ALA 'L-peptide linking'
BIG non-polymer
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
${bonds.map((b) => `${b.join(' ')} sing`).join('\n')}
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_symmetry
covale1 covale A ALA 1 CB 1_555 A ALA 2 CB 1_555
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
${rows.join('\n')}
#
`;
}

function signature(model: Model, props?: Partial<BondComputationProps>) {
  const rows: number[][] = [];
  for (const u of Structure.ofModel(model).units) {
    if (!Unit.isAtomic(u)) continue;
    const { a, b, edgeProps } = computeIntraUnitBonds(u, props);
    for (let i = 0; i < a.length; i++) {
      if (a[i] > b[i]) continue;
      rows.push([u.id, a[i], b[i], edgeProps.order[i], edgeProps.flags[i], edgeProps.key[i]]);
    }
  }
  rows.sort((x, y) => {
    for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) return x[i] - y[i];
    return 0;
  });
  return { count: rows.length, hash: (hashFnv32a(rows.flat()) >>> 0).toString(16) };
}

/** Moves atom `b` along the a–b direction to `distance` from `a`, in place */
function stretch(model: Model, a: number, b: number, distance: number) {
  const { x, y, z } = model.atomicConformation as unknown as { x: Float32Array; y: Float32Array; z: Float32Array };
  const dx = x[b] - x[a];
  const dy = y[b] - y[a];
  const dz = z[b] - z[a];
  const f = distance / Math.hypot(dx, dy, dz);
  x[b] = x[a] + dx * f;
  y[b] = y[a] + dy * f;
  z[b] = z[a] + dz * f;
}

/** Atom index of `atomId` in the residue of 1cbs with `compId` and `seqId` */
function atomOf(model: Model, compId: string, seqId: number, atomId: string) {
  const { atoms, residueAtomSegments, residues } = model.atomicHierarchy;
  for (let r = 0; r < residues._rowCount; r++) {
    if (residues.label_seq_id.value(r) !== seqId) continue;
    for (let i = residueAtomSegments.offsets[r]; i < residueAtomSegments.offsets[r + 1]; i++) {
      if (atoms.label_comp_id.value(i) === compId && atoms.label_atom_id.value(i) === atomId) return i;
    }
  }
  throw new Error(`no ${compId} ${seqId} ${atomId}`);
}

describe('intra-unit bonds', () => {
  it('from distances (1crn, no dictionary)', async () => {
    expect(signature(await parseCif(readExample('1crn.cif')))).toEqual({ count: 337, hash: '3a35ad0c' });
  });

  it('from the dictionary and distances (1cbs, 1grm)', async () => {
    expect(signature(await parseCif(readExample('1cbs_updated.cif')))).toEqual({ count: 1129, hash: 'c9e64fb5' });
    expect(signature(await parseCif(readExample('1grm_updated.cif')))).toEqual({ count: 286, hash: 'f7aea59c' });
  });

  it('from distances only (forceCompute)', async () => {
    const m = await parseCif(readExample('1cbs_updated.cif'));
    expect(signature(m, { forceCompute: true })).toEqual({ count: 1129, hash: 'f0d77f2a' });
  });

  it('coarse-grained (CA only)', async () => {
    expect(signature(await parseCif(caOnly(readExample('1crn.cif'))))).toEqual({ count: 45, hash: 'cf316a07' });
  });

  it('dictionary bonds up to the 4 Å radius, at any geometry', async () => {
    const m = await parseCif(readExample('1cbs_updated.cif'));
    // a stretched PHE ring bond (3.5 Å, found) and a stretched peptide side chain bond (4.5 Å, beyond the radius)
    stretch(m, atomOf(m, 'PHE', 3, 'CB'), atomOf(m, 'PHE', 3, 'CG'), 3.5);
    stretch(m, atomOf(m, 'LYS', 20, 'CD'), atomOf(m, 'LYS', 20, 'CE'), 4.5);
    expect(signature(m)).toEqual({ count: 1128, hash: '4f171223' });
    expect(signature(m, { maxRadius: 2 })).toEqual({ count: 1126, hash: '45f0735a' });
  });

  it('no bonds to atoms without coordinates (NaN), as on master', async () => {
    const m = await parseCif(readExample('1cbs_updated.cif'));
    const { x } = m.atomicConformation as unknown as { x: Float32Array };
    for (let i = 0; i < x.length; i += 41) x[i] = NaN;
    expect(signature(m)).toEqual({ count: 1074, hash: '8be27fcf' });
  });

  it('hydrogens off the dictionary, struct_conn between residues, a residue of 300 atoms (synthetic), as on master', async () => {
    const m = await parseCif(syntheticCif());
    const units = Structure.ofModel(m).units.map((u) => u.elements.length);
    expect(units).toEqual([18, 300]);
    expect(signature(m)).toEqual({ count: 317, hash: '774bf775' });
    expect(signature(m, { forceCompute: true })).toEqual({ count: 315, hash: '13902604' });
  });
});
