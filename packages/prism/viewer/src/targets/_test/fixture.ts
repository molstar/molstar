/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '@molstar/core/data/db';
import { RuntimeContext } from '@molstar/core/task';
import { createModels } from '@molstar/model/formats/structure/basic/parser';
import { BasicSchema, createBasic } from '@molstar/model/formats/structure/basic/schema';
import { ComponentBuilder } from '@molstar/model/formats/structure/common/component';
import { EntityBuilder } from '@molstar/model/formats/structure/common/entity';
import { Structure } from '@molstar/model/model/structure';
import { MoleculeType } from '@molstar/model/model/structure/model/types';

// A structure built without a file: `@molstar/io`, which parses CIF, is not a dependency of this package, and the
// targets must not need the plugin. It has what the targets read: two protein chains (insertion code, a histidine, a
// cysteine, a threonine), a DNA chain, a heme in author chain A (its own label chain, as ligands are), waters in A, and
// a zinc ion in author chain B.
//
//   A   100-119 and 110A   polypeptide (label A, entity 1); residue i is ALA GLY SER THR CYS HIS by (i - 100) % 6
//                          (HIS 105, 111, 117; CYS 104, 110, 116; THR 103, 109, 115); 110A is GLY
//   A   500                HEM, label D, beside CYS 104
//   A   600-601            HOH, label F
//   B   200-209            polypeptide (label B, entity 2)
//   B   301                ZN2, label E
//   C   1-4                DNA (label C, entity 3)

interface Residue {
  auth: string;
  label: string;
  comp: string;
  seq: number;
  ins?: string;
  type: MoleculeType;
  /** Atom names and elements; the first atom (the backbone N, CA or the only atom) sits on the line of residues. */
  atoms: [string, string][];
  /** Offset along x of the residue (A). */
  x: number;
  y?: number;
}

const BACKBONE: [string, string][] = [
  ['N', 'N'],
  ['CA', 'C'],
  ['C', 'C'],
  ['O', 'O'],
];
const SIDE: Record<string, [string, string][]> = {
  ALA: [['CB', 'C']],
  GLY: [],
  SER: [
    ['CB', 'C'],
    ['OG', 'O'],
  ],
  THR: [
    ['CB', 'C'],
    ['OG1', 'O'],
    ['CG2', 'C'],
  ],
  CYS: [
    ['CB', 'C'],
    ['SG', 'S'],
  ],
  HIS: [
    ['CB', 'C'],
    ['CG', 'C'],
    ['ND1', 'N'],
    ['NE2', 'N'],
  ],
};

function protein(auth: string, seq: number, comp: string, x: number, ins?: string): Residue {
  return {
    auth,
    label: auth,
    comp,
    seq,
    ins,
    type: MoleculeType.Protein,
    atoms: [...BACKBONE, ...(SIDE[comp] ?? [])],
    x,
  };
}

function residues(): Residue[] {
  const out: Residue[] = [];
  const comps = ['ALA', 'GLY', 'SER', 'THR', 'CYS', 'HIS'];
  let x = 0;
  for (let seq = 100; seq < 120; seq++) {
    // 100 ALA, 101 GLY, 102 SER, 103 THR, 104 CYS, 105 HIS, 106 ALA, 107 GLY, ...
    out.push(protein('A', seq, comps[(seq - 100) % comps.length] as string, x));
    x += 3.8;
    if (seq === 110) {
      out.push(protein('A', 110, 'GLY', x, 'A'));
      x += 3.8;
    }
  }
  // The heme sits at the CYS 104 (x = 4 * 3.8), 5 A off the line.
  out.push({
    auth: 'A',
    label: 'D',
    comp: 'HEM',
    seq: 500,
    type: MoleculeType.Other,
    atoms: [
      ['FE', 'FE'],
      ['NA', 'N'],
      ['NB', 'N'],
    ],
    x: 4 * 3.8,
    y: 5,
  });
  out.push({
    auth: 'A',
    label: 'F',
    comp: 'HOH',
    seq: 600,
    type: MoleculeType.Water,
    atoms: [['O', 'O']],
    x: 60,
    y: 30,
  });
  out.push({
    auth: 'A',
    label: 'F',
    comp: 'HOH',
    seq: 601,
    type: MoleculeType.Water,
    atoms: [['O', 'O']],
    x: 64,
    y: 30,
  });
  x = 0;
  for (let seq = 200; seq < 210; seq++) {
    out.push(protein('B', seq, comps[(seq - 200) % comps.length] as string, x, undefined));
    out[out.length - 1]!.y = -40;
    x += 3.8;
  }
  out.push({
    auth: 'B',
    label: 'E',
    comp: 'ZN2',
    seq: 301,
    type: MoleculeType.Ion,
    atoms: [['ZN', 'ZN']],
    x: 20,
    y: -60,
  });
  x = 0;
  ['DA', 'DT', 'DG', 'DC'].forEach((comp, i) => {
    out.push({
      auth: 'C',
      label: 'C',
      comp,
      seq: i + 1,
      type: MoleculeType.DNA,
      atoms: [
        ['P', 'P'],
        ["O3'", 'O'],
        ["C3'", 'C'],
        ["O2'", 'O'],
        ["C4'", 'C'],
      ],
      x,
      y: 80,
    });
    x += 3.4;
  });
  return out;
}

export async function tinyStructure(): Promise<Structure> {
  const list = residues();
  const rows: {
    r: Residue;
    entity: string;
    name: string;
    element: string;
    labelSeq: number;
    xyz: [number, number, number];
  }[] = [];
  const entities = new EntityBuilder();
  const polymerSeq = new Map<string, number>();
  for (const r of list) {
    const entity = entities.getEntityId(r.comp, r.type, r.auth);
    const polymer = r.type === MoleculeType.Protein || r.type === MoleculeType.DNA;
    const labelSeq = polymer ? (polymerSeq.get(r.auth) ?? 0) + 1 : 0;
    if (polymer) polymerSeq.set(r.auth, labelSeq);
    r.atoms.forEach(([name, element], i) => {
      rows.push({
        r,
        entity,
        name,
        element,
        labelSeq,
        xyz: [
          r.x + (r.comp === 'HEM' ? 0.5 + i : (i % 4) * 0.4),
          (r.y ?? 0) + (r.comp === 'HEM' ? i : i * 0.3),
          i * 0.2,
        ],
      });
    });
  }
  const n = rows.length;
  const authSeq = Column.ofIntArray(rows.map((a) => a.r.seq));
  const authAtom = Column.ofStringArray(rows.map((a) => a.name));
  const authComp = Column.ofStringArray(rows.map((a) => a.r.comp));
  const atom_site = Table.ofPartialColumns(
    BasicSchema.atom_site,
    {
      id: Column.ofIntArray(rows.map((_, i) => i + 1)),
      auth_asym_id: Column.ofStringArray(rows.map((a) => a.r.auth)),
      auth_atom_id: authAtom,
      auth_comp_id: authComp,
      auth_seq_id: authSeq,
      pdbx_PDB_ins_code: Column.ofStringArray(rows.map((a) => a.r.ins ?? '')),
      Cartn_x: Column.ofFloatArray(Float32Array.from(rows.map((a) => a.xyz[0]))),
      Cartn_y: Column.ofFloatArray(Float32Array.from(rows.map((a) => a.xyz[1]))),
      Cartn_z: Column.ofFloatArray(Float32Array.from(rows.map((a) => a.xyz[2]))),
      label_asym_id: Column.ofStringArray(rows.map((a) => a.r.label)),
      label_atom_id: authAtom,
      label_comp_id: authComp,
      label_seq_id: Column.ofIntArray(rows.map((a) => a.labelSeq)),
      label_entity_id: Column.ofStringArray(rows.map((a) => a.entity)),
      occupancy: Column.ofConst(1, n, Column.Schema.float),
      type_symbol: Column.ofStringArray(rows.map((a) => a.element)),
      pdbx_PDB_model_num: Column.ofConst(1, n, Column.Schema.int),
    },
    n,
  );
  const components = new ComponentBuilder(authSeq, authAtom);
  let first = 0;
  for (const r of list) {
    if (!components.has(r.comp)) components.add(r.comp, first);
    first += r.atoms.length;
  }
  const basic = createBasic({
    entity: entities.getEntityTable(),
    chem_comp: components.getChemCompTable(),
    atom_site,
  });
  const trajectory = await createModels(
    basic,
    { kind: 'test', name: 'tiny', data: undefined },
    RuntimeContext.Synchronous,
  );
  return Structure.ofModel(trajectory.representative);
}
