/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '../../../../structure.js';
import { StructureElement } from '../../../element.js';
import { Bundle } from '../../bundle.js';
import { Loci } from '../../loci.js';
import { Schema, type SchemaItem } from '../../schema.js';
import { IndicesAndSortings } from '../indexing.js';
import { unitIndicesInRanges } from '../to-loci.js';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import '@molstar/model/script/runtime/query/table'; // (registers the symbols of the queries)

/** Chains A, B (entity 1) and C (entity 2) of 12 residues with 4 atoms; every 5th residue has the insertion code A.
 * Chain L is a ligand with 3 atoms and no label_seq_id. */
function mmcif(): string {
  const lines = [
    'data_T',
    'loop_',
    ...'group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'
      .split(' ')
      .map((c) => `_atom_site.${c}`),
  ];
  const names = ['N', 'CA', 'C', 'O'];
  let id = 1;
  for (const [asym, entity] of [
    ['A', '1'],
    ['B', '1'],
    ['C', '2'],
  ]) {
    for (let r = 1; r <= 12; r++) {
      const comp = r % 2 ? 'ALA' : 'GLY';
      names.forEach((name, a) => {
        lines.push(
          `ATOM ${id++} ${name[0]} ${name} . ${comp} ${asym} ${entity} ${r} ${r % 5 === 0 ? 'A' : '?'} ${a} ${r} 0 1 10 ${r + 100} ${comp} ${asym.toLowerCase()} ${name} 1`,
        );
      });
    }
  }
  ['X1', 'X2', 'X3'].forEach((name, a) => {
    lines.push(`HETATM ${id++} C ${name} . LIG L 3 . ? ${a} 50 0 1 10 900 LIG L ${name} 1`);
  });
  return `${lines.join('\n')}\n`;
}

async function createStructure() {
  const parsed = await CIF.parseText(mmcif()).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

function elementList(loci: Loci) {
  const out: Record<number, number[]> = {};
  StructureElement.Loci.forEachLocation(loci, (loc) => {
    (out[loc.unit.id] ??= []).push(loc.element);
  });
  return out;
}

describe('Schema.toLoci with and without the index', async () => {
  const structure = await createStructure();
  const nAtoms = structure.elementCount;
  const viaQuery = (schema: Schema) => Schema.toLoci(structure, schema);
  const viaIndex = (schema: Schema) => Schema.toLoci(structure, schema, { useIndex: true });

  const schemas: [string, Schema][] = [
    ['chain', { label_asym_id: 'B' }],
    ['entity', { label_entity_id: '1' }],
    ['residue', { label_asym_id: 'A', label_seq_id: 7 }],
    ['atom', { label_asym_id: 'C', label_seq_id: 3, label_atom_id: 'CA' }],
    ['auth residue with insertion code', { auth_asym_id: 'a', auth_seq_id: 105, pdbx_PDB_ins_code: 'A' }],
    ['no insertion code', { label_asym_id: 'A', pdbx_PDB_ins_code: '' }],
    ['range', { label_asym_id: 'B', beg_label_seq_id: 3, end_label_seq_id: 8 }],
    ['auth range', { auth_asym_id: 'c', beg_auth_seq_id: 103, end_auth_seq_id: 106 }],
    ['comp id', { label_comp_id: 'GLY', label_asym_id: 'A' }],
    ['element', { type_symbol: 'o', label_asym_id: 'C' }],
    ['atom id', { atom_id: 5 }],
    ['atom index', { atom_index: 6 }],
    ['residue index', { residue_index: 4 }],
    ['ligand', { label_comp_id: 'LIG' }],
    ['nothing', { label_asym_id: 'Z' }],
    ['null bounds', { label_asym_id: 'A', beg_label_seq_id: 4, end_label_seq_id: null as any }],
    ['null begin', { auth_asym_id: 'b', beg_auth_seq_id: null as any, end_auth_seq_id: 103 }],
    ['null scalar items', { label_asym_id: 'C', label_seq_id: null as any, label_atom_id: undefined }],
    ['numbers as text select nothing in both', { label_asym_id: 'A', label_seq_id: '3' as any }],
    ['atom index as text', { atom_index: '4' as any }],
    ['auth residue as text', { auth_asym_id: 'a', auth_seq_id: '105' as any }],
    ['everything', {}],
    [
      'array of items',
      {
        items: [
          { label_asym_id: 'A', label_seq_id: 1 },
          { label_asym_id: 'C', label_seq_id: 12 },
          { label_asym_id: 'A', label_seq_id: 1 },
        ],
      },
    ],
    [
      'prefix and items',
      { prefix: { label_asym_id: 'B' }, items: [{ label_seq_id: 2 }, { label_seq_id: 9, label_atom_id: 'O' }] },
    ],
    [
      'prefix and object of arrays',
      { prefix: { label_asym_id: 'C' }, items: { label_seq_id: [1, 5, 12], label_atom_id: ['N', 'CA', 'O'] } },
    ],
    ['empty items', { items: [] }],
  ];

  for (const [name, schema] of schemas) {
    it(`gives the same loci: ${name}`, () => {
      expect(elementList(viaIndex(schema))).toEqual(elementList(viaQuery(schema)));
    });
  }

  it('is not vacuous: the schemas select different, non-empty parts', () => {
    const counts = schemas.map(([, s]) => StructureElement.Loci.size(viaQuery(s)));
    expect(counts.filter((n) => n > 0).length).toBeGreaterThan(15);
    expect(new Set(counts).size).toBeGreaterThan(8);
    expect(StructureElement.Loci.size(viaQuery({}))).toBe(nAtoms);
  });

  it('the loci of the index is a valid loci of the structure (iterable, bundle, equal to the one of the query)', () => {
    for (const [, schema] of schemas) {
      const a = viaIndex(schema);
      const b = viaQuery(schema);
      expect(a.structure).toBe(structure);
      expect(StructureElement.Loci.size(a)).toBe(StructureElement.Loci.size(b));
      expect(StructureElement.Loci.areEqual(a, b) || StructureElement.Loci.isEmpty(a)).toBe(true);
      expect(Bundle.areEqual(Bundle.fromLoci(a), Bundle.fromLoci(b))).toBe(true);
    }
  });

  it('the elements of the loci are ordered by unit id (as are the units of every structure), like those of a query', () => {
    const schema: Schema = { label_asym_id: 'A' };
    const ids = (loci: Loci) => loci.elements.map((e) => e.unit.id);
    const result = ids(viaIndex({}));
    expect(result).toEqual([...result].sort((a, b) => a - b));
    expect(result.length).toBeGreaterThan(1);
    const reversed = Structure.create([...structure.units].reverse(), { parent: structure.parent });
    const loci = Schema.toLoci(reversed, { label_entity_id: '1' }, { useIndex: true });
    expect(ids(loci)).toEqual(ids(Schema.toLoci(reversed, { label_entity_id: '1' })));
    expect(ids(loci)).toEqual([...ids(loci)].sort((a, b) => a - b));
    expect(ids(Schema.toLoci(structure, schema, { useIndex: true }))).toEqual(ids(viaQuery(schema)));
  });

  it('toBundle, Loci.fromSchema and Bundle.fromSchema take the option too', () => {
    const schema: Schema = { label_asym_id: 'A', beg_label_seq_id: 2, end_label_seq_id: 4 };
    const expected = Bundle.fromLoci(viaQuery(schema));
    expect(Bundle.areEqual(Schema.toBundle(structure, schema, { useIndex: true }), expected)).toBe(true);
    expect(Bundle.areEqual(Bundle.fromSchema(structure, schema, { useIndex: true }), expected)).toBe(true);
    expect(
      StructureElement.Loci.areEqual(Loci.fromSchema(structure, schema, { useIndex: true }), viaQuery(schema)),
    ).toBe(true);
  });

  it('the indices are used only when asked for, and built once', () => {
    const create = vi.spyOn(IndicesAndSortings, 'create');
    try {
      viaQuery({ label_asym_id: 'A' });
      Schema.toLoci(structure, { label_asym_id: 'A' }, { useIndex: false });
      expect(create).not.toHaveBeenCalled(); // (the default is a query)
      viaIndex({ label_asym_id: 'A' });
      viaIndex({ label_asym_id: 'B' });
      expect(create.mock.calls.length).toBeLessThanOrEqual(1); // (once for the model, or earlier by an earlier test)
    } finally {
      create.mockRestore();
    }
  });

  it('the loci of a substructure only has the elements of the substructure', () => {
    const half = structure.subsetBuilder(true);
    for (const unit of structure.units) {
      if (unit.elements.length > 0)
        half.setUnit(unit.id, unit.elements.slice(0, Math.ceil(unit.elements.length / 2)) as any);
    }
    const part = half.getStructure();
    const schema: Schema = { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 12 };
    expect(elementList(Schema.toLoci(part, schema, { useIndex: true }))).toEqual(
      elementList(Schema.toLoci(part, schema)),
    );
    expect(StructureElement.Loci.size(Schema.toLoci(part, schema, { useIndex: true }))).toBeLessThan(
      StructureElement.Loci.size(viaQuery(schema)),
    );
  });

  it('where the index differs from the query: residues without label_seq_id are not matched by a sequence criterion', () => {
    const ligand = 3;
    expect(StructureElement.Loci.size(viaQuery({ label_seq_id: 0 }))).toBe(ligand);
    expect(StructureElement.Loci.size(viaIndex({ label_seq_id: 0 }))).toBe(0);
    expect(StructureElement.Loci.size(viaQuery({ end_label_seq_id: 1 }))).toBe(3 * 4 + ligand);
    expect(StructureElement.Loci.size(viaIndex({ end_label_seq_id: 1 }))).toBe(3 * 4);
  });
});

describe('unitIndicesInRanges', () => {
  const elements = SortedArray.ofSortedArray([2, 3, 4, 10, 11, 20, 30, 31, 32]);
  const make = (ranges: [number, number][]) => ({ from: ranges.map((r) => r[0]), to: ranges.map((r) => r[1]) }) as any;
  const positions = (ranges: [number, number][]) => {
    const set = unitIndicesInRanges(elements as any, make(ranges));
    if (!set) return [];
    return Array.from({ length: OrderedSet.size(set) }, (_, i) => OrderedSet.getAt(set, i));
  };
  const brute = (ranges: [number, number][]) =>
    Array.from(elements).flatMap((e, i) => (ranges.some(([a, b]) => e >= a && e < b) ? [i] : []));

  it('finds the positions of the elements in the ranges, for few ranges and for many', () => {
    for (const ranges of [
      [],
      [[0, 100]],
      [
        [2, 4],
        [10, 12],
      ], // (positions 0-1 and 3-4: one element between them)
      [
        [2, 5],
        [10, 12],
      ], // (positions 0-2 and 3-4: next to each other)
      [
        [2, 4],
        [4, 5],
        [10, 12],
      ], // (touching ranges)
      [[3, 5]],
      [
        [0, 3],
        [4, 12],
        [30, 31],
      ],
      [[5, 9]],
      [
        [0, 1],
        [2, 3],
        [3, 4],
        [10, 12],
        [20, 21],
        [31, 33],
        [40, 50],
        [60, 70],
        [80, 90],
        [100, 110],
      ],
    ] as [number, number][][]) {
      expect(positions(ranges)).toEqual(brute(ranges));
    }
  });
});
