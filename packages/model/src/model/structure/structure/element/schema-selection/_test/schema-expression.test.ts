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
import { Column } from '@molstar/core/data/db';
import { Expression } from '@molstar/query-language/language/expression';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import { formatMolScript } from '@molstar/query-language/language/expression-formatter';
import { Script } from '@molstar/model/script/script';
import { compile } from '@molstar/model/script/runtime/query/compiler';
import { QueryContext, StructureSelection } from '@molstar/model/model/structure';
import { IndicesAndSortings } from '../indexing.js';
import { substructureFromRows } from '../row-structure.js';
import { atomQualifies, getAtomRangesForRows } from '../selections.js';
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

describe('Schema.toLoci, the generator element-schema', async () => {
  const structure = await createStructure();
  const nAtoms = structure.elementCount;
  const model = structure.models[0];

  /** The atoms (of the one unit) which satisfy the schema, by the per-atom check of the item (not through the indices) */
  function bruteForce(schema: Schema): number[] {
    const items: SchemaItem[] = [];
    // (a missing value, `undefined`, `null` or NaN, is not a criterion)
    Schema.forEachItem(schema, (item) =>
      items.push(
        Object.fromEntries(Object.entries(item).filter(([, v]) => v !== undefined && v !== null && !Number.isNaN(v))),
      ),
    );
    const atoms: number[] = [];
    for (let i = 0; i < nAtoms; i++) {
      if (items.some((item) => atomQualifies(model, i as any, item))) atoms.push(i);
    }
    return atoms;
  }
  const lociAtoms = (loci: Loci) =>
    Object.values(elementList(loci))
      .flat()
      .sort((a, b) => a - b);

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
    ['inverted range', { label_asym_id: 'A', beg_label_seq_id: 8, end_label_seq_id: 3 }],
    ['null bounds', { label_asym_id: 'A', beg_label_seq_id: 4, end_label_seq_id: null as any }],
    ['null begin', { auth_asym_id: 'b', beg_auth_seq_id: null as any, end_auth_seq_id: 103 }],
    ['null scalar items', { label_asym_id: 'C', label_seq_id: null as any, label_atom_id: undefined }],
    ['numbers as text select nothing', { label_asym_id: 'A', label_seq_id: '3' as any }],
    ['atom index as text', { atom_index: '4' as any }],
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
    [
      'object of arrays with holes',
      {
        items: {
          label_asym_id: ['A', 'B', 'C', 'A'],
          label_seq_id: [1, undefined, 3, 4],
          label_atom_id: ['CA', 'N', undefined, Number.NaN as any],
        },
      },
    ],
    ['empty items', { items: [] }],
  ];

  for (const [name, schema] of schemas) {
    it(`finds the elements of the schema: ${name}`, () => {
      expect(lociAtoms(Schema.toLoci(structure, schema))).toEqual(bruteForce(schema));
    });
  }

  it('is not vacuous: the schemas select different, non-empty parts', () => {
    const counts = schemas.map(([, s]) => StructureElement.Loci.size(Schema.toLoci(structure, s)));
    expect(counts.filter((n) => n > 0).length).toBeGreaterThan(15);
    expect(new Set(counts).size).toBeGreaterThan(8);
    expect(StructureElement.Loci.size(Schema.toLoci(structure, {}))).toBe(nAtoms);
  });

  it('the expression is a valid expression that survives JSON, in the shape of the schema', () => {
    for (const [, schema] of schemas) {
      const expression = Schema.toExpression(schema);
      expect(expressionValidationIssues(expression)).toBeUndefined();
      expect(Expression.is(expression)).toBe(true);
      expect(JSON.parse(JSON.stringify(expression))).toEqual(expression);
      expect((expression as any).head.name).toBe('structure-query.generator.element-schema');
    }
    // an item is the item with the values it has
    expect(Schema.toExpression({ label_asym_id: 'A', label_seq_id: undefined, auth_seq_id: null as any })).toEqual(
      MS.struct.generator.elementSchema([MS.struct.type.schemaItem({ label_asym_id: 'A' })]),
    );
    // a table of columns is one item, whose lists are the columns (and the prefix is the same in every row), a missing value is null
    expect(
      Schema.toExpression({
        prefix: { label_asym_id: 'A' },
        items: { label_seq_id: [1, undefined, 3], label_atom_id: ['CA', 'N', Number.NaN as any] },
      }),
    ).toEqual(
      MS.struct.generator.elementSchema([
        MS.struct.type.schemaItem({ label_asym_id: 'A', label_seq_id: [1, null, 3], label_atom_id: ['CA', 'N', null] }),
      ]),
    );
    expect(Schema.toExpression({ items: [] })).toEqual(MS.struct.generator.elementSchema([]));
  });

  it('many rows are given as columns, which the compiler does not look into', () => {
    const rows: SchemaItem[] = Array.from({ length: 50 }, (_, i) =>
      i % 3 === 0 ? { label_asym_id: 'A', label_seq_id: 1 + (i % 12) } : { label_asym_id: 'B', label_atom_id: 'CA' },
    );
    const expression = Schema.toExpression({ items: rows });
    const item = (expression as any).args[0];
    expect(item.head.name).toBe('structure-query.type.schema-item');
    expect(Object.keys(item.args).sort()).toEqual(['label_asym_id', 'label_atom_id', 'label_seq_id']);
    expect(item.args.label_asym_id).toHaveLength(50);
    expect(item.args.label_seq_id.filter((v: unknown) => v === null)).toHaveLength(33); // (the holes)
    expect(JSON.parse(JSON.stringify(expression))).toEqual(expression);
    expect(expressionValidationIssues(expression)).toBeUndefined();
    // (a list of values is a literal: compiling it is passing it through)
    expect(Expression.isLiteral(item.args.label_asym_id)).toBe(true);
    expect(lociAtoms(Loci.fromExpression(structure, expression))).toEqual(bruteForce({ items: rows }));
    // fewer rows are given as they are
    expect((Schema.toExpression({ items: rows.slice(0, 9) }) as any).args).toHaveLength(9);
  });

  it('a table of columns selects what the list of its rows does', () => {
    const columns: Schema = {
      prefix: { label_asym_id: 'B' },
      items: { label_seq_id: [1, 5, undefined, 12], label_atom_id: ['N', undefined, 'CA', 'O'] },
    };
    const rows: Schema = {
      items: [
        { label_asym_id: 'B', label_seq_id: 1, label_atom_id: 'N' },
        { label_asym_id: 'B', label_seq_id: 5 },
        { label_asym_id: 'B', label_atom_id: 'CA' },
        { label_asym_id: 'B', label_seq_id: 12, label_atom_id: 'O' },
      ],
    };
    const atoms = lociAtoms(Schema.toLoci(structure, columns));
    expect(atoms).toEqual(lociAtoms(Schema.toLoci(structure, rows)));
    // residue 5 (4 atoms), all CA of chain B (12 - 1 of residue 5 which is already there), N of 1, O of 12
    expect(atoms.length).toBe(4 + 11 + 1 + 1);
  });

  it('a table of columns: no columns or no values are no rows, columns of different lengths are rejected, the same everywhere', () => {
    const entries = (schema: Schema) => [
      StructureElement.Loci.size(Schema.toLoci(structure, schema)),
      Schema.toStructure(structure, schema).elementCount,
    ];
    for (const schema of [
      { items: {} },
      { prefix: { label_asym_id: 'A' }, items: {} },
      { items: { label_seq_id: [] } },
      { prefix: { label_asym_id: 'A' }, items: { label_seq_id: [], label_atom_id: [] } },
    ] as Schema[]) {
      expect(entries(schema), JSON.stringify(schema)).toEqual([0, 0]);
    }
    // (columns which are not fields count as columns: the rows are the same)
    expect(entries({ items: { group_id: ['x', 'y'] } } as any)).toEqual([nAtoms, nAtoms]);
    for (const schema of [
      { items: { label_asym_id: ['A', 'B'], label_seq_id: [1] } },
      { items: { group_id: ['x', 'y', 'z'], label_asym_id: ['A'] } },
    ] as Schema[]) {
      expect(() => Schema.toLoci(structure, schema)).toThrow('different lengths');
      expect(() => Schema.toStructure(structure, schema)).toThrow('different lengths');
    }
  });

  it('the elements are one set', () => {
    const run = (schema: Schema) =>
      (compile<StructureSelection>(Schema.toExpression(schema))(new QueryContext(structure)) as any)
        .structures as unknown[];
    expect(run({ label_asym_id: 'A' }).length).toBe(1);
    expect(
      run({ items: [{ label_asym_id: 'A' }, { label_asym_id: 'B' }, { label_asym_id: 'A', label_seq_id: 1 }] }).length,
    ).toBe(1);
    expect(run({ label_asym_id: 'Z' }).length).toBe(0);
    expect(run({ items: [] }).length).toBe(0);
  });

  it('is a generator like the others: the text language, modifiers and the structure it is run on', () => {
    const text =
      '(sel.atom.element-schema (schema-item :label_asym_id A :label_seq_id 3) (schema-item :label_asym_id C))';
    const parsed = Script.toExpression({ language: 'mol-script', expression: text });
    const schema: Schema = { items: [{ label_asym_id: 'A', label_seq_id: 3 }, { label_asym_id: 'C' }] };
    expect(lociAtoms(Loci.fromExpression(structure, parsed))).toEqual(lociAtoms(Schema.toLoci(structure, schema)));
    // (columns in the text language: lists, and `empty-value`)
    const columns = Script.toExpression({
      language: 'mol-script',
      expression: '(sel.atom.element-schema (schema-item :label_asym_id [A B] :label_seq_id [3 empty-value]))',
    });
    expect(lociAtoms(Loci.fromExpression(structure, columns))).toEqual(
      lociAtoms(Schema.toLoci(structure, { items: [{ label_asym_id: 'A', label_seq_id: 3 }, { label_asym_id: 'B' }] })),
    );
    // modifiers work on it
    const residues = Loci.fromExpression(
      structure,
      MS.struct.modifier.wholeResidues([Schema.toExpression({ label_asym_id: 'A', label_atom_id: 'CA' })]),
    );
    expect(lociAtoms(residues)).toEqual(bruteForce({ label_asym_id: 'A' }));
    // it is run on the structure which the query gets (a part of the structure has fewer elements)
    const half = structure.subsetBuilder(true);
    for (const unit of structure.units) {
      if (unit.elements.length > 0)
        half.setUnit(unit.id, unit.elements.slice(0, Math.ceil(unit.elements.length / 2)) as any);
    }
    const part = half.getStructure();
    const inPart = Schema.toLoci(part, { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 12 });
    expect(StructureElement.Loci.size(inPart)).toBeGreaterThan(0);
    expect(StructureElement.Loci.size(inPart)).toBeLessThan(
      StructureElement.Loci.size(
        Schema.toLoci(structure, { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 12 }),
      ),
    );
    expect(inPart.structure).toBe(part);
  });

  it('rejects what it cannot make sense of instead of selecting everything', () => {
    const run = (expression: Expression) => compile(expression)(new QueryContext(structure));
    const text = (expression: string) => run(Script.toExpression({ language: 'mol-script', expression }));
    expect(() => run(MS.struct.generator.elementSchema([MS.struct.type.schemaItem({ label_asym: 'A' })]))).toThrow(
      "'label_asym' is not a field of a schema item",
    );
    expect(() =>
      run(
        MS.struct.generator.elementSchema([
          MS.struct.type.schemaItem({ label_asym_id: ['A', 'B'], label_seq_id: [1] }),
        ]),
      ),
    ).toThrow('have different lengths');
    // an argument which is not an item (a list of them, a value) is not a row without criteria, which would be everything
    for (const expression of [
      '(sel.atom.element-schema [(schema-item :label_asym_id A)])',
      '(sel.atom.element-schema A)',
      '(sel.atom.element-schema 1)',
      '(sel.atom.element-schema [A B])',
    ]) {
      expect(() => text(expression), expression).toThrow('must be schema-item');
    }
    expect(
      expressionValidationIssues(MS.struct.generator.elementSchema([MS.struct.type.schemaItem({ label_asym: 'A' })])),
    ).toBeDefined();
    expect(expressionValidationIssues(MS.struct.generator.elementSchema([['A']] as any))).toBeDefined();
  });

  it('an item without fields is everything, no items are nothing', () => {
    const count = (expression: string) =>
      StructureElement.Loci.size(
        Loci.fromExpression(structure, Script.toExpression({ language: 'mol-script', expression })),
      );
    expect(count('(sel.atom.element-schema (schema-item))')).toBe(nAtoms);
    expect(count('(sel.atom.element-schema)')).toBe(0);
    expect(count('(sel.atom.element-schema (schema-item :label_asym_id A) (schema-item))')).toBe(nAtoms);
    expect(compile(MS.struct.type.schemaItem())(new QueryContext(structure))).toEqual([{}]);
  });

  it('is run in a query which can time out', () => {
    const schema: Schema = { label_asym_id: 'A' };
    const context = new QueryContext(structure, { timeoutMs: 1 });
    const expired = vi.spyOn(context, 'throwIfTimedOut').mockImplementation(() => {
      throw new Error('The query took too long to execute');
    });
    expect(() => compile(Schema.toExpression(schema))(context)).toThrow('too long');
    expect(expired).toHaveBeenCalled();
  });

  it('the text of an expression (also with columns and missing values) is the expression', () => {
    const rows: SchemaItem[] = Array.from({ length: 12 }, (_, i) =>
      i % 3 === 0
        ? { label_asym_id: 'A', label_seq_id: 1 + i }
        : i % 3 === 1
          ? { label_asym_id: 'B', label_atom_id: 'CA' }
          : { label_asym_id: 'C', auth_seq_id: 100 + i },
    );
    for (const schema of [{ items: rows }, { items: rows.slice(0, 4) }] as Schema[]) {
      const expression = Schema.toExpression(schema);
      const text = formatMolScript(expression);
      expect(text).toContain(
        Array.isArray((schema as any).items) && (schema as any).items.length > 9 ? 'empty-value' : 'schema-item',
      );
      const parsed = Script.toExpression({ language: 'mol-script', expression: text });
      expect(lociAtoms(Loci.fromExpression(structure, parsed))).toEqual(lociAtoms(Schema.toLoci(structure, schema)));
      expect(lociAtoms(Schema.toLoci(structure, schema)).length).toBeGreaterThan(0);
    }
    // a missing value is not the string null
    expect(formatMolScript(MS.struct.type.schemaItem({ label_asym_id: [null, 'A'] }))).toContain('[empty-value A]');
  });

  it('the loci is a valid loci of the structure (iterable, bundle, ordered by unit id)', () => {
    const ids = (loci: Loci) => loci.elements.map((e) => e.unit.id);
    for (const [, schema] of schemas) {
      const loci = Schema.toLoci(structure, schema);
      expect(loci.structure).toBe(structure);
      expect(Bundle.areEqual(Bundle.fromLoci(loci), Schema.toBundle(structure, schema))).toBe(true);
      expect(ids(loci)).toEqual([...ids(loci)].sort((a, b) => a - b));
    }
    const schema: Schema = { label_asym_id: 'A', beg_label_seq_id: 2, end_label_seq_id: 4 };
    const expected = Schema.toLoci(structure, schema);
    expect(Bundle.areEqual(Bundle.fromSchema(structure, schema), Bundle.fromLoci(expected))).toBe(true);
    expect(StructureElement.Loci.areEqual(Loci.fromSchema(structure, schema), expected)).toBe(true);
  });

  it('a residue without label_seq_id is not matched by a sequence criterion, a missing value is not 0', () => {
    expect(StructureElement.Loci.size(Schema.toLoci(structure, { label_seq_id: 0 }))).toBe(0);
    expect(StructureElement.Loci.size(Schema.toLoci(structure, { end_label_seq_id: 1 }))).toBe(3 * 4);
    const ligand = 3;
    expect(StructureElement.Loci.size(Schema.toLoci(structure, { label_asym_id: 'L' }))).toBe(ligand);
  });

  describe('values which are missing or of another type, at every entry', () => {
    const viaRows = (row: SchemaItem) => substructureFromRows(structure, [row]).elementCount;
    const viaRanges = (row: SchemaItem) => {
      const ranges = getAtomRangesForRows(
        [row],
        model,
        { instanceId: '1_555', name: '1_555' },
        IndicesAndSortings.get(model),
      );
      let n = 0;
      for (let i = 0; i < ranges.from.length; i++) n += ranges.to[i] - ranges.from[i];
      return n;
    };
    const entries = (row: SchemaItem) => [
      viaRows(row),
      viaRanges(row),
      Schema.toStructure(structure, row).elementCount,
      StructureElement.Loci.size(Schema.toLoci(structure, row)),
    ];
    const everywhere = (row: SchemaItem) => {
      const counts = entries(row);
      // (substructureFromRows, getAtomRangesForRows, Schema.toStructure, Schema.toLoci)
      expect(
        counts,
        JSON.stringify(row, (_, v) => (typeof v === 'number' && Number.isNaN(v) ? 'NaN' : v)),
      ).toEqual(counts.map(() => counts[0]));
      return counts[0];
    };

    it('a value which is of another type than its field selects nothing, whichever criterion finds the candidates', () => {
      for (const row of [
        { label_asym_id: 'A', end_auth_seq_id: '110' },
        { label_asym_id: 'A', end_auth_seq_id: '110', label_seq_id: 3 },
        { label_asym_id: 'A', beg_auth_seq_id: '101', residue_index: 2 },
        { label_asym_id: 'A', atom_index: 3, end_auth_seq_id: '110' },
        { label_asym_id: 'A', label_seq_id: '3' },
        { atom_index: '4' },
        { atom_id: '4' },
        { type_symbol: 6 },
        { label_asym_id: 1 },
        { label_seq_id: true },
        { label_asym_id: 'A', beg_auth_ins_code: 5, beg_auth_seq_id: 101 },
        { label_asym_id: 'A', pdbx_PDB_ins_code: 1 },
      ] as unknown as SchemaItem[]) {
        expect(everywhere(row), JSON.stringify(row)).toBe(0);
      }
      // (a list or an object cannot be a value in an expression, where a list is a column)
      expect(() => Schema.toExpression({ label_asym_id: 'A', label_comp_id: ['ALA'] as any })).toThrow(
        "the value of 'label_comp_id' must be a string or a number",
      );
      expect(() => Schema.toExpression({ items: { label_comp_id: [{}] as any } })).toThrow(
        'must be a string or a number',
      );
      // (the same rows with the types of the fields select something)
      expect(everywhere({ label_asym_id: 'A', end_auth_seq_id: 110 })).toBe(10 * 4);
      expect(everywhere({ label_asym_id: 'A', end_auth_seq_id: 110, label_seq_id: 3 })).toBe(4);
    });

    it('undefined, null and NaN are no criterion, not a value which nothing has', () => {
      const chainA = 12 * 4;
      for (const missing of [undefined, null, Number.NaN]) {
        const m = missing as any;
        expect(everywhere({ label_asym_id: 'A', label_seq_id: m })).toBe(chainA);
        expect(
          everywhere({ label_asym_id: 'A', auth_seq_id: m, beg_label_seq_id: m, residue_index: m, atom_id: m }),
        ).toBe(chainA);
        expect(everywhere({ label_asym_id: 'A', label_atom_id: m, type_symbol: m, label_comp_id: m })).toBe(chainA);
        expect(everywhere({ label_asym_id: m, auth_asym_id: m, label_entity_id: m, label_seq_id: 3 })).toBe(3 * 4);
        expect(everywhere({ label_asym_id: 'A', end_label_seq_id: 5, beg_label_seq_id: m })).toBe(5 * 4);
        // insertion codes of a range
        const range = { auth_asym_id: 'b', beg_auth_seq_id: 103, end_auth_seq_id: 106, end_auth_ins_code: '' };
        expect(everywhere({ ...range, beg_auth_ins_code: m })).toBe(everywhere(range));
        expect(everywhere({ ...range, beg_auth_ins_code: m })).toBe(4 * 4); // (103, 104, 105A, 106)
      }
    });
  });

  it('toStructure is the substructure of the same elements, the structure itself when all are selected', () => {
    for (const [, schema] of schemas) {
      const result = Schema.toStructure(structure, schema);
      expect(result.units.flatMap((u) => Array.from(u.elements))).toEqual(bruteForce(schema));
    }
    expect(Schema.toStructure(structure, {})).toBe(structure);
  });
});

void Column;
