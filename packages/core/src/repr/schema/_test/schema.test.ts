/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CategoryColumn, hashNode, type Node, pathToPointer } from '../../tree/index.js';
import {
  checkLanguage,
  defineLanguage,
  type Language,
  LanguageIssues,
  languageToMarkdown,
  normalize,
  SchemaIssues,
  signature,
  type UnitsHook,
  validate,
  validateAndNormalize,
  withDefaults,
} from '../index.js';
import { Mini } from './fixtures.js';

const spec = (...data: Node[]): Node => ({
  kind: 'spec',
  params: { version: '0.0.0' },
  children: [{ kind: 'data', id: 'data', children: data }],
});
const p53: Node = { kind: 'sequence', id: 'p53', params: { alphabet: 'protein', letters: 'MEEPQ' } };
const domains = (rows: Record<string, unknown>, extra: Record<string, unknown> = {}): Node => ({
  kind: 'ranges',
  id: 'domains',
  params: { target: '/data/p53', rows: rows as never, ...(extra as object) },
});
const codes = (root: Node, language: Language = Mini) => validate(language, root).map((i) => i.code);

describe('repr/schema: a language described as data', () => {
  it('checks a description, and throws on a broken one written in code', () => {
    expect(checkLanguage(Mini)).toEqual([]);
    const broken: Language = {
      id: 'b',
      version: '0',
      doc: '',
      root: 'nope',
      kinds: {
        Bad: { doc: '', parents: ['missing'] },
        'ext:ok': { doc: '', parents: [], params: { n: { type: 'integer', unit: 'A', doc: '' } } },
        t: {
          doc: '',
          parents: [],
          params: {
            c: { type: 'column', of: 'int32', missing: true, doc: '' },
            r: { type: 'table', key: 'k', columns: { x: { type: 'string', doc: '' } }, doc: '' },
            u: {
              type: 'union',
              of: [
                { type: 'record', fields: {} },
                { type: 'table', columns: {} },
              ],
              doc: '',
            },
            d: { type: 'number', min: 1, default: 0, doc: '' },
          },
        },
      },
    };
    const issues = checkLanguage(broken);
    expect(issues.map((i) => i.code)).toEqual([
      'repr.schema.language.root',
      'repr.schema.language.name',
      'repr.schema.language.parent',
      'repr.schema.language.unit',
      'repr.schema.language.column',
      'repr.schema.language.table',
      'repr.schema.language.table',
      'repr.schema.language.union',
      'repr.schema.language.default',
    ]);
    for (const i of issues) expect(LanguageIssues.codes).toHaveProperty([i.code]);
    expect(() => defineLanguage(broken)).toThrow(/repr.schema.language.root/);
  });

  it('asks a units hook whether a unit is one', () => {
    const units: UnitsHook = { isUnit: (u) => u === 'nm', quantity: () => undefined };
    expect(checkLanguage(Mini, { units }).map((i) => [i.code, i.address])).toEqual([
      ['repr.schema.language.unit', '/values/scale'],
    ]);
  });
});

describe('repr/schema: validation', () => {
  it('passes a valid spec, plain arrays in columns and typed ones alike', () => {
    expect(codes(spec(p53, domains({ start: [0, 10], end: [5, 20], id: ['a', 'b'], label: ['x', 'y'] })))).toEqual([]);
    const typed = domains({
      start: new Int32Array([0]),
      end: new Int32Array([5]),
      id: new CategoryColumn(['a'], new Int32Array([0])),
    });
    expect(codes(spec(p53, typed))).toEqual([]);
  });

  it('reports every misuse, by code, at its path and address, with its severity', () => {
    const bad: Node = {
      kind: 'spec',
      children: [
        {
          kind: 'data',
          id: 'data',
          children: [
            { kind: 'sequence', id: 'p53', params: { alphabet: 'klingon', extra: 1 } },
            { kind: 'sequence', id: 'p53', params: { alphabet: 'dna', cyclic: 'yes' } },
            { kind: 'sequence', params: { alphabet: 'dna' } },
            { kind: 'data', id: 'nested' },
            { kind: 'values', id: 'v', params: { target: 'p53', values: [1, 'x'], scale: '2nm' } },
            { kind: 'nope', id: 'n' },
            { kind: 'ext:track', id: 'e' },
            { kind: 'field', args: [] },
          ],
        },
      ],
    };
    const issues = validate(Mini, bad);
    expect(issues.map((i) => [i.code, pathToPointer(i.path ?? []), i.address])).toEqual([
      ['repr.schema.param.missing', '/params/version', '/'],
      ['repr.schema.param.unknown', '/children/0/children/0/params/extra', '/data/p53'],
      ['repr.schema.param.enum', '/children/0/children/0/params/alphabet', '/data/p53'],
      ['repr.schema.id.duplicate', '/children/0/children/1/id', '/data'],
      ['repr.schema.param.type', '/children/0/children/1/params/cyclic', '/data/p53'],
      ['repr.schema.id.missing', '/children/0/children/2', undefined],
      ['repr.schema.parent.invalid', '/children/0/children/3/kind', '/data/nested'],
      ['repr.schema.ref.address', '/children/0/children/4/params/target', '/data/v'],
      ['repr.schema.column.type', '/children/0/children/4/params/values/1', '/data/v'],
      ['repr.schema.param.type', '/children/0/children/4/params/scale', '/data/v'],
      ['repr.schema.kind.unknown', '/children/0/children/5/kind', '/data/n'],
      ['repr.schema.kind.extension', '/children/0/children/6/kind', '/data/e'],
      ['repr.schema.parent.invalid', '/children/0/children/7/kind', undefined],
      ['repr.schema.args.count', '/children/0/children/7/args', undefined],
    ]);
    for (const i of issues) {
      expect(SchemaIssues.codes).toHaveProperty([i.code]);
      expect(i.severity).toBe(SchemaIssues.codes[i.code as keyof typeof SchemaIssues.codes].severity);
    }
    expect(issues.find((i) => i.code === 'repr.schema.kind.extension')?.severity).toBe('warning');
  });

  it('checks columns and tables: element types, lengths, extra columns, the key', () => {
    const cases: [Record<string, unknown>, Record<string, unknown>, string[]][] = [
      [{ start: [0], end: [1, 2] }, {}, ['repr.schema.table.length']],
      [{ start: [0], end: [1], note: 3 }, {}, ['repr.schema.table.column']],
      [{ start: [0], end: [1], note: ['kept'] }, {}, []],
      [{ start: [0, 1], end: [1, 2], id: ['a', 'a'] }, {}, ['repr.schema.table.key']],
      [{ start: new Float32Array([0]), end: [1] }, {}, ['repr.schema.column.type']],
      [{ start: [0.5], end: [1] }, {}, ['repr.schema.column.type']],
      [{ start: [2 ** 31], end: [1] }, {}, ['repr.schema.column.type']],
      [{ start: [0], end: [1], id: [null] }, {}, ['repr.schema.column.type']],
      [{ start: [0] }, {}, ['repr.schema.param.missing']],
      [{ start: [0], end: [1] }, { color: { kind: 'sequence', id: 'x' } }, ['repr.schema.slot.kind']],
      [{ start: [0], end: [1] }, { color: { kind: 'field', args: [1] } }, ['repr.schema.param.type']],
      [{ start: [0], end: [1] }, { color: { kind: 'field', args: ['label'] } }, []],
    ];
    for (const [rows, extra, want] of cases)
      expect(codes(spec(p53, domains(rows, extra))), JSON.stringify(rows)).toEqual(want);
  });

  it('reads a float column with null and NaN as missing', () => {
    const v: Node = { kind: 'values', id: 'v', params: { target: '/data/p53', values: [0.5, null, Number.NaN] } };
    expect(codes(spec(p53, v))).toEqual([]);
  });

  it('reads quantities through the units hook only', () => {
    const units: UnitsHook = {
      isUnit: () => true,
      quantity: (text, unit) => {
        const m = /^([0-9.]+)nm$/.exec(text);
        if (!m) return undefined;
        return unit === 'A' ? { value: Number(m[1]) * 10 } : { message: 'This number takes no unit' };
      },
    };
    const v = (scale: unknown): Node => ({
      kind: 'values',
      id: 'v',
      params: { target: '/data/p53', scale: scale as never },
    });
    expect(validate(Mini, spec(p53, v('2nm')), { units })).toEqual([]);
    expect(validate(Mini, spec(p53, v('2parsec')), { units }).map((i) => i.code)).toEqual(['repr.schema.param.type']);
    expect(validate(Mini, spec(p53, v('2nm'))).map((i) => i.code)).toEqual(['repr.schema.param.type']);
  });

  it('stops at maxIssues with a warning', () => {
    const many = spec(...Array.from({ length: 5 }, (_, i): Node => ({ kind: 'sequence', id: `s${i}` })));
    const issues = validate(Mini, many, { maxIssues: 3 });
    expect(issues).toHaveLength(4);
    expect(issues[3]).toMatchObject({ code: 'repr.schema.limit_exceeded', severity: 'warning' });
  });
});

describe('repr/schema: normalization on read', () => {
  const plain = spec(
    p53,
    domains(
      { start: [0, 10], end: [5, 20], id: ['a', 'b'], label: ['x', 'y'], note: [1, 2] },
      { color: { kind: 'field', args: ['id'] } },
    ),
    { kind: 'values', id: 'v', params: { target: '/data/p53', values: [0.5, null] } },
  );

  it('types every column param as declared, and keeps what is not one', () => {
    const n = normalize(Mini, plain);
    const rows = n.children![0].children![1].params!.rows as Record<string, unknown>;
    expect(rows.start).toEqual(new Int32Array([0, 10]));
    expect(rows.end).toBeInstanceOf(Int32Array);
    expect(rows.id).toBeInstanceOf(CategoryColumn);
    expect(rows.label).toEqual(['x', 'y']);
    expect(rows.note).toEqual([1, 2]);
    const values = n.children![0].children![2].params!.values as Float32Array;
    expect(values).toBeInstanceOf(Float32Array);
    expect(Number.isNaN(values[1])).toBe(true);
    // the unchanged subtree is shared, and a normalized tree comes back as itself
    expect(n.children![0].children![0]).toBe(p53);
    expect(normalize(Mini, n)).toBe(n);
    expect(hashNode(normalize(Mini, plain))).toBe(hashNode(n));
    expect(hashNode(n)).not.toBe(hashNode(plain));
  });

  it('validates, then normalizes when there is no error', () => {
    expect(validateAndNormalize(Mini, plain).node).toEqual(normalize(Mini, plain));
    const bad = validateAndNormalize(Mini, spec(p53, domains({ start: [0] })));
    expect(bad.node).toBeUndefined();
    expect(bad.issues.map((i) => i.code)).toEqual(['repr.schema.param.missing']);
  });
});

describe('repr/schema: generated from a description', () => {
  it('fills defaults, and gives signatures and reference docs', () => {
    const filled = withDefaults(Mini, spec(p53));
    expect(filled.children![0].children![0].params).toEqual({ alphabet: 'protein', letters: 'MEEPQ', cyclic: false });
    expect(signature(Mini, 'sequence')).toBe(
      '(sequence :alphabet "protein" | "dna" | "rna" [:letters string] [:cyclic boolean])',
    );
    const md = languageToMarkdown(Mini);
    expect(md).toContain('## `ranges`');
    expect(md).toContain('| `rows.start` | int32 column | required | Starts, 0-based. |');
  });
});
