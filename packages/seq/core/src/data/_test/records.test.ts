/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { encodeMsgPack } from '@molstar/core/data/msgpack';
import { MsgPackCodec } from '@molstar/core/repr/codecs/msgpack';
import { CategoryColumn, type Column, categoryColumn, type Issue, type Value } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { FastaRuntime } from '../fasta.js';
import {
  DefaultOpenLimits,
  type FormatRuntime,
  type OpenContext,
  type RangesList,
  type RegionAnswer,
  type SourceEntry,
  type ValuesPositions,
} from '../protocol.js';
import {
  RecordsIssues,
  type RecordsInput,
  RecordsRuntime,
  type RecordsSource,
  type RecordsTransform,
  readRecords,
} from '../records.js';
import { checkCoverage, seqToValue } from '../seq-format.js';
import { openSeq, SeqRuntime, SourceOpenError } from '../seq-runtime.js';
import { checkAnswer } from '../tables.js';

// The records runtime (_spec/README.md "Data", Records): the formats' conformance suite, one table in every encoding,
// the transform's property, typed columns kept by identity, each issue, fuzz without a throw, bounds before parsing.

const signal = new AbortController().signal;
const limits = DefaultOpenLimits;

function seeded(seed: number) {
  let s = seed;
  return () => {
    s = (s * 1103515245 + 12345) & 0x7fffffff;
    return s / 0x7fffffff;
  };
}

function read(input: RecordsInput, transform: RecordsTransform = {}, options = { system: 's', limits }) {
  return readRecords(input, transform, options);
}

/** The ranges table of a system in read data. */
function ranges(r: ReturnType<typeof read>, system = 's'): RangesList {
  expect(r.data, JSON.stringify(r.issues)).toBeDefined();
  const pair = r.data!.pairs.find((p) => p.request.system === system)!;
  return (pair.answer as Extract<RegionAnswer, { ok: true }>).tables.ranges as RangesList;
}

function values(r: ReturnType<typeof read>, system = 's'): ValuesPositions {
  expect(r.data, JSON.stringify(r.issues)).toBeDefined();
  const pair = r.data!.pairs.find((p) => p.request.system === system)!;
  return (pair.answer as Extract<RegionAnswer, { ok: true }>).tables.values as ValuesPositions;
}

const plain = (c: Column | undefined) =>
  c === undefined
    ? undefined
    : c instanceof CategoryColumn
      ? Array.from({ length: c.length }, (_, i) => c.get(i))
      : [...c];

/** A ranges table as plain arrays, its fields by name. */
function summary(t: RangesList) {
  return {
    start: [...t.start],
    end: [...t.end],
    strand: plain(t.strand),
    id: plain(t.id),
    row: plain(t.row),
    fields: Object.fromEntries(
      (t.fields ?? []).map((f) => [f.name, { type: f.column.constructor.name, v: plain(f.column) }]),
    ),
  };
}

const codes = (issues: readonly Issue[]) => issues.map((i) => i.code);

function ctx(files: Record<string, string | Uint8Array> = {}, l = limits): OpenContext {
  const bytes = (uri: string) => {
    const f = files[uri];
    return typeof f === 'string' ? new TextEncoder().encode(f) : f;
  };
  return {
    signal,
    limits: l,
    assets: {
      readText: async (uri) => (uri in files ? new TextDecoder().decode(bytes(uri)) : Promise.reject(new Error('no'))),
      readBytes: async (uri) => (uri in files ? bytes(uri) : Promise.reject(new Error('no'))),
    },
  };
}

// ---- a random table and its encodings ----

interface Row {
  start: number;
  end: number;
  strand: string;
  id: string;
  label: string;
  score: number;
  value: number;
}

const Words = ['kinase', 'a, b', 'say "hi"', ' padded ', 'tab\there', 'line\nbreak', 'é☃', '', '#x', '01'];

function randomTable(rand: () => number, n = 1 + Math.floor(rand() * 20)): Row[] {
  return Array.from({ length: n }, (_, k) => {
    const start = 3 + Math.floor(rand() * 200);
    // one row in seven ends before it starts in every convention: left out alike
    const dropped = rand() < 1 / 7;
    return {
      start,
      end: dropped ? start - 2 : start + Math.floor(rand() * 30),
      strand: ['+', '-', '.'][Math.floor(rand() * 3)],
      id: `f${k}`,
      label: Words[Math.floor(rand() * Words.length)],
      score: Math.floor(rand() * 64) / 4,
      value: Math.floor(rand() * 64) / 8,
    };
  });
}

const Columns = ['start', 'end', 'strand', 'id', 'label', 'score', 'value'] as const;

function delimited(rows: readonly Row[], d: string): string {
  const cell = (v: string | number) => {
    const s = String(v);
    return s.includes(d) || /["\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [Columns.join(d), ...rows.map((r) => Columns.map((c) => cell(r[c])).join(d))].join(
    rand2() < 0.5 ? '\n' : '\r\n',
  );
}
const rand2 = seeded(99);

const soa = (rows: readonly Row[]) => Object.fromEntries(Columns.map((c) => [c, rows.map((r) => r[c])]));

/** Every encoding of one table, each with the transform that reads it as the TSV is read. */
function encodings(rows: readonly Row[]): { name: string; input: RecordsInput; transform: RecordsTransform }[] {
  const typed: RecordsTransform = { types: { score: 'float64' } };
  const zero = { ...typed, base: 0, ends: 'exclusive' } as const;
  const node = {
    kind: 'table',
    params: {
      start: Int32Array.from(rows, (r) => r.start - 1),
      end: Int32Array.from(rows, (r) => r.end),
      strand: categoryColumn(rows.map((r) => r.strand)),
      id: categoryColumn(rows.map((r) => r.id)),
      label: categoryColumn(rows.map((r) => r.label)),
      score: Float64Array.from(rows, (r) => r.score),
      value: Float32Array.from(rows, (r) => r.value),
    },
  };
  return [
    { name: 'tsv', input: { format: 'tsv', text: delimited(rows, '\t') }, transform: typed },
    { name: 'csv', input: { format: 'csv', text: delimited(rows, ',') }, transform: typed },
    { name: 'json aos', input: { format: 'json', text: JSON.stringify(rows) }, transform: typed },
    {
      name: 'json soa',
      input: { format: 'json', text: JSON.stringify({ features: soa(rows) }) },
      transform: { ...typed, records: 'features' },
    },
    { name: 'msgpack aos', input: { format: 'msgpack', bytes: encodeMsgPack(rows) }, transform: typed },
    { name: 'msgpack soa', input: { format: 'msgpack', bytes: encodeMsgPack(soa(rows)) }, transform: typed },
    {
      name: 'msgpack typed soa',
      input: { format: 'msgpack', bytes: MsgPackCodec.encode(node) },
      transform: { ...zero, records: 'params' },
    },
    {
      name: 'inline aos',
      input: { value: rows.map((r) => ({ ...r, start: r.start - 1 })) as unknown as Value },
      transform: typed,
    },
    { name: 'inline typed soa', input: { value: node.params as unknown as Value }, transform: typed },
  ];
}

// ---- the conformance suite: every runtime's contract ----

/** The format readers' shared contract: a coverage that checks, answers that check, empty and edge regions, cancellation. */
async function conformance(runtime: FormatRuntime, entry: SourceEntry, context: OpenContext) {
  const s = await runtime.open(entry, context);
  expect(checkCoverage(s.coverage)).toBeUndefined();
  for (const { system, shape } of s.coverage.shapes) {
    const length = s.coverage.systems.find((x) => x.name === system)?.length ?? 1024;
    for (const [start, end] of [
      [0, length],
      [0, 0],
      [Math.max(0, length - 1), length],
      [length, length + 7],
    ]) {
      const request = { system, start, end, level: 1, shapes: [shape] };
      const answer = await s.read(request, signal);
      expect(checkAnswer(answer, request), `${system} ${shape} ${start}-${end}`).toBeUndefined();
    }
  }
  const unknown = await s.read({ system: '\u0000none', start: 0, end: 1, level: 1, shapes: ['ranges'] }, signal);
  expect(unknown).toMatchObject({ ok: false, kind: 'not-served' });
  const stop = new AbortController();
  stop.abort();
  const first = s.coverage.shapes[0];
  await expect(
    s.read({ system: first.system, start: 0, end: 1, level: 1, shapes: [first.shape] }, stop.signal),
  ).rejects.toMatchObject({ name: 'AbortError' });
  s.dispose();
}

describe('the formats’ conformance suite', () => {
  const rows = randomTable(seeded(3), 12);
  for (const e of encodings(rows)) {
    it(`records: ${e.name}`, async () => {
      const files: Record<string, string | Uint8Array> = {};
      let entry: SourceEntry;
      const format = 'format' in e.input ? e.input.format : 'json';
      if ('value' in e.input) entry = { address: '/data/t', format, value: e.input.value };
      else {
        files.f = 'text' in e.input ? e.input.text : e.input.bytes;
        entry = { address: '/data/t', format, uri: 'f' };
      }
      await conformance(RecordsRuntime, { ...entry, ...e.transform } as SourceEntry, ctx(files));
    });
  }

  it('records with letters, values and other systems', async () => {
    const doc = {
      sequence: { value: 'MEEPQSDPSV' },
      features: [
        { begin: '2', end: '4', type: 'Domain' },
        { begin: 9, end: 12, type: 'Region' },
        { begin: 3, end: 3, type: 'Site', target: 'P1-2' },
      ],
    };
    const entry = {
      address: '/data/p',
      format: 'json',
      uri: 'p.json',
      records: 'features',
      letters: 'sequence.value',
      fields: { start: 'begin', category: 'type' },
    } as SourceEntry;
    await conformance(RecordsRuntime, entry, ctx({ 'p.json': JSON.stringify(doc) }));
    const v = { address: '/data/v', format: 'tsv', uri: 'v.tsv' };
    await conformance(RecordsRuntime, v, ctx({ 'v.tsv': 'position\tvalue\n1\t0.5\n4\t1\n' }));
  });

  it('the seq and FASTA runtimes', async () => {
    const data = readRecords({ format: 'tsv', text: 'start\tend\n1\t5\n3\t9\n' }, {}, { system: 's', limits }).data!;
    await conformance(SeqRuntime, { address: '/data/q', format: 'seq', value: seqToValue(data) }, ctx());
    await conformance(
      FastaRuntime,
      { address: '/data/f', format: 'fasta', uri: 'f.fa' },
      ctx({ 'f.fa': '>a\nACGT\nAC\n>b\nGG\n' }),
    );
  });
});

// ---- one table in every encoding ----

describe('one table in every encoding', () => {
  it('reads equal tables from TSV, CSV, JSON and MessagePack in both layouts, and inline (a property, seeded)', () => {
    const rand = seeded(11);
    for (let t = 0; t < 30; t++) {
      const rows = randomTable(rand);
      const all = encodings(rows).map((e) => ({ name: e.name, r: read(e.input, e.transform) }));
      // the rows left out are said alike: one cause, one count
      const dropped = rows.filter((r) => r.end < r.start).length;
      for (const { name, r } of all) {
        expect(codes(r.issues), name).toEqual(dropped ? ['seq.source.records.order'] : []);
        if (dropped) expect(r.issues[0].message, name).toMatch(new RegExp(`: ${dropped} rows? \\(`));
      }
      const expected = summary(ranges(all[0].r));
      expect(expected.start.length).toBe(rows.length - dropped);
      expect(expected.fields.value.type).toBe('Float32Array');
      // rows sorted by start, each keeping its place in the file
      for (let k = 1; k < expected.start.length; k++)
        expect(expected.start[k]).toBeGreaterThanOrEqual(expected.start[k - 1]);
      for (const { name, r } of all.slice(1)) expect(summary(ranges(r)), name).toEqual(expected);
    }
  });

  it('keeps text cells as written: no trimming, `01` stays `01`, numbers parse past spaces', () => {
    const r = read({ format: 'csv', text: 'start,end,label,code\n 3 ,5, x ,01\n' });
    const t = ranges(r);
    expect([t.start[0], t.end[0]]).toEqual([2, 5]);
    expect(summary(t).fields).toMatchObject({ label: { v: [' x '] }, code: { v: ['01'] } });
  });
});

// ---- the transform ----

describe('the transform', () => {
  it('reads a table written 1/inclusive and 0/exclusive as the same rows (a property, seeded)', () => {
    const rand = seeded(5);
    for (let t = 0; t < 30; t++) {
      const rows = randomTable(rand);
      const expected = summary(ranges(read({ format: 'tsv', text: delimited(rows, '\t') })));
      // the same rows written in each of the four conventions: every field read alike
      for (const base of [0, 1] as const)
        for (const ends of ['inclusive', 'exclusive'] as const) {
          const written = rows.map((r) => ({
            ...r,
            start: r.start - 1 + base,
            end: r.end - 1 + base + (ends === 'exclusive' ? 1 : 0),
          }));
          const r = read({ format: 'tsv', text: delimited(written, '\t') }, { base, ends });
          expect(summary(ranges(r)), `${base} ${ends}`).toEqual(expected);
        }
      // and inline, whose defaults are 0 and exclusive
      const shifted = rows.map((r) => ({ ...r, start: r.start - 1 }));
      const inline = read({ value: shifted as unknown as Value }, { types: { score: 'category' } });
      const fromInline = summary(ranges(inline));
      expect({ ...fromInline, fields: { ...fromInline.fields, score: expected.fields.score } }).toEqual(expected);
    }
  });

  it('maps roles and fields by column or path; default fields are the scalar members', () => {
    const doc = {
      sequence: { value: 'MEEPQSDPSVEPPLSQ' },
      features: [
        {
          type: 'Domain',
          description: 'TAD',
          featureId: 'F1',
          location: { start: { value: 3, modifier: 'EXACT' }, end: { value: 7, modifier: 'OUTSIDE' } },
          evidences: [{ code: 'ECO' }],
        },
        {
          type: 'Region',
          description: 'Disordered',
          featureId: 'F2',
          location: { start: { value: 1, modifier: 'EXACT' }, end: { value: 2, modifier: 'EXACT' } },
        },
      ],
    };
    const r = read(
      { format: 'json', text: JSON.stringify(doc) },
      {
        records: 'features',
        letters: 'sequence.value',
        system: 'P04637',
        fields: {
          start: 'location.start.value',
          end: 'location.end.value',
          id: 'featureId',
          category: 'type',
          start_modifier: 'location.start.modifier',
          end_modifier: 'location.end.modifier',
        },
      },
    );
    expect(r.issues).toEqual([]);
    expect(r.data!.coverage.systems).toEqual([{ name: 'P04637', length: 16 }]);
    const t = summary(ranges(r, 'P04637'));
    expect(t.start).toEqual([0, 2]);
    expect(t.end).toEqual([2, 7]);
    expect(t.id).toEqual(['F2', 'F1']);
    expect(t.row).toEqual([1, 0]);
    // mapped columns are not fields again; objects and lists (location, evidences) are not fields
    expect(Object.keys(t.fields).sort()).toEqual(['category', 'description', 'end_modifier', 'start_modifier']);
    expect(t.fields.end_modifier.v).toEqual(['EXACT', 'OUTSIDE']);
  });

  it('puts rows on their target’s system, an empty target on the source’s own', () => {
    const r = read(
      { format: 'tsv', text: 'chrom\tstart\tend\nchr2\t5\t6\n\t1\t2\nchr2\t1\t1\n' },
      { fields: { target: 'chrom' } },
    );
    expect(r.data!.coverage.systems.map((s) => s.name)).toEqual(['chr2', 's']);
    expect(summary(ranges(r, 'chr2')).start).toEqual([0, 4]);
    expect(summary(ranges(r, 's')).start).toEqual([0]);
  });

  it('lists the source’s own system only when a row names no target, it holds letters, or no row names any', () => {
    // every row with a target: no empty own system, and one system keeps a typed column as it is
    const start = Int32Array.of(0, 4);
    const all = read(
      { value: { chrom: ['chrX', 'chrX'], start, end: Int32Array.of(3, 6) } as unknown as Value },
      {
        fields: { target: 'chrom' },
      },
    );
    expect(all.data!.coverage.systems.map((s) => s.name)).toEqual(['chrX']);
    expect(ranges(all, 'chrX').start).toBe(start);
    const letters = read(
      { format: 'json', text: '{"seq":"ACGT","rows":[{"t":"x","start":1,"end":2}]}' },
      { records: 'rows', letters: 'seq', fields: { target: 't' } },
    );
    expect(letters.data!.coverage.systems).toEqual([{ name: 's', length: 4 }, { name: 'x' }]);
    const empty = read({ format: 'tsv', text: 'start\tend\n' });
    expect(empty.data!.coverage.systems).toEqual([{ name: 's' }]);
  });

  it('reads values sparse by position or a row each, `__proto__` and `constructor` ordinary names', () => {
    const sparse = values(read({ format: 'csv', text: 'position,value\n3,0.5\n5,NaN\n6,2\n' }));
    expect(sparse.start).toBe(2);
    expect([...sparse.value]).toEqual([0.5, Number.NaN, Number.NaN, 2]);
    const dense = values(read({ format: 'json', text: '{"value":[1,2,3]}' }));
    expect([dense.start, ...dense.value]).toEqual([0, 1, 2, 3]);
    const odd = read({ format: 'json', text: '[{"start":1,"end":2,"__proto__":"a","constructor":"b"}]' });
    const oddFields = ranges(odd).fields ?? [];
    expect(oddFields.map((f) => [f.name, plain(f.column)])).toEqual([
      ['__proto__', ['a']],
      ['constructor', ['b']],
    ]);
    // a missing `constructor` member is no field: it is not read from the prototype
    const none = read({ format: 'json', text: '[{"start":1,"end":2}]' }, { fields: { label: 'constructor' } });
    expect(codes(none.issues)).toEqual(['seq.source.records.field']);
    // nor is a transform's map read through its prototype: a column named `constructor`, with `types` and `fields`
    for (const transform of [{ types: {} }, { fields: {} }, { types: {}, fields: { label: 'constructor' } }]) {
      const t = read({ format: 'tsv', text: 'start\tend\tconstructor\n1\t2\tx\n' }, transform as RecordsTransform);
      expect(t.issues).toEqual([]);
      expect(plain(ranges(t).fields?.[0].column)).toEqual(['x']);
    }
  });

  it('values without `position` are each at its row’s own index; with `start` and `end`, a field only', () => {
    const v = values(read({ format: 'tsv', text: 'value\n1\nabc\n3\n' }));
    expect([v.start, ...v.value]).toEqual([0, 1, Number.NaN, 3]);
    const r = read({ format: 'tsv', text: 'start\tend\tvalue\n1\t2\t0.5\n' });
    expect(r.data!.pairs[0].request.shapes).toEqual(['ranges']);
    expect(r.data!.coverage.shapes.map((x) => x.shape)).toEqual(['ranges']);
    expect(summary(ranges(r)).fields.value.v).toEqual([0.5]);
  });

  it('types fields as declared; strands from their spellings', () => {
    const r = read(
      { format: 'tsv', text: 'start\tend\tstrand\tn\tf\tc\n1\t2\t+\t7\t0.5\tx\n3\t4\t-1\t8\t1\t9\n5\t6\t?\t9\t2\ty\n' },
      { types: { n: 'int32', f: 'float32', c: 'string' } },
    );
    const t = ranges(r);
    expect([...t.strand!]).toEqual([1, -1, 0]);
    expect(summary(t).fields).toEqual({
      n: { type: 'Int32Array', v: [7, 8, 9] },
      f: { type: 'Float32Array', v: [0.5, 1, 2] },
      c: { type: 'CategoryColumn', v: ['x', '9', 'y'] },
    });
  });
});

// ---- typed columns kept by identity ----

describe('typed columns', () => {
  it('a typed aligned SoA column in the protocol’s conventions is kept by identity', () => {
    const start = Int32Array.of(0, 4, 9);
    const end = Int32Array.of(3, 6, 12);
    const id = categoryColumn(['a', 'b', 'c']);
    const value = Float32Array.of(0.5, Number.NaN, 2);
    const t = ranges(read({ value: { start, end, id, value } as unknown as Value }));
    expect(t.start).toBe(start);
    expect(t.end).toBe(end);
    expect(t.id).toBe(id);
    expect(t.fields?.find((f) => f.name === 'value')?.column).toBe(value);
    // read by region through the join: the same columns
    const s = openSeq(read({ value: { start, end } as unknown as Value }).data!);
    return s.read({ system: 's', start: 0, end: 1024, level: 1, shapes: ['ranges'] }, signal).then((a) => {
      expect(a.ok && (a.tables.ranges as RangesList).start).toBe(start);
    });
  });

  it('a values column of every row is kept; one that needs work is copied once', () => {
    const v = Float32Array.of(1, Number.NaN, 3);
    expect(values(read({ value: { value: v } as unknown as Value })).value).toBe(v);
    // unsorted, another base, a row left out, another type: a copy
    const start = Int32Array.of(4, 0);
    const end = Int32Array.of(6, 3);
    expect(ranges(read({ value: { start, end } as unknown as Value })).start).not.toBe(start);
    const sorted = Int32Array.of(0, 4);
    expect(
      ranges(read({ value: { start: sorted, end: Int32Array.of(3, 6) } as unknown as Value }, { base: 1 })).start,
    ).not.toBe(sorted);
    expect(ranges(read({ value: { start: sorted, end: Int32Array.of(3, 2) } as unknown as Value })).start).not.toBe(
      sorted,
    );
    const wide = Float64Array.of(0, 4);
    const t = ranges(read({ value: { start: wide, end: Float64Array.of(3, 6) } as unknown as Value }));
    expect(t.start).toBeInstanceOf(Int32Array);
    expect([...t.start]).toEqual([0, 4]);
  });

  it('a MessagePack column aligned in its document is a view into the file, kept; unaligned it is one copy', () => {
    const transform = { records: 'params', base: 0, ends: 'exclusive' } as const;
    // a longer kind moves the column's payload: some length aligns it to 4 bytes, the next does not
    const views: boolean[] = [];
    for (let pad = 0; pad < 4; pad++) {
      const node = {
        kind: 'x'.repeat(pad + 1),
        params: { start: Int32Array.of(0, 2, 5), end: Int32Array.of(1, 3, 9) },
      };
      const bytes = MsgPackCodec.encode(node);
      const t = ranges(read({ format: 'msgpack', bytes }, transform));
      expect([...t.start]).toEqual([0, 2, 5]);
      views.push(t.start.buffer === bytes.buffer);
    }
    expect(views).toContain(true);
    expect(views).toContain(false);
  });
});

// ---- issues ----

describe('issues', () => {
  const issue = (r: ReturnType<typeof read>, code: string) => r.issues.find((i) => i.code === code);

  it('malformed: a file that cannot be read as its format', () => {
    for (const [input, transform, why] of [
      [{ format: 'json', text: '{"a": [1, 2' }, {}, /does not parse/],
      [{ format: 'json', text: '{"a":[1,2],"b":[1]}' }, {}, /unequal lengths/],
      [{ format: 'json', text: '{"a":{"b":1}}' }, {}, /not a list/],
      [{ format: 'json', text: '{"a":3}' }, { records: 'a' }, /reaches neither/],
      [{ format: 'json', text: '[1, 2]' }, {}, /not an object/],
      [{ format: 'tsv', text: 'a\ta\n1\t2\n' }, {}, /twice/],
      [{ format: 'tsv', text: '' }, {}, /no header/],
      [{ format: 'csv', text: 'a\n"x\n' }, {}, /unterminated/i],
      [{ format: 'msgpack', bytes: Uint8Array.of(0x92, 0x01) }, {}, /does not parse/],
      [{ format: 'json', text: '[{"start":1,"end":2}]' }, { records: Array(40).fill('a').join('.') }, /segments/],
      [
        { format: 'json', text: '{"l":3,"rows":[{"start":1,"end":2}]}' },
        { letters: 'l', records: 'rows' },
        /no string/,
      ],
    ] as const) {
      const r = read(input as RecordsInput, transform as RecordsTransform);
      expect(r.data).toBeUndefined();
      expect(codes(r.issues), String(why)).toEqual(['seq.source.records.malformed']);
      expect(r.issues[0].message).toMatch(why);
    }
  });

  it('malformed: a MessagePack document deeper than the stack (its decoder recurses), never a crash', () => {
    const bytes = new Uint8Array(200_001).fill(0x91);
    bytes[200_000] = 0xc0;
    const r = read({ format: 'msgpack', bytes });
    expect(codes(r.issues)).toEqual(['seq.source.records.malformed']);
    expect(r.issues[0].message).toMatch(/nested too deeply/);
  });

  it('malformed: a deeply nested JSON document (JSON.parse does not recurse) is records of no object', () => {
    const deep = `${'['.repeat(200_000)}${']'.repeat(200_000)}`;
    const r = read({ format: 'json', text: deep });
    expect(codes(r.issues)).toEqual(['seq.source.records.malformed']);
    expect(r.issues[0].message).toMatch(/not an object/);
  });

  it('value: an end past int32 leaves the row out', () => {
    const r = read({ format: 'tsv', text: 'start\tend\n1\t2147483647\n1\t2\n' }, { base: 0 });
    expect(summary(ranges(r)).start).toEqual([1]);
    expect(issue(r, 'seq.source.records.value')?.message).toMatch(/past int32: 1 row \(line 2\)/);
  });

  it('one severity a code: every issue has its catalog’s (EC-87)', () => {
    const r = read({
      format: 'json',
      text: '[{"start":1,"end":2,"strand":"x","l":{}},{"start":"a","end":2},{"start":3,"end":1},{"end":2}]',
    });
    expect(new Set(codes(r.issues))).toEqual(
      new Set([
        'seq.source.records.cell',
        'seq.source.records.value',
        'seq.source.records.order',
        'seq.source.records.missing',
      ]),
    );
    for (const i of r.issues)
      expect(i.severity, i.code).toBe(RecordsIssues.codes[i.code as keyof typeof RecordsIssues.codes].severity);
    const twice = read({ format: 'tsv', text: 'position\tvalue\n1\t1\n1\t2\n' });
    expect(codes(twice.issues)).toEqual(['seq.source.records.cell']);
  });

  it('field: a column found in no record, listing the columns there are', () => {
    const r = read({ format: 'tsv', text: 'first\tlast\n1\t2\n' }, { fields: { start: 'begin', end: 'last' } });
    expect(codes(r.issues)).toEqual(['seq.source.records.field']);
    expect(r.issues[0].message).toMatch(/“start” reads “begin”.*there are: first, last/);
    expect(r.issues[0].hint).toBeDefined();
    expect(codes(read({ format: 'tsv', text: 'a\tb\n1\t2\n' }).issues)).toEqual(['seq.source.records.field']);
    expect(codes(read({ format: 'tsv', text: 'start\n1\n' }).issues)).toEqual(['seq.source.records.field']);
  });

  it('row-level misfits: each left out, one issue a cause with its count and first ten places', () => {
    const lines = ['start\tend\tvalue', '1\t2\t0.5', 'x\t3\t1', '\t4\t1', '5\t3\t1', '1\t2\tabc'];
    for (let i = 0; i < 12; i++) lines.push('y\t1\t1');
    const r = read({ format: 'tsv', text: `${lines.join('\n')}\n` });
    expect(summary(ranges(r)).start).toEqual([0]);
    const notInt = r.issues.find((i) => i.message.startsWith('A position is not an integer'))!;
    expect(notInt.code).toBe('seq.source.records.value');
    expect(notInt.message).toMatch(/13 rows \(line 3, line 7, .*line 15, …\)/);
    expect(issue(r, 'seq.source.records.missing')?.message).toMatch(/1 row \(line 4\)/);
    expect(issue(r, 'seq.source.records.order')?.message).toMatch(/1 row \(line 5\)/);
    expect(r.issues.find((i) => i.message.startsWith('A value is not a number'))?.message).toMatch(/line 6/);
    // UniProt's UNKNOWN: a null position
    const unknown = read({ format: 'json', text: '[{"start":null,"end":4},{"start":1,"end":2}]' });
    expect(issue(unknown, 'seq.source.records.missing')?.message).toMatch(/record 0/);
  });

  it('range: past the letters left out or clipped, counted in one issue, hinting the other base', () => {
    const doc = (rows: [number, number][]) =>
      JSON.stringify({ seq: 'ACGTACGTAC', rows: rows.map(([start, end]) => ({ start, end })) });
    const r = read(
      {
        format: 'json',
        text: doc([
          [1, 4],
          [9, 12],
          [11, 12],
        ]),
      },
      { records: 'rows', letters: 'seq' },
    );
    expect(summary(ranges(r)).end).toEqual([4, 10]);
    const range = issue(r, 'seq.source.records.range')!;
    expect(range.message).toMatch(/2 rows/);
    expect(range.hint).toBeUndefined();
    // every row one past the end: a 0-based file read as 1-based would never do that; a 1-based one read as 0 would
    const shifted = read(
      {
        format: 'json',
        text: doc([
          [1, 4],
          [8, 10],
        ]),
      },
      { records: 'rows', letters: 'seq', base: 0, ends: 'inclusive' },
    );
    expect(issue(shifted, 'seq.source.records.range')?.hint).toMatch(/base: 1/);
    // a first position 0 with base 1
    const zero = read({ format: 'tsv', text: 'start\tend\n0\t3\n4\t5\n' });
    expect(issue(zero, 'seq.source.records.range')?.hint).toMatch(/base: 0/);
  });

  it('cells: a ragged row read with empty cells, said once', () => {
    const r = read({ format: 'csv', text: 'start,end,label\n1,2\n3,4,x,extra\n5,6,y\n' });
    expect(summary(ranges(r)).fields.label.v).toEqual(['', 'x', 'y']);
    expect(issue(r, 'seq.source.records.cells')?.message).toMatch(/2 rows \(line 2, line 3\)/);
  });

  it('encoding: replacement characters said once; a path to an object a null field, said', () => {
    const r = read({ format: 'tsv', text: 'start\tend\tlabel\n1\t2\tcaf�\n' });
    expect(issue(r, 'seq.source.records.encoding')?.severity).toBe('info');
    const o = read(
      { format: 'json', text: '[{"start":1,"end":2,"l":{"a":1}},{"start":2,"end":3,"l":"x"}]' },
      { fields: { label: 'l' } },
    );
    expect(summary(ranges(o)).fields.label.v).toEqual([null, 'x']);
    expect(issue(o, 'seq.source.records.cell')?.severity).toBe('info');
  });

  it('the runtime: refusals as SourceOpenError with the address, row issues on the opened source', async () => {
    await expect(
      RecordsRuntime.open({ address: '/data/t', format: 'json', uri: 'f' }, ctx({ f: '{' })),
    ).rejects.toMatchObject({ issues: [{ code: 'seq.source.records.malformed', address: '/data/t' }] });
    await expect(RecordsRuntime.open({ address: '/data/t', format: 'json', uri: 'none' }, ctx())).rejects.toMatchObject(
      {
        issues: [{ code: 'seq.source.unreadable' }],
      },
    );
    await expect(RecordsRuntime.open({ address: '/data/t', format: 'json' }, ctx())).rejects.toBeInstanceOf(
      SourceOpenError,
    );
    const s = (await RecordsRuntime.open(
      { address: '/data/t', format: 'csv', uri: 'f' },
      ctx({ f: 'start,end\n1,2\n4,3\n' }),
    )) as RecordsSource;
    expect(s.issues).toMatchObject([{ code: 'seq.source.records.order', address: '/data/t' }]);
    // the source's id names its own system
    expect(s.coverage.systems[0].name).toBe('t');
  });
});

// ---- fuzz ----

describe('fuzz', () => {
  it('reads mutated text and documents without a throw', () => {
    const rand = seeded(17);
    const rows = randomTable(rand, 8);
    const sources = encodings(rows).filter((e) => !('value' in e.input));
    const alphabet = ['"', ',', '\t', '\n', '\r', '{', '}', '[', ']', ':', '0', '-', 'a', '�', ' '];
    for (let t = 0; t < 400; t++) {
      const e = sources[t % sources.length];
      // the transform varied too: swapped roles, types of any name (`constructor` among them), paths
      const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
      const names = ['start', 'end', 'label', 'score', 'constructor', '__proto__', 'toString', 'x.y', '0'];
      const transform: RecordsTransform = {
        ...e.transform,
        ...(rand() < 0.3 && { fields: { start: 'end', end: 'start', [pick(names)]: pick(names) } }),
        ...(rand() < 0.3 && {
          types: { [pick(names)]: pick(['int32', 'float32', 'float64', 'category', 'string'] as const) },
        }),
        ...(rand() < 0.2 && { records: pick(['features', 'params', '0', 'constructor', 'a.b.c', '']) }),
        ...(rand() < 0.2 && { letters: pick(['0.label', 'constructor', 'features.0', 'x']) }),
      };
      const input = e.input;
      if ('text' in input) {
        const chars = [...input.text];
        for (let m = 0; m < 1 + rand() * 5; m++) {
          const at = Math.floor(rand() * (chars.length + 1));
          if (rand() < 0.5) chars.splice(at, 1);
          else chars.splice(at, 0, alphabet[Math.floor(rand() * alphabet.length)]);
        }
        expect(() => read({ format: input.format, text: chars.join('') }, transform)).not.toThrow();
      } else if ('bytes' in e.input) {
        const bytes = Uint8Array.from(e.input.bytes);
        for (let m = 0; m < 1 + rand() * 5; m++) bytes[Math.floor(rand() * bytes.length)] = Math.floor(rand() * 256);
        const cut = rand() < 0.3 ? bytes.subarray(0, Math.floor(rand() * bytes.length)) : bytes;
        expect(() => read({ format: 'msgpack', bytes: cut }, transform)).not.toThrow();
      }
    }
  });
});

// ---- bounds ----

describe('bounds', () => {
  it('a file past maxBytes is refused before it is parsed', async () => {
    const big = `[${'{"start":1,"end":2},'.repeat(100)}`; // not even JSON: refused by its size, never parsed
    await expect(
      RecordsRuntime.open(
        { address: '/data/t', format: 'json', uri: 'f' },
        ctx({ f: big }, { maxRows: 10, maxBytes: 100 }),
      ),
    ).rejects.toMatchObject({ issues: [{ code: 'seq.source.too_large', hint: expect.stringMatching(/region/) }] });
  });

  it('rows past maxRows are refused before columns are built', () => {
    const text = `start\tend\n${'1\t2\n'.repeat(11)}`;
    const r = readRecords({ format: 'tsv', text }, {}, { system: 's', limits: { maxRows: 10, maxBytes: 1e9 } });
    expect(codes(r.issues)).toEqual(['seq.source.too_large']);
    const span = readRecords(
      { format: 'tsv', text: 'position\tvalue\n1\t1\n100\t1\n' },
      {},
      { system: 's', limits: { maxRows: 10, maxBytes: 1e9 } },
    );
    expect(codes(span.issues)).toEqual(['seq.source.too_large']);
  });
});

// ---- the open's target ----

describe('the open', () => {
  it('reads 10,000 TSV rows', () => {
    const n = 10_000;
    let text = 'chrom\tstart\tend\tstrand\tname\tscore\n';
    for (let i = 0; i < n; i++)
      text += `\t${((i * 7919) % 100_000) + 1}\t${((i * 7919) % 100_000) + 50}\t+\tf${i}\t${i % 100}\n`;
    const t = ranges(read({ format: 'tsv', text }, { fields: { target: 'chrom' } }));
    expect(t.start.length).toBe(n);
  });
});
