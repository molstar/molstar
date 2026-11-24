/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { encodeMsgPack } from '@molstar/core/data/msgpack';
import { MsgPackCodec } from '@molstar/core/repr/codecs/msgpack';
import { CategoryColumn, type Column, type Value } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits, type OpenContext, type RangesList, type RegionAnswer } from '../protocol.js';
import { type RecordsInput, RecordsRuntime, type RecordsSource, readMatrix, readRecords } from '../records.js';
import { checkCoverage } from '../seq-format.js';

// The records runtime's GFF3 preset and its matrix read (_spec/README.md "Data", Records).

const limits = DefaultOpenLimits;
const signal = new AbortController().signal;

const plain = (c: Column | undefined) =>
  c === undefined
    ? undefined
    : c instanceof CategoryColumn
      ? Array.from({ length: c.length }, (_, i) => c.get(i))
      : [...c];

const ranges = (r: ReturnType<typeof readRecords>, system = 'ring'): RangesList => {
  expect(r.data, JSON.stringify(r.issues)).toBeDefined();
  const pair = r.data!.pairs.find((p) => p.request.system === system)!;
  return (pair.answer as Extract<RegionAnswer, { ok: true }>).tables.ranges as RangesList;
};

const field = (t: RangesList, name: string) => plain(t.fields?.find((f) => f.name === name)?.column);

const gff = (text: string, transform = {}) =>
  readRecords({ format: 'gff3', text }, transform, { system: 'own', limits });

const Sample = [
  '##gff-version 3',
  '##sequence-region ring 1 100',
  '# a comment',
  'ring\tseqx\tregion\t1\t100\t.\t+\t.\tID=region;Name=ring;Is_circular=true',
  'ring\tseqx\tCDS\t11\t20\t5.5\t-\t0\tID=g1;Name=gene%2C one;product=a%3Bb;Parent=r1,r2',
  'ring\tseqx\tgene\t5\t9\t.\t.\t.\tName=no id',
  'other\tseqx\tCDS\t3\t4\t.\t+\t.\tID=o1',
  '',
  '##FASTA',
  '>ring',
  'ACGTACGT',
  'ring\tseqx\tCDS\t1\t2\t.\t+\t.\tID=after-fasta',
].join('\n');

describe('the GFF3 preset', () => {
  it('reads nine columns from 1, ends inclusive: comments, directives and what follows ##FASTA read past', () => {
    const r = gff(Sample);
    expect(r.issues).toEqual([]);
    const t = ranges(r);
    // sorted by start; each row keeps its place among the rows kept
    expect([...t.start]).toEqual([0, 4, 10]);
    expect([...t.end]).toEqual([100, 9, 20]);
    expect(plain(t.id)).toEqual(['region', null, 'g1']);
    expect(plain(t.strand)).toEqual([1, 0, -1]);
    expect(r.data!.coverage.systems.map((s) => s.name).sort()).toEqual(['other', 'ring']);
  });

  it('puts a feature on its seqid’s system, one with none on the source’s own', () => {
    const r = gff('\tseqx\tCDS\t1\t5\t.\t+\t.\tID=a\nring\tseqx\tCDS\t1\t5\t.\t+\t.\tID=b');
    expect(plain(ranges(r, 'own').id)).toEqual(['a']);
    expect(plain(ranges(r, 'ring').id)).toEqual(['b']);
  });

  it('makes the other columns and the attributes fields, percent-decoded, `.` missing, multi-valued as text', () => {
    const t = ranges(gff(Sample));
    expect(field(t, 'type')).toEqual(['region', 'gene', 'CDS']);
    expect(field(t, 'source')).toEqual(['seqx', 'seqx', 'seqx']);
    expect(field(t, 'score')).toEqual([null, null, '5.5']);
    expect(field(t, 'phase')).toEqual([null, null, '0']);
    expect(field(t, 'Name')).toEqual(['ring', 'no id', 'gene, one']);
    expect(field(t, 'product')).toEqual([null, null, 'a;b']);
    // `Parent` is kept as a field: no hierarchy is assembled
    expect(field(t, 'Parent')).toEqual([null, null, 'r1,r2']);
    expect(field(t, 'Is_circular')).toEqual(['true', null, null]);
    // the roles are not fields under their own names
    expect(field(t, 'ID')).toBeUndefined();
    expect(field(t, 'seqid')).toBeUndefined();
  });

  it('maps fields as any other source does, and fixes its own base and ends', () => {
    const t = ranges(gff(Sample, { fields: { category: 'type', label: 'Name' }, base: 0, ends: 'exclusive' }));
    expect(field(t, 'category')).toEqual(['region', 'gene', 'CDS']);
    expect(field(t, 'label')).toEqual(['ring', 'no id', 'gene, one']);
    // the preset's conventions stand: the transform's `base` and `ends` are not read
    expect([...t.start]).toEqual([0, 4, 10]);
  });

  it('leaves a row with a missing position out, keeps a ragged one with its cells said once, CRLF and a BOM read', () => {
    const r = gff(
      '﻿##gff-version 3\r\nring\ts\tCDS\t.\t5\t.\t+\t.\tID=a\r\nring\ts\tCDS\t3\t4\r\nring\ts\tCDS\t6\t8\t.\t+\t.\tID=c\r\n',
    );
    const t = ranges(r);
    expect(plain(t.id)).toEqual([null, 'c']);
    expect(r.issues.map((i) => i.code).sort()).toEqual(['seq.source.records.cells', 'seq.source.records.missing']);
  });

  it('a GFF3 with no feature reads no row, as an empty table does', () => {
    const r = gff('##gff-version 3\n');
    expect(r.issues).toEqual([]);
    expect(ranges(r, 'own').start.length).toBe(0);
  });

  it('is read by the runtime from a file, and from inline text', async () => {
    const files: Record<string, string> = { 'a.gff3': Sample };
    const context: OpenContext = {
      signal,
      limits,
      assets: {
        readText: async (u) => files[u],
        readBytes: async (u) => new TextEncoder().encode(files[u]),
      },
    };
    const entry = { address: '/data/s', format: 'gff3', uri: 'a.gff3' };
    for (const source of [entry, { address: '/data/s', format: 'gff3', value: Sample }]) {
      const opened = (await RecordsRuntime.open(source, context)) as RecordsSource;
      expect(checkCoverage(opened.coverage)).toBeUndefined();
      const a = await opened.read({ system: 'ring', start: 0, end: 100, level: 1, shapes: ['ranges'] }, signal);
      expect(a.ok && (a.tables.ranges as RangesList).start.length).toBe(3);
    }
  });
});

describe('the GFF3 preset, the format’s own cases', () => {
  it('a CDS over rows sharing an ID is one item, its rows its pieces; a strand `?` is read as 0', () => {
    const r = gff(
      [
        'ring\tt\tCDS\t11\t20\t.\t-\t0\tID=cds1;Name=a',
        'ring\tt\tCDS\t41\t50\t.\t-\t1\tID=cds1;Name=a',
        'ring\tt\tgene\t1\t5\t.\t?\t.\tID=g',
      ].join('\n'),
    );
    expect(r.issues).toEqual([]);
    const t = ranges(r);
    expect(plain(t.id)).toEqual(['g', 'cds1', 'cds1']);
    expect([...t.start]).toEqual([0, 10, 40]);
    expect([...t.end]).toEqual([5, 20, 50]);
    expect(plain(t.strand)).toEqual([0, -1, -1]);
    // each row keeps its place in the file
    expect(plain(t.row)).toEqual([2, 0, 1]);
  });

  it('a repeated attribute key keeps the last, and an attribute named as a column leaves the column its own: each said once, info', () => {
    const r = gff(
      [
        'ring\tt\tCDS\t1\t5\t.\t+\t.\tID=a;Name=x;Name=y;type=other',
        'ring\tt\tCDS\t6\t9\t.\t+\t.\tID=b;Name=p;Name=q;type=more',
      ].join('\n'),
    );
    const t = ranges(r);
    expect(field(t, 'Name')).toEqual(['y', 'q']);
    expect(field(t, 'type')).toEqual(['CDS', 'CDS']);
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([
      ['seq.source.records.cell', 'info'],
      ['seq.source.records.cell', 'info'],
    ]);
    expect(r.issues.map((i) => i.message).join(' ')).toMatch(/2 rows/);
  });
});

describe('the matrix read, as a source holds it', () => {
  it('parses once at open and checks the size on each call; its issues carry the source’s address', async () => {
    let reads = 0;
    const context: OpenContext = {
      signal,
      limits,
      assets: {
        readText: async () => '',
        readBytes: async () => {
          reads++;
          return new TextEncoder().encode('{"m":[[1,2],[3,4]]}');
        },
      },
    };
    const opened = (await RecordsRuntime.open(
      { address: '/data/pae', format: 'json', uri: 'm.json', shape: 'matrix', records: 'm' },
      context,
    )) as RecordsSource;
    expect(reads).toBe(1);
    const a = opened.matrix!({ rows: 2, cols: 2 });
    const b = opened.matrix!({ rows: 3, cols: 3 });
    expect(a.matrix!.values.length).toBe(4);
    expect(b.issues.map((i) => [i.code, i.address])).toEqual([['seq.source.records.size', '/data/pae']]);
    expect(reads).toBe(1);
    // a path that reaches nothing: the same issue on every call, with the address
    const lost = (await RecordsRuntime.open(
      { address: '/data/pae', format: 'json', uri: 'm.json', shape: 'matrix', records: 'none' },
      context,
    )) as RecordsSource;
    expect(lost.matrix!().issues.map((i) => [i.code, i.address])).toEqual([['seq.source.records.field', '/data/pae']]);
  });

  it('is flat only as a typed array or an array of numbers and nulls', () => {
    const flat = readMatrix({ value: [null, 1, 2, 3] as unknown as Value }, {}, { limits });
    expect([flat.matrix!.rows, flat.matrix!.cols]).toEqual([2, 2]);
    // a column of one-element rows is rows, not a flat array
    const column = readMatrix({ value: [[1], [2], [3]] as unknown as Value }, {}, { limits });
    expect([column.matrix!.rows, column.matrix!.cols]).toEqual([3, 1]);
    // a mix is no flat array and no array of arrays
    const mixed = readMatrix({ value: [1, [2]] as unknown as Value }, {}, { limits });
    expect(mixed.issues.map((i) => i.code)).toEqual(['seq.source.records.size']);
  });
});

describe('the matrix read', () => {
  const aoa = [
    [0, 1, 2],
    [3, null, 5],
  ];
  const size = { rows: 2, cols: 3 };

  it('reads an array of arrays, JSON or MessagePack, by the records path, `null` as NaN', () => {
    const inputs: [RecordsInput, string?][] = [
      [{ format: 'json', text: JSON.stringify(aoa) }],
      [
        { format: 'json', text: JSON.stringify({ pae: { predicted_aligned_error: aoa } }) },
        'pae.predicted_aligned_error',
      ],
      [{ format: 'json', text: JSON.stringify([{ predicted_aligned_error: aoa }]) }, '0.predicted_aligned_error'],
      [{ format: 'msgpack', bytes: encodeMsgPack({ m: aoa }) }, 'm'],
      [{ value: aoa as unknown as Value }],
    ];
    for (const [input, records] of inputs) {
      const r = readMatrix(input, records ? { records } : {}, { limits, expected: size });
      expect(r.issues).toEqual([]);
      expect([r.matrix!.rows, r.matrix!.cols]).toEqual([2, 3]);
      expect([...r.matrix!.values].map((v) => (Number.isNaN(v) ? 'nan' : v))).toEqual([0, 1, 2, 3, 'nan', 5]);
    }
  });

  it('checks the size against the axes’ total: an issue on a mismatch, the matrix not served', () => {
    const wrong = readMatrix({ value: aoa as unknown as Value }, {}, { limits, expected: { rows: 3, cols: 3 } });
    expect(wrong.matrix).toBeUndefined();
    expect(wrong.issues.map((i) => i.code)).toEqual(['seq.source.records.size']);
    expect(wrong.issues[0].message).toContain('2 × 3');
    expect(wrong.issues[0].message).toContain('3 × 3');
    // rows of unequal lengths
    const ragged = readMatrix({ value: [[1, 2], [3]] as unknown as Value }, {}, { limits });
    expect(ragged.issues.map((i) => i.code)).toEqual(['seq.source.records.size']);
  });

  it('keeps a typed flat Float32Array as it is, its dimensions the axes’; a flat plain array is copied', () => {
    const flat = Float32Array.from([1, 2, 3, 4, 5, 6]);
    const node = { kind: 'm', params: { values: flat } };
    const r = readMatrix({ value: { v: flat } as unknown as Value }, { records: 'v' }, { limits, expected: size });
    expect(r.matrix!.values).toBe(flat);
    expect([r.matrix!.rows, r.matrix!.cols]).toEqual([2, 3]);
    // through MessagePack's typed column: a view into the document, equal
    const packed = readMatrix(
      { format: 'msgpack', bytes: MsgPackCodec.encode(node) },
      { records: 'params.values' },
      { limits, expected: size },
    );
    expect(packed.issues).toEqual([]);
    expect([...packed.matrix!.values]).toEqual([...flat]);
    // a square when no size is given; a flat array that is not one is refused
    const square = readMatrix({ value: { v: [1, 2, 3, 4] } as unknown as Value }, { records: 'v' }, { limits });
    expect([square.matrix!.rows, square.matrix!.cols]).toEqual([2, 2]);
    const bad = readMatrix({ value: { v: [1, 2, 3] } as unknown as Value }, { records: 'v' }, { limits });
    expect(bad.issues.map((i) => i.code)).toEqual(['seq.source.records.size']);
    // a flat array of the wrong size for the axes
    const off = readMatrix(
      { value: { v: flat } as unknown as Value },
      { records: 'v' },
      { limits, expected: { rows: 4, cols: 4 } },
    );
    expect(off.issues.map((i) => i.code)).toEqual(['seq.source.records.size']);
  });

  it('refuses a path that reaches no array and a document that does not parse', () => {
    expect(readMatrix({ format: 'json', text: '{"a": 1}' }, { records: 'a' }, { limits }).issues[0].code).toBe(
      'seq.source.records.malformed',
    );
    expect(readMatrix({ format: 'json', text: '{"a": 1}' }, { records: 'b' }, { limits }).issues[0].code).toBe(
      'seq.source.records.field',
    );
    expect(readMatrix({ format: 'json', text: '[[1,' }, {}, { limits }).issues[0].code).toBe(
      'seq.source.records.malformed',
    );
    // an empty array is an empty matrix
    const empty = readMatrix({ value: [] as unknown as Value }, {}, { limits });
    expect(empty.matrix).toEqual({ rows: 0, cols: 0, values: new Float32Array(0) });
  });

  it('a source of `shape: matrix` serves its matrix to a grid, and no region', async () => {
    const context: OpenContext = {
      signal,
      limits,
      assets: { readText: async () => '', readBytes: async () => encodeMsgPack({ pae: aoa }) },
    };
    const opened = (await RecordsRuntime.open(
      { address: '/data/pae', format: 'msgpack', uri: 'pae.msgpack', shape: 'matrix', records: 'pae' },
      context,
    )) as RecordsSource;
    expect(opened.coverage.shapes).toEqual([]);
    expect(opened.matrix!(size).matrix!.values.length).toBe(6);
    expect(opened.matrix!({ rows: 9, cols: 9 }).issues[0].code).toBe('seq.source.records.size');
    // a table is no matrix
    await expect(
      RecordsRuntime.open({ address: '/data/pae', format: 'tsv', value: 'a\n1', shape: 'matrix' }, context),
    ).rejects.toMatchObject({ issues: [{ code: 'seq.source.no_data' }] });
  });
});
