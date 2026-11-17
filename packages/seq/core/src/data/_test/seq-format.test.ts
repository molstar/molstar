/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync } from 'node:fs';
import { categoryColumn, hashNode, issueOf } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { SourceIssues, StoreIssues } from '../issues.js';
import type { RegionAnswer, RegionRequest, Tables } from '../protocol.js';
import { decodeSeq, encodeSeq, inlineSeq, readSeq, type SeqData, seqToNode, seqToValue } from '../seq-format.js';
import { openSeq, SeqRuntime, SourceOpenError } from '../seq-runtime.js';
import { answerBytes, checkTables, isEmptyAnswer, sliceTable } from '../tables.js';

const transcript = { kind: 'sequence', id: 'T1', params: { name: 'T1', length: 300 } };

/** Every shape and form, odd values, nodes and a failure: what a file must carry back. */
const data: SeqData = {
  coverage: {
    systems: [
      { name: 'chr1', length: 4000, aliases: ['1'] },
      { name: 'p53', alphabet: 'protein', cyclic: false },
    ],
    shapes: [
      { system: 'chr1', shape: 'ranges', fields: [{ name: 'score', type: 'float32' }] },
      { system: 'chr1', shape: 'values' },
      { system: 'p53', shape: 'letters' },
    ],
    levels: [1, 10],
    nodes: [{ key: 'T1', node: transcript }],
  },
  pairs: [
    {
      request: { system: 'chr1', start: 0, end: 4000, level: 1, shapes: ['ranges', 'values'] },
      answer: {
        ok: true,
        tables: {
          ranges: {
            form: 'rows',
            start: Int32Array.of(100, 1000, 2500),
            end: Int32Array.of(400, 2100, 2600),
            strand: Int8Array.of(1, -1, 0),
            id: categoryColumn(['a', 'b', null]),
            fields: { score: Float32Array.of(0.5, Number.NaN, -0) },
          },
          values: { form: 'positions', start: 0, value: Float32Array.from({ length: 4000 }, (_, i) => i % 7) },
        },
        nodes: [{ key: 'T1', node: transcript }],
      },
    },
    {
      request: { system: 'chr1', start: 0, end: 4000, level: 10, shapes: ['ranges', 'values'] },
      answer: {
        ok: true,
        tables: {
          ranges: { form: 'density', bin: Int32Array.of(10, 100, 250), count: Float32Array.of(1, 2, 1) },
          values: {
            form: 'bins',
            start: 0,
            min: new Float32Array(400),
            max: new Float32Array(400).fill(6),
            mean: new Float32Array(400).fill(3),
            count: new Float32Array(400).fill(10),
          },
        },
      },
    },
    {
      request: { system: 'p53', start: 0, end: 8, level: 1, shapes: ['letters'] },
      answer: { ok: true, tables: { letters: { form: 'letters', start: 0, letters: 'MEEPQSDP' } } },
    },
    {
      request: { system: 'p53', start: 8, end: 12, level: 1, shapes: ['letters'] },
      answer: {
        ok: false,
        kind: 'failed',
        issue: issueOf(SourceIssues, 'seq.source.failed', 'Not there', { path: ['x', 1] }),
      },
    },
  ],
};

describe('the seq format', () => {
  it('round-trips through JSON: typed columns, NaN, -0, categories, nodes and failures', () => {
    const text = encodeSeq(data);
    const back = decodeSeq(text);
    expect(back.issues).toEqual([]);
    expect(hashNode(seqToNode(back.data!))).toBe(hashNode(seqToNode(data)));
    const ranges = back.data!.pairs[0].answer.ok && back.data!.pairs[0].answer.tables.ranges;
    expect(ranges && ranges.form === 'rows' && ranges.fields?.score).toBeInstanceOf(Float32Array);
    expect(back.data!.pairs[3].answer).toMatchObject({
      ok: false,
      kind: 'failed',
      issue: { code: 'seq.source.failed' },
    });
    expect(encodeSeq(back.data!)).toBe(text);
  });

  it('reads inline data as the record a spec value holds', () => {
    const r = readSeq(seqToValue(data));
    expect(r.issues).toEqual([]);
    expect(hashNode(seqToNode(r.data!))).toBe(hashNode(seqToNode(data)));
  });

  it('reports the first misfit at its path, never throws', () => {
    const at = (v: Parameters<typeof readSeq>[0]) => readSeq(v).issues[0];
    const value = seqToValue(data) as Record<string, any>;
    expect(at({ ...value, coverage: { ...value.coverage, levels: [10, 1] } })).toMatchObject({
      code: 'seq.source.malformed',
      severity: 'error',
      path: ['coverage', 'levels', 1],
    });
    expect(
      at({ ...value, coverage: { ...value.coverage, shapes: [{ system: 'chrX', shape: 'ranges' }] } }).path,
    ).toEqual(['coverage', 'shapes', 0, 'system']);
    const pair = {
      request: { system: 'chr1', start: 5, end: 1, level: 1, shapes: [] },
      answer: { ok: true, tables: {} },
    };
    expect(at({ ...value, pairs: [pair] }).path).toEqual(['pairs', 0, 'request', 'end']);
    const unequal = { form: 'rows', start: Int32Array.of(1, 2), end: Int32Array.of(3) };
    const bad = { request: data.pairs[0].request as any, answer: { ok: true, tables: { ranges: unequal } } };
    expect(at({ ...value, pairs: [bad] }).path).toEqual(['pairs', 0, 'answer', 'tables']);
    expect(at({ kind: 'other' }).path).toEqual(['kind']);
    expect(decodeSeq('{"kind": "seq", "params": {').issues[0].code).toMatch(/^repr\.codec\.json\./);
  });
});

describe('tables', () => {
  const rows: Tables['ranges'] = {
    form: 'rows',
    start: Int32Array.of(0, 900, 1500),
    end: Int32Array.of(100, 1100, 1600),
    id: categoryColumn(['a', 'b', 'c']),
  };

  it('are cut to a region: ranges by overlap, positions by offset', () => {
    const cut = sliceTable(rows!, 1024, 2048, 1);
    expect(cut.form === 'rows' && Array.from(cut.start)).toEqual([900, 1500]);
    expect(cut.form === 'rows' && cut.id?.dict).toEqual(['b', 'c']);
    const letters = sliceTable({ form: 'letters', start: 10, letters: 'ABCDEFGH' }, 12, 15, 1);
    expect(letters).toEqual({ form: 'letters', start: 12, letters: 'CDE' });
    const bins = sliceTable(
      {
        form: 'bins',
        start: 0,
        min: new Float32Array(8),
        max: new Float32Array(8),
        mean: new Float32Array(8),
        count: Float32Array.of(0, 1, 2, 3, 4, 5, 6, 7),
      },
      20,
      50,
      10,
    );
    expect(bins.start).toBe(20);
    expect(Array.from(bins.count)).toEqual([2, 3, 4]);
    expect(sliceTable(rows!, 0, 5000, 1)).toBe(rows);
  });

  it('are checked for their forms, types and lengths', () => {
    expect(checkTables({ ranges: rows })).toBeUndefined();
    expect(checkTables({ ranges: { ...rows!, end: Int32Array.of(1) } as any })).toMatch(/unequal/);
    expect(checkTables({ values: { form: 'rows' } as any })).toMatch(/no form/);
    expect(checkTables({ letters: { form: 'letters', start: 0, letters: Uint8Array.of(1) } })).toMatch(/alphabet/);
    expect(checkTables({ matrix: {} } as any)).toMatch(/shape/);
  });

  it('cost their bytes, and are empty with no rows and no nodes', () => {
    expect(answerBytes({ ranges: rows })).toBe(12 + 12 + 4 * 3 + 6);
    expect(isEmptyAnswer({ ranges: sliceTable(rows!, 200, 800, 1) }, undefined)).toBe(true);
    expect(isEmptyAnswer({}, [{ key: 'T1', node: transcript }])).toBe(false);
  });
});

describe('the seq runtime', () => {
  const source = openSeq(data);
  const signal = new AbortController().signal;
  const read = (r: Partial<RegionRequest>): Promise<RegionAnswer> =>
    source.read({ system: 'chr1', start: 0, end: 1024, level: 1, shapes: ['ranges'], ...r }, signal);

  it('answers a tile from the pair that holds it, cut to it and to the shapes asked', async () => {
    const a = await read({ system: '1', start: 1024, end: 2048 });
    expect(a.ok && Object.keys(a.tables)).toEqual(['ranges']);
    expect(a.ok && a.tables.ranges?.form === 'rows' && Array.from(a.tables.ranges.start)).toEqual([1000]);
    expect(a.ok && a.nodes?.[0].key).toBe('T1');
    const coarse = await read({ start: 0, end: 10_240, level: 10, shapes: ['values'] });
    expect(coarse.ok && coarse.tables.values?.form).toBe('bins');
  });

  it('is empty where no pair is, not served off its coverage, and passes a failure on', async () => {
    expect(await read({ system: 'p53', start: 100, end: 200, shapes: ['letters'] })).toEqual({ ok: true, tables: {} });
    expect(await read({ system: 'chrX' })).toMatchObject({ ok: false, kind: 'not-served' });
    expect(await read({ level: 100 })).toMatchObject({ ok: false, kind: 'not-served' });
    expect(await read({ system: 'p53', start: 8, end: 12, shapes: ['letters'] })).toMatchObject({ kind: 'failed' });
    const aborted = new AbortController();
    aborted.abort();
    await expect(
      source.read({ system: 'chr1', start: 0, end: 1, level: 1, shapes: [] }, aborted.signal),
    ).rejects.toThrow(/aborted/);
  });

  it('opens inline data, a file through the asset port, and refuses neither', async () => {
    const ctx = (texts: Record<string, string>) => ({
      signal,
      assets: { readText: async (uri: string) => texts[uri] ?? Promise.reject(new Error('no such file')) },
    });
    const inline = await SeqRuntime.open({ address: '/data/a', format: 'seq', value: seqToValue(data) }, ctx({}));
    expect(inline.coverage.levels).toEqual([1, 10]);
    const file = await SeqRuntime.open(
      { address: '/data/b', format: 'seq', uri: 'a.json' },
      ctx({ 'a.json': encodeSeq(data) }),
    );
    expect(file.coverage.systems.length).toBe(2);
    const none = SeqRuntime.open({ address: '/data/c', format: 'seq' }, ctx({}));
    await expect(none).rejects.toBeInstanceOf(SourceOpenError);
    await expect(none).rejects.toMatchObject({ issues: [{ code: 'seq.source.no_data', address: '/data/c' }] });
  });

  it('serves inline data as one region over the whole system', async () => {
    const inline = inlineSeq({ name: 'p53' }, { letters: { form: 'letters', start: 0, letters: 'MEEPQSDPSV' } });
    expect(inline.pairs[0].request).toEqual({ system: 'p53', start: 0, end: 10, level: 1, shapes: ['letters'] });
    const a = await openSeq(inline).read({ system: 'p53', start: 4, end: 1028, level: 1, shapes: ['letters'] }, signal);
    expect(a.ok && a.tables.letters?.letters).toBe('QSDPSV');
  });
});

describe('the issue catalogs', () => {
  it('keep their codes in their namespaces', () => {
    for (const c of [SourceIssues, StoreIssues])
      for (const code of Object.keys(c.codes)) expect(code.startsWith(`${c.namespace}.`)).toBe(true);
  });

  it('catalog every code the sources use', () => {
    const known = new Set([...Object.keys(SourceIssues.codes), ...Object.keys(StoreIssues.codes)]);
    const dir = new URL('../', import.meta.url);
    const used = new Set<string>();
    for (const f of readdirSync(dir).filter((n) => n.endsWith('.ts')))
      for (const m of readFileSync(new URL(f, dir), 'utf8').matchAll(/'(seq\.(?:source|store)\.[a-z_.]+)'/g))
        used.add(m[1]);
    expect(used.size).toBeGreaterThan(4);
    for (const code of used) expect(known, code).toContain(code);
  });
});
