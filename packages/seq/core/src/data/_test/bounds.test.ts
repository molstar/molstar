/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The data bounds (_spec/README.md, Data, Bounds): a column's length and a grid or matrix block's values, refused when
// an answer or a file is read, before its tables are cut, joined or kept.

import { describe, expect, it } from 'vitest';
import { dataKey } from '../keys.js';
import type { RegionRequest, Table } from '../protocol.js';
import { decodeSeq, readSeq } from '../seq-format.js';
import { DataStore } from '../store.js';
import {
  checkAnswer,
  checkTableBounds,
  DefaultAnswerBounds,
  DefaultMaxBlockValues,
  DefaultMaxColumnLength,
} from '../tables.js';
import { FakeClock, rangesAnswer, SyntheticSource, syntheticCoverage } from './synthetic.js';

/** A table that claims `rows` rows without holding them: only the bounds' arithmetic reads it. */
const claims = (form: 'positions' | 'list', rows: number): Table =>
  form === 'positions'
    ? ({ form, start: 0, value: { length: rows } } as unknown as Table)
    : ({ form, start: { length: rows }, end: { length: rows } } as unknown as Table);

describe('the data bounds', () => {
  it('default to int32 columns and 2^24 values a block', () => {
    expect(DefaultMaxColumnLength).toBe(2 ** 31 - 1);
    expect(DefaultMaxBlockValues).toBe(16_777_216);
    expect(DefaultAnswerBounds).toEqual({ maxColumnLength: 2 ** 31 - 1, maxBlockValues: 2 ** 24 });
  });

  it('pass a block of 2^24 values and refuse one more', () => {
    expect(checkTableBounds(claims('positions', 2 ** 24))).toBeUndefined();
    expect(checkTableBounds(claims('positions', 2 ** 24 + 1))).toMatch(/block of more than 16777216/);
    // a list's values are its rows times its columns (start and end)
    expect(checkTableBounds(claims('list', 2 ** 23))).toBeUndefined();
    expect(checkTableBounds(claims('list', 2 ** 23 + 1))).toMatch(/block/);
  });

  it('pass a column of the bound’s length and refuse one more, with the block bound raised past it', () => {
    const bounds = { maxColumnLength: 100, maxBlockValues: 1e9 };
    expect(checkTableBounds(claims('positions', 100), bounds)).toBeUndefined();
    expect(checkTableBounds(claims('positions', 101), bounds)).toMatch(/column of more than 100/);
  });

  it('are checked on an answer before its extent is walked', () => {
    const request: RegionRequest = { system: 'chr1', start: 0, end: 1000, level: 1, shapes: ['ranges'] };
    const column = (n: number) => {
      const a = new Int32Array(1);
      Object.defineProperty(a, 'length', { value: n });
      return a;
    };
    // a column that claims 2^31 rows: the extent loop would never end, the bound answers first
    const answer = { ok: true, tables: { ranges: { form: 'list', start: column(2 ** 31), end: column(2 ** 31) } } };
    expect(checkAnswer(answer, request, Number.MAX_SAFE_INTEGER)).toMatch(/column of more than 2147483647/);
    const small = checkAnswer(rangesAnswer(request), request, 1_000_000, { maxColumnLength: 5, maxBlockValues: 1e9 });
    expect(small).toMatch(/ranges table: a column of more than 5/);
    expect(checkAnswer(rangesAnswer(request), request, 1_000_000, { maxColumnLength: 10, maxBlockValues: 1e9 })).toBe(
      undefined,
    );
  });

  it('refuse a file that declares a longer system than a column can be', () => {
    const coverage = (length: number) => ({
      systems: [{ name: 's', length }],
      shapes: [{ system: 's', shape: 'letters' }],
      levels: [1],
    });
    const ok = readSeq({ coverage: coverage(2 ** 31 - 1), pairs: [] } as never);
    expect(ok.issues).toEqual([]);
    const over = readSeq({ coverage: coverage(2 ** 31), pairs: [] } as never);
    expect(over.issues.map((i) => i.code)).toEqual(['seq.source.malformed']);
    expect(over.issues[0].path).toEqual(['coverage', 'systems', 0, 'length']);
    const text = JSON.stringify({ kind: 'seq', params: { coverage: coverage(2 ** 40), pairs: [] } });
    expect(decodeSeq(text).issues.map((i) => i.code)).toEqual(['seq.source.malformed']);
  });
});

describe('the store holds the bounds a composition sets', () => {
  it('refuses an answer past its block bound as malformed, nothing kept', async () => {
    const clock = new FakeClock();
    const source = new SyntheticSource(syntheticCoverage(), clock, { delayMs: 10 });
    const store = new DataStore({ clock, settings: { maxBlockValues: 3 } });
    const S = '/data/genes';
    store.attach(S, source);
    store.want([{ source: S, system: 'chr1', start: 0, end: 1024, level: 1 }]);
    await clock.advance(10);
    expect(store.get(dataKey(S, 'chr1', 1, 0))).toMatchObject({
      status: 'failed',
      kind: 'malformed',
      issue: { code: 'seq.source.malformed' },
    });
  });
});
