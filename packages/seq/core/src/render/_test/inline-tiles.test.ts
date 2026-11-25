/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Inline data over a long system (_spec/README.md "Data", The seq format; "Render", Data): an inline ranges node is
// served at the level whose one tile covers its system too, so a frame at the fit of 83 Mb reads a tile of it, not
// 81,307 at level 1; a reader past the store's bound on a want's tiles builds no key a tile and is pending, never empty.

import { inlineSeq, openSeq, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it, vi } from 'vitest';
import type { RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const Length = 83_257_441;
const Starts = [1_000_000, 40_000_000, 80_000_000];

/** A chromosome-long system (its length a source's coverage) and an inline ranges node on it, three items far apart. */
function spec() {
  const b = createSeqSpec({ title: 'inline' });
  const value = seqToValue({
    coverage: { systems: [{ name: 'chr', length: Length }], shapes: [], levels: [1] },
    pairs: [],
  });
  const source = b.data.source('chr-source', { format: 'seq', value: value as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source, name: 'chr' });
  const genes = b.data.ranges('genes', {
    target: chr,
    rows: { start: Starts, end: Starts.map((s) => s + 20_000), label: ['a', 'b', 'c'] },
  });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr }).track('genes').layer('genes', { type: 'ranges', data: genes });
  });
  return b.build();
}

describe('inline data over a long system', () => {
  it('a ranges list is served at its whole level too, the same rows with one identity at both', async () => {
    const tables = {
      ranges: { form: 'list' as const, start: Int32Array.of(5, 2_000_000), end: Int32Array.of(9, 2_000_100) },
    };
    const data = inlineSeq({ name: 'chr', length: Length }, tables);
    expect(data.coverage.levels).toEqual([1, 16 ** 5]);
    const source = openSeq(data);
    const signal = new AbortController().signal;
    const read = async (level: number) => {
      const a = await source.read({ system: 'chr', start: 0, end: Length, level, shapes: ['ranges'] }, signal);
      return a.ok && a.tables.ranges?.form === 'list' ? Array.from(a.tables.ranges.row ?? []) : undefined;
    };
    expect(await read(16 ** 5)).toEqual(await read(1));
    // a short system's inline data stays one level
    expect(inlineSeq({ name: 's', length: 300 }, tables).coverage.levels).toEqual([1]);
  });

  it('a file’s list rows without ids keep their file-wide numbering at every level', async () => {
    // two levels of other rows: a row numbered after the rows of every pair before it, whatever its level
    const list = (n: number) => ({
      form: 'list' as const,
      start: Int32Array.from({ length: n }, (_, i) => 10 * i),
      end: Int32Array.from({ length: n }, (_, i) => 10 * i + 5),
    });
    const pair = (level: number, n: number) => ({
      request: { system: 'chr', start: 0, end: 100, level, shapes: ['ranges' as const] },
      answer: { ok: true as const, tables: { ranges: list(n) } },
    });
    const source = openSeq({
      coverage: {
        systems: [{ name: 'chr', length: 100 }],
        shapes: [{ system: 'chr', shape: 'ranges' }],
        levels: [1, 16],
      },
      pairs: [pair(1, 3), pair(16, 2)],
    });
    const signal = new AbortController().signal;
    const rows = async (level: number) => {
      const a = await source.read({ system: 'chr', start: 0, end: 100, level, shapes: ['ranges'] }, signal);
      return a.ok && a.tables.ranges?.form === 'list' ? Array.from(a.tables.ranges.row ?? []) : undefined;
    };
    expect(await rows(1)).toEqual([0, 1, 2]);
    expect(await rows(16)).toEqual([3, 4]);
  });

  it('a frame at the fit of 83 Mb reads a tile of an inline node, its items drawn', async () => {
    const s = await settle(spec(), [SeqRuntime]);
    const keysOf = vi.spyOn(s.store, 'keysOf');
    const f = await s.frame({ size: { width: 1440, height: 80 }, dpr: 2, viewport: { start: 0, end: Length } });
    let keys = 0;
    for (const r of keysOf.mock.results) keys += (r.value as string[]).length;
    expect(keys).toBeLessThan(64);
    expect(Math.max(...f.wants.map((w) => s.store.keysOf(w).length))).toBeLessThanOrEqual(2);
    expect(s.store.issues.value).toEqual([]);
    const lay = layerLay<RangesLay>(f, 0);
    expect(lay.pending).toEqual([]);
    expect(lay.items).toBe(Starts.length);
  });

  it('a view past a shorter source’s end is drawn empty, not pending; past the tile bound it is pending', async () => {
    // a source of ranges on the chromosome that knows it 2,000 positions long: past that, nothing to read
    const b = createSeqSpec({ title: 'short' });
    const value = seqToValue({
      coverage: { systems: [{ name: 'chr', length: Length }], shapes: [], levels: [1] },
      pairs: [],
    });
    const chr = b.data.sequence('chr', {
      alphabet: 'dna',
      source: b.data.source('chr-source', { format: 'seq', value: value as never }),
      name: 'chr',
    });
    const short = seqToValue(
      inlineSeq(
        { name: 'chr', length: 2000 },
        { ranges: { form: 'list', start: Int32Array.of(100), end: Int32Array.of(200) } },
      ),
    );
    const sites = b.data.ranges('sites', {
      target: chr,
      source: b.data.source('short', { format: 'seq', value: short as never }),
    });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('sites').layer('sites', { type: 'ranges', data: sites });
    });
    const s = await settle(b.build(), [SeqRuntime]);
    const past = await s.frame({ size: { width: 1440, height: 80 }, viewport: { start: 5000, end: 6000 } });
    expect(layerLay<RangesLay>(past, 0).pending).toEqual([]);
    // the store's bound is told apart: a want over it is pending
    const want = { source: '/data/short', system: 'chr', start: 0, end: 2000, level: 1 };
    expect(s.store.tooMany(want)).toBe(false);
    const big = { source: '/data/chr-source', system: 'chr', start: 0, end: Length, level: 1 };
    expect(s.store.tooMany(big)).toBe(true);
    expect(s.store.keysOf(big)).toEqual([]);
  });
});
