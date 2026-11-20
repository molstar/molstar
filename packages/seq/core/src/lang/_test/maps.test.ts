/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Language 0.5.0's map params (_spec/README.md, Language): a map's `segments` and `structure` methods and its
// `priority`, a section's `via`, a numbering's `offset` and `segments` methods; what each method takes, checked.

import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';

const codes = (n: Node) => checkSeqSpec(n).map((i) => i.code);

function mapSpec(params: Record<string, unknown>): Node {
  const b = createSeqSpec();
  const a = b.data.sequence('a', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
  const t = b.data.sequence('t', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
  b.data.map('m', { from: a, to: t, ...params } as never);
  return b.build();
}

describe('maps: methods and their params', () => {
  it('segments and structure take their runs; structure its provenance and status; any map a priority', () => {
    const segments = { from: [0, 5], to: [0, 6], length: [4, 3], strand: ['+', '-'] };
    expect(codes(mapSpec({ method: 'segments', segments, priority: -2 }))).toEqual([]);
    const structure = { ...segments, provenance: ['sifts', 'alignment'] };
    const status = { start: [1], end: [2], status: ['conflict'] };
    expect(codes(mapSpec({ method: 'structure', segments: structure, status, priority: 3 }))).toEqual([]);
    expect(codes(mapSpec({ method: 'offset', by: 1, priority: 1 }))).toEqual([]);
  });

  it('params another method takes, or a missing one, are `seq.lang.map.params`', () => {
    const segments = { from: [0], to: [0], length: [4] };
    for (const params of [
      { method: 'segments' },
      { method: 'segments', segments, by: 1 },
      { method: 'segments', segments, status: { start: [1], end: [2], status: ['conflict'] } },
      { method: 'segments', segments: { ...segments, provenance: ['sifts'] } },
      { method: 'structure', segments, strand: '+' },
      { method: 'offset', by: 1, segments },
    ])
      expect(codes(mapSpec(params))).toEqual(['seq.lang.map.params']);
  });

  it('runs that overlap on `from`, hold nothing or read no strand, and a status that is none, are `seq.lang.map.segments`', () => {
    for (const segments of [
      { from: [0, 3], to: [0, 5], length: [4, 2] },
      { from: [0], to: [0], length: [0] },
      { from: [0], to: [0], length: [2], strand: ['x'] },
    ])
      expect(codes(mapSpec({ method: 'segments', segments }))).toEqual(['seq.lang.map.segments']);
    const status = { start: [3], end: [2], status: ['conflict'] };
    const runs = { from: [0], to: [0], length: [4] };
    expect(codes(mapSpec({ method: 'structure', segments: runs, status }))).toEqual(['seq.lang.map.segments']);
    const word = { start: [1], end: [2], status: ['paired'] };
    expect(codes(mapSpec({ method: 'structure', segments: runs, status: word }))).toEqual(['seq.lang.map.segments']);
  });
});

describe('a section’s via', () => {
  function viaSpec(via: (m: { a: string; b: string; s: string }) => unknown): Node {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEP' });
    const t = b.data.sequence('t', { alphabet: 'protein', letters: 'MEEP' });
    const a = b.data.map('a', { from: s, to: t, method: 'offset', by: 0 });
    const m = b.data.map('b', { from: s, to: t, method: 'offset', by: 1 });
    b.views.view('v', (v) => {
      v.section('sec', { of: t, via: via({ a, b: m, s }) as never });
    });
    return b.build();
  }

  it('names a map, or a list of maps', () => {
    expect(codes(viaSpec(({ a }) => a))).toEqual([]);
    expect(codes(viaSpec(({ a, b }) => [a, b]))).toEqual([]);
  });

  it('names nothing but maps (`seq.lang.ref.target`, at the item)', () => {
    const [one] = checkSeqSpec(viaSpec(({ s }) => s));
    expect([one.code, one.path?.slice(-1)]).toEqual(['seq.lang.ref.target', ['via']]);
    const [item] = checkSeqSpec(viaSpec(({ a, s }) => [a, s]));
    expect([item.code, item.path?.slice(-2)]).toEqual(['seq.lang.ref.target', ['via', 1]]);
  });
});

describe('numberings: offset and segments', () => {
  function numbered(params: Record<string, unknown>): Node {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'protein', letters: 'MEEPQS' }, (s) => s.numbering('auth', params as never));
    return b.build();
  }

  it('an offset takes `by`; segments their runs, an insertion code on a run of one', () => {
    expect(codes(numbered({ method: 'offset', by: 94 }))).toEqual([]);
    const segments = { start: [0, 3, 4], end: [3, 4, 6], first: [1, 3, 4], insertion: [null, 'A', null] };
    expect(codes(numbered({ method: 'segments', segments }))).toEqual([]);
  });

  it('params that do not fit the method are `seq.lang.numbering.content`; bad runs `seq.lang.numbering.segments`', () => {
    for (const params of [
      { method: 'offset' },
      { method: 'offset', by: 1, labels: ['1', '2', '3', '4', '5', '6'] },
      { method: 'segments' },
      { method: 'labels', by: 1, labels: ['1', '2', '3', '4', '5', '6'] },
    ])
      expect(codes(numbered(params))).toEqual(['seq.lang.numbering.content']);
    for (const segments of [
      { start: [3], end: [1], first: [1] },
      { start: [0, 2], end: [3, 5], first: [1, 10] },
      { start: [0], end: [2], first: [1], insertion: ['A'] },
    ])
      expect(codes(numbered({ method: 'segments', segments }))).toEqual(['seq.lang.numbering.segments']);
  });
});
