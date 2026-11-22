/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover (EC-212): an old and a new assembly named alike (`chr1`, `chr2` on both: the map's source is sided), the
// map's runs (a straight stretch, a deletion, an inversion, a duplication and a translocation) served by region from
// either side, the old assembly's genes inline, and the figure: the old section, a bridge over the map, the new
// section with the genes drawn through it. The truths are the runs' arithmetic, computed here.

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import type { PairsList, SeqData, SeqPair } from '@molstar/seqx-core/data';
import { binPairs, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';

export const OldChromosomes = [
  { name: 'chr1', length: 100_000 },
  { name: 'chr2', length: 60_000 },
];
export const NewChromosomes = [
  { name: 'chr1', length: 80_000 },
  { name: 'chr2', length: 80_000 },
];

/** A run: `[from, from + length)` on old `a` to new `b` from `to`, read backwards when `strand` is −1. */
export interface LiftRun {
  readonly a: string;
  readonly from: number;
  readonly b: string;
  readonly to: number;
  readonly length: number;
  readonly strand: 1 | -1;
}

export const Runs: readonly LiftRun[] = [
  { a: 'chr1', from: 0, b: 'chr1', to: 0, length: 20_000, strand: 1 },
  // old chr1 20–30 kb deleted (no run)
  { a: 'chr1', from: 30_000, b: 'chr1', to: 20_000, length: 20_000, strand: 1 },
  // an inversion
  { a: 'chr1', from: 50_000, b: 'chr1', to: 40_000, length: 20_000, strand: -1 },
  // a duplication: one old stretch in two places
  { a: 'chr1', from: 70_000, b: 'chr1', to: 60_000, length: 20_000, strand: 1 },
  { a: 'chr1', from: 70_000, b: 'chr2', to: 60_000, length: 20_000, strand: 1 },
  // a translocation: old chr1's tail at the head of new chr2
  { a: 'chr1', from: 90_000, b: 'chr2', to: 0, length: 10_000, strand: 1 },
  { a: 'chr2', from: 0, b: 'chr2', to: 10_000, length: 50_000, strand: 1 },
];

/** The old genes: one straight, one half deleted, one inverted, one duplicated, one translocated, one lost. */
export const Genes = {
  id: ['straight', 'half', 'inverted', 'twice', 'moved', 'other', 'lost'],
  target: ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr1'],
  start: [5_000, 25_000, 55_000, 72_000, 92_000, 20_000, 22_000],
  end: [15_000, 35_000, 60_000, 78_000, 95_000, 30_000, 28_000],
};

/** The images of old `[s, e)` on `a`: each run's part of it, on the new assembly, sorted. */
export function imagesOf(a: string, s: number, e: number): { b: string; start: number; end: number }[] {
  const out: { b: string; start: number; end: number }[] = [];
  for (const r of Runs) {
    if (r.a !== a) continue;
    const lo = Math.max(s, r.from);
    const hi = Math.min(e, r.from + r.length);
    if (hi <= lo) continue;
    if (r.strand === 1) out.push({ b: r.b, start: r.to + lo - r.from, end: r.to + hi - r.from });
    else out.push({ b: r.b, start: r.to + r.length - (hi - r.from), end: r.to + r.length - (lo - r.from) });
  }
  return out.sort((x, y) => x.b.localeCompare(y.b) || x.start - y.start);
}

/** The runs as a pairs list: `a` on the old side, `b` on the new, a row each. */
export function runsList(runs: readonly LiftRun[] = Runs): PairsList {
  const n = runs.length;
  return {
    form: 'list',
    aSystem: categoryColumn(runs.map((r) => r.a)),
    bSystem: categoryColumn(runs.map((r) => r.b)),
    aStart: Int32Array.from(runs, (r) => r.from),
    aEnd: Int32Array.from(runs, (r) => r.from + r.length),
    bStart: Int32Array.from(runs, (r) => r.to),
    bEnd: Int32Array.from(runs, (r) => r.to + r.length),
    aStrand: new Int8Array(n).fill(1),
    bStrand: Int8Array.from(runs, (r) => r.strand),
    row: Int32Array.from(runs, (_, i) => i),
  };
}

const Levels = [1, 1024];

/** The runs served by region on each side (side `a` the old assembly's names, `b` the new's), bins above 1. */
export function liftSource(runs: readonly LiftRun[] = Runs): SeqData {
  const list = runsList(runs);
  const pairs: SeqPair[] = [];
  for (const [side, chrs] of [
    ['a', OldChromosomes],
    ['b', NewChromosomes],
  ] as const)
    for (const c of chrs)
      for (const level of Levels)
        pairs.push({
          request: { system: c.name, start: 0, end: c.length, level, shapes: ['pairs'], side },
          answer: { ok: true, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        });
  return {
    coverage: {
      systems: [
        ...OldChromosomes.map((c) => ({ ...c, side: 'a' as const })),
        ...NewChromosomes.map((c) => ({ ...c, side: 'b' as const })),
      ],
      shapes: ['chr1', 'chr2'].map((system) => ({ system, shape: 'pairs' as const })),
      levels: Levels,
    },
    pairs,
  };
}

export const OldSection = '/views/v/old';
export const NewSection = '/views/v/new';
export const LiftMap = '/data/lift';
export const GenesData = '/data/genes';

/** The figure: the old section with its genes, a bridge over the map, the new section with the genes lifted. */
export function liftSpec(runs: readonly LiftRun[] = Runs): Node {
  const b = createSeqSpec({ title: 'a liftover' });
  const old = b.data.assembly('old', assembly('Old', OldChromosomes, { order: 'source' } as never));
  const neu = b.data.assembly('new', assembly('New', NewChromosomes, { order: 'source' } as never));
  const src = b.data.source('chain', { format: 'seq', value: seqToValue(liftSource(runs)) as never });
  const map = b.data.map('lift', { from: old, to: neu, method: 'segments', source: src } as never);
  const genes = b.data.ranges('genes', { target: old, rows: Genes } as never);
  b.views.view('v', (v) => {
    v.section('old', { of: old } as never, (s) => s.track('genes').layer('genes', { type: 'ranges', data: genes }));
    v.bridge('br', { a: OldSection, b: NewSection } as never, (br) =>
      br.layer('runs', { type: 'pairs', data: map } as never),
    );
    v.section('new', { of: neu } as never, (s) => s.track('genes').layer('genes', { type: 'ranges', data: genes }));
  });
  return b.build();
}
