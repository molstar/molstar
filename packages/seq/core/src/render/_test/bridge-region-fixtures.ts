/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge over a source's pairs read by region (generated): alignments of blocks between two chromosomes `A` and `B`
// (1 Mb each, levels 1, 1 kb, 16 kb), each alignment's blocks one `row`, and two sections with a bridge drawing them.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { PairsList, SeqData } from '@molstar/seqx-core/data';
import { binPairs, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';

export const AlnLength = 1_000_000;
const Levels = [1, 1024, 16_384];

/** An alignment: where it starts on each side, its strand on `B`, its blocks as `[length, gap on A, gap on B]`. */
export interface Alignment {
  readonly a: number;
  readonly b: number;
  readonly strand?: 1 | -1;
  readonly blocks: readonly (readonly [number, number, number])[];
}

/** The alignments' blocks as a pairs list (a block a row, an alignment's blocks one `row`), in order. */
export function alignmentList(alignments: readonly Alignment[]): PairsList {
  const as: number[] = [];
  const ae: number[] = [];
  const bs: number[] = [];
  const be: number[] = [];
  const strand: number[] = [];
  const row: number[] = [];
  alignments.forEach((al, r) => {
    let a = al.a;
    let b = al.b;
    for (const [len, ga, gb] of al.blocks) {
      as.push(a);
      ae.push(a + len);
      bs.push(b);
      be.push(b + len);
      strand.push(al.strand ?? 1);
      row.push(r);
      a += len + ga;
      b += len + gb;
    }
  });
  const n = as.length;
  return {
    form: 'list',
    aSystem: categoryColumn(Array.from({ length: n }, () => 'A')),
    bSystem: categoryColumn(Array.from({ length: n }, () => 'B')),
    aStart: Int32Array.from(as),
    aEnd: Int32Array.from(ae),
    bStart: Int32Array.from(bs),
    bEnd: Int32Array.from(be),
    aStrand: new Int8Array(n).fill(1),
    bStrand: Int8Array.from(strand),
    row: Int32Array.from(row),
  };
}

/** `count` alignments of 2 to 4 blocks (200 to 600 long, gaps under 100), both strands, seeded. */
export function randomAlignments(count: number, seed: number): Alignment[] {
  const r = new PCG(seed);
  const out: Alignment[] = [];
  for (let i = 0; i < count; i++) {
    const blocks: [number, number, number][] = [];
    const m = 2 + (r.int() % 3);
    for (let j = 0; j < m; j++) blocks.push([200 + (r.int() % 400), 1 + (r.int() % 99), 1 + (r.int() % 99)]);
    out.push({
      a: r.int() % (AlnLength - 5000),
      b: r.int() % (AlnLength - 5000),
      strand: r.int() % 3 ? 1 : -1,
      blocks,
    });
  }
  return out;
}

/** The alignments served by region from either end: the whole list a system at each level (bins above 1). */
export function alignmentSource(alignments: readonly Alignment[]): SeqData {
  const list = alignmentList(alignments);
  const whole = (system: string, level: number) => ({
    system,
    start: 0,
    end: AlnLength,
    level,
    shapes: ['pairs' as const],
  });
  return {
    coverage: {
      systems: [
        { name: 'A', length: AlnLength },
        { name: 'B', length: AlnLength },
      ],
      shapes: [
        { system: 'A', shape: 'pairs' },
        { system: 'B', shape: 'pairs' },
      ],
      levels: Levels,
    },
    pairs: ['A', 'B'].flatMap((system) =>
      Levels.map((level) => ({
        request: whole(system, level),
        answer: { ok: true as const, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
      })),
    ),
  };
}

/**
 * Two sections (`top` of A, `bottom` of B) and a bridge between them drawing the alignments read by region; with
 * `third`, a third section (`third`, of A) below them, joined to nothing.
 */
export function regionBridgeSpec(
  alignments: readonly Alignment[],
  options: { third?: boolean; pieces?: number; format?: string; mirror?: boolean } = {},
) {
  const b = createSeqSpec({ title: 'alignments across' });
  // `format`: a source of a runtime of the test's (`uri` none), the alignments its to serve
  const src = options.format
    ? b.data.source('src', { format: options.format, uri: 'test' } as never)
    : b.data.source('src', { format: 'seq', value: seqToValue(alignmentSource(alignments)) as never });
  const A = b.data.sequence('A', { alphabet: 'dna', source: src, name: 'A' });
  const B = b.data.sequence('B', { alphabet: 'dna', source: src, name: 'B' });
  const aln = b.data.pairs('aln', { a: A, b: B, source: src } as never);
  // `pieces`: the top shows every other position of A's first `2 · pieces` (single-position pieces, EC-33)
  const n = options.pieces ?? 0;
  const pieces =
    n > 0
      ? b.data.ranges('pieces', {
          target: A,
          rows: {
            start: Int32Array.from({ length: n }, (_, i) => 2 * i),
            end: Int32Array.from({ length: n }, (_, i) => 2 * i + 1),
          },
        } as never)
      : undefined;
  b.views.view('v', (v) => {
    v.section('top', { of: A, ...(pieces ? { pieces } : {}) } as never)
      .track('ruler')
      .layer('ruler', { type: 'ruler' });
    v.bridge('br', { a: '/views/v/top', b: '/views/v/bottom' } as never, (br) =>
      br.layer('aln', { type: 'pairs', data: aln } as never),
    );
    v.section('bottom', { of: B }).track('ruler').layer('ruler', { type: 'ruler' });
    if (options.third) v.section('third', { of: A }).track('ruler').layer('ruler', { type: 'ruler' });
  });
  return b.build();
}
