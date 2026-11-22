/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Bridge fixtures (seeded, generated): two sections of two sequences, or of two sets of blocks (forward, reversed,
// scaled), a bridge between them with a pairs layer of inline pairs or a segments map.

import type { Node } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { dna } from './genome-fixtures.js';

/** A pair as written: its two ranges and its `b` end's strand. */
export interface BridgePair {
  readonly a: readonly [number, number];
  readonly b: readonly [number, number];
  readonly strand?: '+' | '-';
  /** Its systems: `x` (the top's) and `y` (the bottom's) by default. */
  readonly on?: readonly ['x' | 'y' | 'z', 'x' | 'y' | 'z'];
}

/** A block of a side: its range, and whether it is reversed or scaled. */
export interface SideBlock {
  readonly range: readonly [number, number];
  readonly reversed?: boolean;
  readonly scale?: number;
}

export interface BridgeSpecOptions {
  readonly pairs: readonly BridgePair[];
  /** The top section's and the bottom section's blocks (absent: the whole sequence, no blocks). */
  readonly top?: readonly SideBlock[];
  readonly bottom?: readonly SideBlock[];
  /** A side's direction (`reversed`). */
  readonly topReversed?: boolean;
  readonly bottomReversed?: boolean;
  readonly bridge?: Record<string, unknown>;
  readonly layer?: Record<string, unknown>;
  /** The bridge's `a` side is the bottom section (its sides written the other way round). */
  readonly swapped?: boolean;
  /** Lengths of x and y. */
  readonly lengths?: readonly [number, number];
  /** The pairs written as a segments map's runs from x to y (each its `a` length, `b` its image), the layer's data. */
  readonly map?: boolean;
  /** One section only (`top`, of x), with a section of x below it too: one genome against itself. */
  readonly self?: boolean;
}

/** Two sections (`top` of x, `bottom` of y) and a bridge between them drawing the pairs. */
export function bridgeSpec(o: BridgeSpecOptions): Node {
  const b = createSeqSpec({ title: 'a bridge' });
  const [lx, ly] = o.lengths ?? [400, 300];
  const x = b.data.sequence('x', { alphabet: 'dna', letters: dna(lx, 3), name: 'chrX' });
  const y = b.data.sequence('y', { alphabet: 'dna', letters: dna(ly, 5), name: 'chrY' });
  const z = b.data.sequence('z', { alphabet: 'dna', letters: dna(100, 7), name: 'chrZ' });
  const sys = { x, y, z };
  const rows = {
    a_start: o.pairs.map((p) => p.a[0]),
    a_end: o.pairs.map((p) => p.a[1]),
    b_start: o.pairs.map((p) => p.b[0]),
    b_end: o.pairs.map((p) => p.b[1]),
    b_strand: o.pairs.map((p) => p.strand ?? '+'),
    id: o.pairs.map((_, i) => `p${i}`),
  };
  const on = o.pairs[0]?.on ?? ['x', 'y'];
  const p = o.map
    ? b.data.map('runs', {
        from: x,
        to: y,
        method: 'segments',
        segments: {
          from: o.pairs.map((q) => q.a[0]),
          to: o.pairs.map((q) => (q.strand === '-' ? q.b[0] + (q.a[1] - q.a[0]) - 1 : q.b[0])),
          length: o.pairs.map((q) => q.a[1] - q.a[0]),
          strand: o.pairs.map((q) => q.strand ?? '+'),
        },
      })
    : b.data.pairs('p', { a: sys[on[0]], b: sys[on[1]], rows });
  b.views.view('v', (v) => {
    const side = (id: string, of: typeof x, blocks: readonly SideBlock[] | undefined, reversed: boolean | undefined) =>
      v.section(id, { ...(blocks ? {} : { of }), ...(reversed ? { direction: 'reversed' } : {}) } as never, (s) => {
        s.track('letters').layer('letters', { type: 'letters' });
        blocks?.forEach((bl, i) =>
          s.block(`b${i}`, {
            of,
            range: bl.range as [number, number],
            ...(bl.reversed ? { direction: 'reversed' } : {}),
            ...(bl.scale ? { scale: bl.scale } : {}),
          } as never),
        );
      });
    side('top', x, o.top, o.topReversed);
    const ends = o.swapped ? { a: '/views/v/bottom', b: '/views/v/top' } : { a: '/views/v/top', b: '/views/v/bottom' };
    v.bridge('br', { ...ends, ...o.bridge } as never, (br) =>
      br.layer('ribbons', { type: 'pairs', data: p, ...o.layer } as never),
    );
    side('bottom', o.self ? x : y, o.bottom, o.bottomReversed);
  });
  return b.build();
}

/** `count` random pairs between two sequences of `N` (20,000), of both strands. */
export function randomBridgePairs(count: number, seed: number, N = 20_000): BridgePair[] {
  const r = new PCG(seed);
  const out: BridgePair[] = [];
  for (let i = 0; i < count; i++) {
    const a = r.int() % (N - 200);
    const b = r.int() % (N - 200);
    out.push({ a: [a, a + 1 + (r.int() % 150)], b: [b, b + 1 + (r.int() % 150)], strand: r.int() % 4 ? '+' : '-' });
  }
  return out;
}
