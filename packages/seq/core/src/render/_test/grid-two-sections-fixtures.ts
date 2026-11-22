/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid over two sections of blocks (_spec/README.md "Resolved view", A grid over blocks; EC-209): blocks `s0…` on a
// section `a`, blocks `t0…` on a section `b`, a grid a pair `(s_i, t_j)` (values `i·100 + j`, x-major), the plot
// following `a` across and `b` down; optionally pairs between two blocks as dots and diagonals.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { letters } from './pairs-fixtures.js';

export interface TwoSectionsOptions {
  /** The blocks' lengths on `a` (x) and on `b` (y). */
  readonly xs: readonly number[];
  readonly ys: readonly number[];
  /** `b`'s blocks reversed, by index (EC-145). */
  readonly reversed?: readonly number[];
  /** `b` a stack (refused). */
  readonly stack?: boolean;
  /** No grids: the pairs alone (a dot plot). */
  readonly noGrids?: boolean;
  /** Pairs from block `s[x]` to `t[y]`: each `[aStart, aEnd, bStart, bEnd, strand]` (strand −1: against). */
  readonly pairs?: { readonly x: number; readonly y: number; readonly rows: readonly (readonly number[])[] };
}

export function twoSectionsSpec(o: TwoSectionsOptions): Node {
  const b = createSeqSpec();
  const xs = o.xs.map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'dna', letters: letters(n, k + 1) }));
  const ys = o.ys.map((n, k) => b.data.sequence(`t${k}`, { alphabet: 'dna', letters: letters(n, k + 11) }));
  const grids: string[] = [];
  if (!o.noGrids)
    o.xs.forEach((nx, x) =>
      o.ys.forEach((ny, y) => {
        const v = new Float32Array(nx * ny);
        for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) v[i * ny + j] = x * 1000 + y * 100 + i * 10 + j;
        grids.push(b.data.grid(`g${x}${y}`, { x: xs[x], y: ys[y], values: v }));
      }),
    );
  const p = o.pairs;
  const pairs =
    p &&
    b.data.pairs('p', {
      a: xs[p.x],
      b: ys[p.y],
      rows: {
        a_start: p.rows.map((r) => r[0]),
        a_end: p.rows.map((r) => r[1]),
        b_start: p.rows.map((r) => r[2]),
        b_end: p.rows.map((r) => r[3]),
        a_strand: p.rows.map(() => '+'),
        b_strand: p.rows.map((r) => (r[4] === -1 ? '-' : '+')),
        id: p.rows.map((_, k) => `r${k}`),
      },
    });
  b.views.view('v', (v) => {
    const a = v.section('a', {});
    xs.forEach((q, k) => a.block(`b${k}`, { of: q, label: `s${k}` }));
    a.track('letters').layer('letters', { type: 'letters' });
    const bs = v.section('b', o.stack ? { arrange: 'stack' } : {});
    ys.forEach((q, k) =>
      bs.block(`b${k}`, { of: q, label: `t${k}`, ...(o.reversed?.includes(k) ? { direction: 'reversed' } : {}) }),
    );
    bs.track('letters').layer('letters', { type: 'letters' });
    const g = v.section('plot', { layout: 'grid', follow: a.address, follow_y: bs.address });
    g.track('plot', {}, (t) => {
      grids.forEach((d, k) => t.layer(`g${k}`, { type: 'grid', data: d }));
      if (pairs) t.layer('dots', { type: 'pairs', data: pairs });
    });
  });
  return b.build();
}
