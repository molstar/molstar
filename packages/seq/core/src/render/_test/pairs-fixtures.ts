/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Specs of pairs and grids for the render tests: a sequence with pairs over it as arcs (tracks, a ring, blocks of two
// chains), and a dot plot (a grid section following a tracks section) with a grid of values and the pairs as dots.

import type { Node } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec } from '@molstar/seqx-core/lang';

export const Letters = 'ACGU';

export function letters(n: number, seed = 7): string {
  const r = new PCG(seed);
  let s = '';
  for (let i = 0; i < n; i++) s += Letters[r.int() % 4];
  return s;
}

/** Pairs of single positions `[i, j]` as a pairs node's rows, each `type` its kind. */
export function pairRows(pairs: readonly (readonly [number, number])[], types?: readonly string[]) {
  return {
    a_start: pairs.map((p) => p[0]),
    a_end: pairs.map((p) => p[0] + 1),
    b_start: pairs.map((p) => p[1]),
    b_end: pairs.map((p) => p[1] + 1),
    id: pairs.map((p) => `${p[0]}-${p[1]}`),
    ...(types ? { type: types } : {}),
  };
}

/** Random pairs of distinct positions in `[0, n)`, each position once at most (a secondary structure's shape). */
export function randomPairs(n: number, count: number, seed: number): [number, number][] {
  const r = new PCG(seed);
  const used = new Set<number>();
  const out: [number, number][] = [];
  for (let t = 0; out.length < count && t < 50 * count; t++) {
    const a = r.int() % n;
    const b = r.int() % n;
    if (a === b || used.has(a) || used.has(b)) continue;
    used.add(a);
    used.add(b);
    out.push([Math.min(a, b), Math.max(a, b)]);
  }
  return out;
}

export interface ArcsOptions {
  readonly n?: number;
  readonly pairs: readonly (readonly [number, number])[];
  readonly types?: readonly string[];
  readonly crossing?: 'same' | 'split';
  readonly side?: 'above' | 'below';
  readonly layout?: 'rings' | 'text';
  readonly height?: number;
}

/** A sequence, its letters and its pairs as arcs in one section. */
export function arcsSpec(o: ArcsOptions): Node {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'rna', letters: letters(o.n ?? 76) });
  const p = b.data.pairs('p', { a: s, rows: pairRows(o.pairs, o.types) });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s, ...(o.layout ? { layout: o.layout } : {}) });
    sec.track('letters').layer('letters', { type: 'letters' });
    sec.track('pairs', { ...(o.height ? { height: o.height } : {}) }).layer('arcs', {
      type: 'pairs',
      data: p,
      ...(o.crossing ? { crossing: o.crossing } : {}),
      ...(o.side ? { side: o.side } : {}),
    });
  });
  return b.build();
}

export interface GridOptions {
  readonly n: number;
  readonly values?: Float32Array;
  readonly pairs?: readonly (readonly [number, number])[];
  readonly mirror?: boolean;
  readonly triangle?: 'both' | 'upper' | 'lower';
  readonly aggregate?: 'mean' | 'max' | 'min';
  /** The grid's y another system of `ny` positions (its own window: `follow_y: none`). */
  readonly ny?: number;
}

/** A value per cell from its two positions: `i * 1000 + j` (each cell its own, exact in float32 up to 16,000). */
export function cellValues(nx: number, ny: number): Float32Array {
  const v = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) v[i * ny + j] = i * 1000 + j;
  return v;
}

/** A tracks section of the sequence, and a grid section following it with a grid layer and the pairs as dots. */
export function gridSpec(o: GridOptions): Node {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'rna', letters: letters(o.n) });
  const y = o.ny ? b.data.sequence('y', { alphabet: 'rna', letters: letters(o.ny, 9) }) : undefined;
  const g = b.data.grid('g', { x: s, ...(y ? { y } : {}), values: o.values ?? cellValues(o.n, o.ny ?? o.n) });
  const p = o.pairs ? b.data.pairs('p', { a: s, rows: pairRows(o.pairs) }) : undefined;
  b.views.view('v', (v) => {
    const top = v.section('s', { of: s });
    top.track('letters').layer('letters', { type: 'letters' });
    v.section('plot', {
      layout: 'grid',
      of: s,
      follow: top.address,
      ...(y ? { of_y: y, follow_y: 'none' } : {}),
    }).track('plot', {}, (t) => {
      t.layer('grid', {
        type: 'grid',
        data: g,
        ...(o.triangle ? { triangle: o.triangle } : {}),
        ...(o.aggregate ? { aggregate: o.aggregate } : {}),
      });
      if (p)
        t.layer('dots', {
          type: 'pairs',
          data: p,
          ...(o.mirror ? { mirror: true } : {}),
          ...(o.triangle ? { triangle: o.triangle } : {}),
        });
    });
  });
  return b.build();
}
