/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';

// The tRNA dataset: yeast tRNA-Phe as in 1EHZ (76 nucleotides, modified bases written as their parents), its base
// pairs from the structure (the four stems, nested, and the tertiary pairs that cross them: 8–14, 15–48 the Levitt
// pair, 18–55, 19–56, 26–44, 54–58), each of a kind (canonical, wobble, non-canonical); and a base-pair probability
// grid on the same 76 positions: generated from the structure's pairs and the sequence's complementarity with a seed
// (a stand-in until a computed covariation is recorded), so a grid layer draws on its own systems only.

export const title = 'tRNA-Phe';

/** Yeast tRNA-Phe (1EHZ), 5′ to 3′. */
export const TrnaLetters = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA';

/** The stems' pairs, 1-based as the structure numbers them. */
const Stems: readonly (readonly [number, number])[] = [
  [1, 72],
  [2, 71],
  [3, 70],
  [4, 69],
  [5, 68],
  [6, 67],
  [7, 66],
  [10, 25],
  [11, 24],
  [12, 23],
  [13, 22],
  [27, 43],
  [28, 42],
  [29, 41],
  [30, 40],
  [31, 39],
  [49, 65],
  [50, 64],
  [51, 63],
  [52, 62],
  [53, 61],
];

/** The tertiary pairs, 1-based, and their kinds: they cross the stems (pseudoknot-like in an arc diagram). */
const Tertiary: readonly (readonly [number, number, string])[] = [
  [8, 14, 'non-canonical'],
  [15, 48, 'non-canonical'],
  [18, 55, 'non-canonical'],
  [19, 56, 'canonical'],
  [26, 44, 'non-canonical'],
  [54, 58, 'non-canonical'],
];

const Canonical = new Set(['GC', 'CG', 'AU', 'UA']);
const Wobble = new Set(['GU', 'UG']);

/** A pair's kind by its letters: Watson–Crick, G·U, else non-canonical. */
function kindOf(i: number, j: number): string {
  const p = TrnaLetters[i] + TrnaLetters[j];
  return Canonical.has(p) ? 'canonical' : Wobble.has(p) ? 'wobble' : 'non-canonical';
}

/** The base pairs, 0-based and half-open: a row a pair, its kind, its id (`i-j`, 1-based). */
export function trnaPairs() {
  const all = [
    ...Stems.map(([i, j]) => [i - 1, j - 1, kindOf(i - 1, j - 1)] as const),
    ...Tertiary.map(([i, j, k]) => [i - 1, j - 1, k] as const),
  ];
  return {
    a_start: all.map((p) => p[0]),
    a_end: all.map((p) => p[0] + 1),
    b_start: all.map((p) => p[1]),
    b_end: all.map((p) => p[1] + 1),
    type: all.map((p) => p[2]),
    id: all.map((p) => `${p[0] + 1}-${p[1] + 1}`),
  };
}

/**
 * The base-pair probabilities, 76 × 76, x-major and symmetric: 0.95 at a stem's pair, 0.5 at a tertiary one, a small
 * value from a seed where two bases could pair (Watson–Crick or G·U, at least four apart), missing elsewhere.
 */
export function trnaProbabilities(): Float32Array {
  const n = TrnaLetters.length;
  const v = new Float32Array(n * n).fill(Number.NaN);
  const r = new PCG(76);
  for (let i = 0; i < n; i++)
    for (let j = i + 4; j < n; j++) {
      const p = TrnaLetters[i] + TrnaLetters[j];
      if (!Canonical.has(p) && !Wobble.has(p)) continue;
      const x = 0.02 + (r.int() % 1000) / 10000;
      v[i * n + j] = v[j * n + i] = x;
    }
  for (const [i, j] of Stems) v[(i - 1) * n + (j - 1)] = v[(j - 1) * n + (i - 1)] = 0.95;
  for (const [i, j] of Tertiary) v[(i - 1) * n + (j - 1)] = v[(j - 1) * n + (i - 1)] = 0.5;
  return v;
}

export interface TrnaOptions {
  readonly crossing?: 'same' | 'split';
  readonly layout?: 'rings';
  readonly dotPlot?: boolean;
}

/** The tRNA's sequence, its pairs and its probabilities, and a view of them. */
export function build(options: TrnaOptions = {}): Node {
  const b = createSeqSpec({ title });
  const trna = b.data.sequence('trna', { alphabet: 'rna', letters: TrnaLetters, name: 'tRNA-Phe' });
  const pairs = b.data.pairs('pairs', { a: trna, rows: trnaPairs() });
  const probability = b.data.grid('probability', { x: trna, values: trnaProbabilities() });
  b.views.view('main', (v) => {
    const s = v.section('trna', { of: trna, ...(options.layout ? { layout: options.layout } : {}) });
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    s.track('letters', { label: 'tRNA-Phe' }).layer('letters', { type: 'letters' });
    s.track('pairs', { label: 'Base pairs' }).layer('arcs', {
      type: 'pairs',
      data: pairs,
      ...(options.crossing ? { crossing: options.crossing } : {}),
      ...(options.layout === 'rings' ? { side: 'below' } : {}),
    });
    if (!options.dotPlot) return;
    // the dot plot: x and y the tRNA's, y following x; the probabilities above the diagonal, the structure's pairs
    // below it (mirrored: the builder sets `mirror`, a plain `false` otherwise)
    v.section('dots', { layout: 'grid', of: trna, follow: s.address }).track('plot', { label: 'Dot plot' }, (t) => {
      t.layer('probability', {
        type: 'grid',
        data: probability,
        triangle: 'upper',
        color: theme.field('value', { domain: [0, 1], palette: theme.continuous({ colors: 'blues' }) }),
      });
      t.layer('structure', { type: 'pairs', data: pairs, mirror: true, triangle: 'lower' });
    });
  });
  return b.build();
}

/** Crossing pairs on the other side: the largest nested set above, the tertiary pairs that cross it below. */
export const split = () => build({ crossing: 'split' });

/** The tRNA on a ring: its pairs as chords through the hole. */
export const ring = () => build({ layout: 'rings' });

/** The arcs over the letters, and the base-pair dot plot under them following their window. */
export const dotPlot = () => build({ dotPlot: true });
