/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A dot plot's documents, for its test and heap guard: a sequence with a tandem repeat against itself, its word
// matches as pairs and a grid (|i − j| per pair of positions, missing on the diagonal) on one track, under a ruler.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

/** A repeat of 24, a linker, the repeat again: words of the repeat match at an offset of 30. */
export const REPEAT = 'MKTAYIAKQRQISFVKSHFSRQLE';
export const SEQ = `${REPEAT}GSGSGS${REPEAT}`;

/** Self matches of `k` letters (i ≠ j), as pairs of ranges merged along each diagonal. */
export function wordMatches(seq: string, k = 4) {
  const cols = { target: [] as string[], start: [] as number[], stop: [] as number[], b_start: [] as number[] };
  const b_stop: number[] = [];
  const n = seq.length;
  const word = (i: number, d: number) =>
    i >= 0 && i + d >= 0 && i + k <= n && i + d + k <= n && seq.slice(i, i + k) === seq.slice(i + d, i + d + k);
  for (let d = 1 - n; d < n; d++) {
    if (d === 0) continue;
    for (let i = Math.max(0, -d); i + k <= n && i + d + k <= n; i++) {
      if (!word(i, d) || word(i - 1, d)) continue;
      let e = i;
      while (word(e + 1, d)) e++;
      cols.target.push('s');
      cols.start.push(i);
      cols.stop.push(e + k);
      cols.b_start.push(i + d);
      b_stop.push(e + k + d);
    }
  }
  return { ...cols, b_target: cols.target, b_stop };
}

/** |i − j| per pair of positions of s, missing on the diagonal. */
export function distances(n: number) {
  const value = Array.from({ length: n * n }, (_, x) => {
    const [i, j] = [Math.floor(x / n), x % n];
    return i === j ? null : Math.abs(i - j);
  });
  return { target: ['s'], b_target: ['s'], value: [value] };
}

/** |i − j| per pair of positions of s as one Float32Array, NaN on the diagonal: a large grid without a JSON list. */
export function typedDistances(n: number): Float32Array {
  const value = new Float32Array(n * n);
  for (let x = 0; x < n * n; x++)
    value[x] = Math.floor(x / n) === x % n ? Number.NaN : Math.abs(Math.floor(x / n) - (x % n));
  return value;
}

/** The fixture; `typed`: its grid's cells one Float32Array (written into the document after the builder). */
export function dotPlotDoc(seq = SEQ, typed = false): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'dot plot' });
  const s = b.data.sequence({ ref: 's', value: seq });
  const cells = typed ? { target: ['s'], b_target: ['s'], value: [[]] } : distances(seq.length);
  const grid = b.data.annotation({ ref: 'dist', label: 'Distance', geometry: 'grid', value: cells });
  const words = b.data.annotation({ ref: 'words', geometry: 'pairs', value: wordMatches(seq) });
  const v = b.view({ ref: 'v', label: 'V' });
  const sec = v.section({ ref: 'v.s' }).block({ sequence: s });
  sec.track({ ref: 't-ruler', label: '' }).ruler();
  sec
    .track({ ref: 't-plot', label: 'Plot' })
    .dot_plot({ ref: 'l-dist', data: grid })
    .dot_plot({ ref: 'l-words', data: words, color: { theme: 'uniform', params: { color: '#d1495b' } } });
  const d = b.getDocument();
  if (typed) {
    const g = (d.root.children ?? []).flatMap((n) => n.children ?? []).find((n) => n.ref === 'dist');
    if (g?.params) (g.params.value as { value: unknown[] }).value = [typedDistances(seq.length)];
  }
  return d;
}
