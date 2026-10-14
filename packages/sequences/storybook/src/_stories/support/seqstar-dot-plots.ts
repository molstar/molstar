/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot story's document (Seq*/Views and layouts/Dot plots): human calmodulin against itself. Its two lobes
// are a duplication, each of two EF hands, so the plot shows parallel diagonals off the main one: the windowed
// identity as a grid (the share of equal residues over 9 positions centered on each pair), and the 4-letter words both
// copies share as segments over it. Computed here, as data: the format only draws them.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

/** Human calmodulin (UniProt P0DP23), the initiator methionine included. */
export const CALMODULIN =
  'MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK';
/** UniProt's EF hands, 1-based inclusive. */
const EF_HANDS = [
  [8, 43],
  [44, 79],
  [81, 116],
  [117, 149],
];

/** The share of equal residues over the window of `w` centered on each pair (i, j), a-major; null below a third. */
function windowIdentity(s: string, w = 9): (number | null)[] {
  const n = s.length;
  const h = w >> 1;
  return Array.from({ length: n * n }, (_, x) => {
    const [i, j] = [Math.floor(x / n), x % n];
    let same = 0;
    let seen = 0;
    for (let d = -h; d <= h; d++)
      if (i + d >= 0 && j + d >= 0 && i + d < n && j + d < n) {
        seen++;
        if (s[i + d] === s[j + d]) same++;
      }
    // Below a third, no value: the plot stays the surface there, in any theme.
    return seen > 0 && same / seen >= 1 / 3 ? same / seen : null;
  });
}

/** Words of `k` letters found at two places (i ≠ j), runs along a diagonal merged: pairs of ranges. */
function sharedWords(s: string, k = 4) {
  const cols = { start: [] as number[], stop: [] as number[], b_start: [] as number[], b_stop: [] as number[] };
  const n = s.length;
  const word = (i: number, d: number) =>
    i >= 0 && i + d >= 0 && i + d + k <= n && s.slice(i, i + k) === s.slice(i + d, i + d + k);
  for (let d = 1 - n; d < n; d++)
    for (let i = Math.max(0, -d); d !== 0 && i + k <= n; i++) {
      if (!word(i, d) || word(i - 1, d)) continue;
      let e = i;
      while (word(e + 1, d)) e++;
      cols.start.push(i);
      cols.stop.push(e + k);
      cols.b_start.push(i + d);
      cols.b_stop.push(e + k + d);
    }
  return cols;
}

export function calmodulinDotPlot(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Calmodulin against itself' });
  const s = b.data.sequence({ ref: 'cam', label: 'Calmodulin', value: CALMODULIN });
  const ef = b.data.annotation({
    ref: 'ef-hands',
    geometry: 'ranges',
    target: s,
    value: {
      start: EF_HANDS.map((r) => r[0] - 1),
      stop: EF_HANDS.map((r) => r[1]),
      label: ['EF 1', 'EF 2', 'EF 3', 'EF 4'],
    },
  });
  const identity = b.data.annotation({
    ref: 'identity',
    label: 'Identity over 9',
    geometry: 'grid',
    value: { target: ['cam'], b_target: ['cam'], value: [windowIdentity(CALMODULIN)] },
  });
  const words = b.data.annotation({
    ref: 'words',
    label: 'Shared words',
    geometry: 'pairs',
    target: s,
    value: sharedWords(CALMODULIN),
  });
  const sec = b.view({ ref: 'self', label: 'Against itself' }).section({ ref: 'self.cam' }).block({ sequence: s });
  sec.track({ label: '' }).ruler();
  sec.track({ label: 'Residues' }).letters();
  sec.track({ label: 'EF hands' }).ranges({ data: ef, labels: true });
  const plot = sec.track({ label: 'Calmodulin × calmodulin' });
  plot.dot_plot({
    data: identity,
    color: {
      theme: 'field_palette',
      params: {
        field: 'value',
        palette: {
          kind: 'continuous',
          mode: 'absolute',
          colors: [
            ['#a9bdd6', 1 / 3],
            ['#2b5c8a', 1],
          ],
        },
      },
    },
  });
  plot.dot_plot({ data: words, color: { theme: 'uniform', params: { color: '#d1495b' } } });
  return b.getDocument();
}
