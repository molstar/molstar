/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// A Manhattan-like plot (EC-216, EC-217; the values layer's reference lines): SNPs as dots at their height, -log10 p,
// over one chromosome of a toy association study, a dashed line at the genome-wide threshold (7.3, p = 5e-8) and a
// dotted one at the suggestive one (5). The data stays under the threshold (a peak to about 6): over an open domain
// the lines are part of the scale, so the plot shows that no SNP reached it. Inline and generated from a fixed seed
// (a small LCG), never a file: the whole is under 100 KB.

/** The chromosome's length, positions, and a SNP's spacing. */
export const ManhattanLength = 20_000;
export const ManhattanSpacing = 20;
/** The association's peak: its centre and its width, positions, and the height its SNPs rise to (-log10 p). */
export const ManhattanPeak = { at: 12_000, width: 160, height: 6 } as const;
/** The thresholds, -log10 p: genome-wide (5e-8) and suggestive. */
export const GenomeWide = 7.3;
export const Suggestive = 5;

/** The SNPs' heights, a value a SNP and NaN between (missing, never zero): the null's tail plus the peak's rise. */
export function manhattanValues(seed = 7): Float32Array {
  let s = seed >>> 0;
  const next = () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return (s + 1) / 0x100000001;
  };
  const out = new Float32Array(ManhattanLength).fill(Number.NaN);
  for (let p = ManhattanSpacing / 2; p < ManhattanLength; p += ManhattanSpacing) {
    // under the null, p is uniform: -log10 p is exponential with a mean of 0.43
    let v = -Math.log10(next());
    const d = (p - ManhattanPeak.at) / ManhattanPeak.width;
    v = Math.max(v, ManhattanPeak.height * Math.exp(-d * d) * (0.75 + 0.25 * next()));
    out[p] = Math.min(v, ManhattanPeak.height);
  }
  return out;
}

export function manhattan(): Node {
  const b = createSeqSpec({ title: 'A Manhattan plot' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'A'.repeat(ManhattanLength) });
  const gwas = b.data.values('gwas', { target: chr, values: manhattanValues() });
  b.views.view('main', (v) => {
    const s = v.section('chr', { of: chr });
    s.track('ruler', { label: 'chr' }).layer('ruler', { type: 'ruler', data: chr });
    s.track('gwas', { label: '-log10 p', height: 120 }).layer('values', {
      type: 'values',
      data: gwas,
      encoding: 'points',
      reference: [
        { at: GenomeWide, label: 'genome-wide', style: 'dashed' },
        { at: Suggestive, label: 'suggestive', style: 'dotted' },
      ],
    });
  });
  return b.build();
}
