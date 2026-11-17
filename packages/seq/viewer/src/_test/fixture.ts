/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A protein on its gene, small enough to state every answer by hand: 4 residues coded by two exons, [2, 7) and
// [12, 19) of a 24-nucleotide gene, so the second codon is split by the intron. On the plus strand residue 1 is
// nucleotides 5, 6 and 12; on the minus strand residue 2 is 12, 6 and 5.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

export const Protein = '/data/p';
export const Gene = '/data/g';
export const View = '/views/main';
export const ProteinSection = '/views/main/protein';
export const TwinSection = '/views/main/twin';
export const GeneSection = '/views/main/gene';

export interface FixtureOptions {
  readonly strand?: '+' | '-';
  /** The protein sections' range. */
  readonly range?: readonly [number, number];
  /** A second protein section on the same axis (a scroll group). */
  readonly twin?: boolean;
  readonly title?: string;
  /** The gene's length (24). */
  readonly geneLength?: number;
}

export function proteinOnGene(options: FixtureOptions = {}): Node {
  const b = createSeqSpec({ title: options.title ?? 'p on g' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(6).slice(0, options.geneLength ?? 24) });
  b.data.map('m', {
    from: p,
    to: g,
    method: 'codon',
    exons: { start: [2, 12], end: [7, 19] },
    strand: options.strand ?? '+',
  });
  b.views.view('main', (v) => {
    const range = options.range ? { range: [...options.range] } : {};
    const protein = v.section('protein', { of: p, ...range });
    protein.track('letters').layer('l', { type: 'letters' });
    if (options.twin) v.section('twin', { of: p, ...range });
    v.section('gene', { of: g, follow: protein }).track('letters').layer('l', { type: 'letters' });
  });
  return b.build();
}

/** A Location's parts as plain arrays, to compare. */
export function partsOf(loc: { parts: readonly { system: string; ranges: Int32Array }[] } | undefined) {
  return loc?.parts.map((p) => [p.system, [...p.ranges]]);
}
