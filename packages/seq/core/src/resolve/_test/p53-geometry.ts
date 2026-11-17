/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The geometry of human TP53 (P04637, 393 residues) on chr17, GRCh38, minus strand: its ten coding exons (exons 2 to
// 11 of NM_000546, the stop codon included), 0-based and half-open, in genomic order. Positions only: no letters are
// needed to place a protein on its gene. The lengths add up to 1,182 nucleotides, 3 × 393 + a stop.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** The coding exons, genomic order (exon 11 first, exon 2 last). */
export const TP53_EXONS = {
  start: [7669608, 7670608, 7673534, 7673700, 7674180, 7674858, 7675052, 7675993, 7676381, 7676520],
  end: [7669690, 7670715, 7673608, 7673837, 7674290, 7674971, 7675236, 7676272, 7676403, 7676594],
};

export const TP53_LENGTH = 393;

/** p53 (placeholder letters: the length is what counts) on chr17 read from a source, a domains annotation, a view. */
export function p53GeometrySpec(): Node {
  const b = createSeqSpec({ title: 'p53 on chr17' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: `M${'X'.repeat(TP53_LENGTH - 1)}` });
  const chr17 = b.data.source('grch38', { format: 'fasta', uri: 'GRCh38.fa', index: 'GRCh38.fa.fai' });
  const gene = b.data.sequence('chr17', { alphabet: 'dna', source: chr17, name: '17' });
  b.data.map('p53-on-chr17', { from: p53, to: gene, method: 'codon', exons: TP53_EXONS, strand: '-' });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 93, 323]), end: new Int32Array([61, 292, 355]), id: ['TAD', 'DBD', 'TET'] },
  });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('domains').layer({ type: 'ranges', data: domains });
    const g = v.section('gene', { of: gene, follow: protein });
    g.track('domains').layer({ type: 'ranges', data: domains });
  });
  return b.build();
}
