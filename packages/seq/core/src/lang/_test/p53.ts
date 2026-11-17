/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '../builder.js';

/** The first residues of human TP53 (UniProt P04637). */
export const P53_LETTERS = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN';

/** Its gene's exons on the minus strand of a toy region: 50 + 40 positions, 30 codons. */
export const EXONS = { start: [100, 160], end: [150, 200] };

/** The p53 example of the language description, built through the builder. */
export function p53Spec(): Node {
  const b = createSeqSpec({ title: 'p53' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const chr17 = b.data.source('chr17', {
    format: 'fasta',
    uri: 'data/chr17-region.fa',
    index: 'data/chr17-region.fa.fai',
  });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', source: chr17 });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: EXONS, strand: '-' });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 17]), end: new Int32Array([14, 30]), label: ['Disordered', 'TAD'] },
  });
  b.data.values('conservation', { target: p53, values: new Float32Array(P53_LETTERS.length).fill(0.5) });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('letters').layer({ type: 'letters' });
    protein.track('domains', { label: 'Domains' }).layer({ type: 'ranges', data: domains });
    v.section('gene', { of: gene, follow: protein });
  });
  return b.build();
}
