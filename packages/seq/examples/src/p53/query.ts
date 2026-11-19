/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, conservation, GeneExons, geneLetters, P53 } from './data.js';

// p53 for the query stories (`Features/Query`): the protein with its domains and conservation, and its gene reading
// along it through the codon map. All data is inline, so a query resolves at once. The addresses the stories ask are
// `/data/p53`, `/data/domains`, `/data/conservation` and `/data/tp53`.

export const title = 'p53, queried';

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GeneExons, strand: '-' });
  const domains = addDomains(b, p53);
  const values = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    protein.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    protein.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: values });
    const g = v.section('gene', { of: gene, follow: protein, direction: 'follow', pieces: 'coding' });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}
