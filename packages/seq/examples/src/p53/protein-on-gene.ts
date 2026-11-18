/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, GeneExons, GeneFile, P53 } from './data.js';

// The protein on its gene: the gene's section follows the protein's through the codon map (minus strand, three exons),
// the domains cut into their exons. The gene's letters are read from the committed file, by location: relative to the
// page (Storybook) or to the package folder (Node), through the assets port.

export const title = 'p53 on its gene';

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const file = b.data.source('tp53-gene', { format: 'seq', uri: GeneFile });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', source: file });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GeneExons, strand: '-' });
  const domains = addDomains(b, p53);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    protein.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    const g = v.section('gene', { of: gene, follow: protein });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters' });
    g.track('domains', { label: 'Domains on the gene' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}
