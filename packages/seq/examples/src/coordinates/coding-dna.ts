/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { translateCodon } from '@molstar/seqx-core/resolve';

// A protein over its coding DNA (the old Seq*'s case 14): one gene section on the plus strand in genomic order, the
// whole fragment shown (the UTRs and the intron not collapsed), the protein's letters a codon wide through a codon map
// whose exons are the two coding pieces, and a domain through the map, cut at the intron and joined across it by a
// connector. Synthetic: a short made-up gene, its protein translated from it (the stop excluded).

export const title = 'A protein over its coding DNA';

/** The gene fragment: 5′ UTR, exon 1's coding part, the intron, exon 2's coding part and its stop, 3′ UTR. */
export const Gene =
  'GAATTCGCCACC' + // 5′ UTR (0–12)
  'ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC' + // exon 1's coding part (12–46)
  'GTAAGTCTTTCCTAACTTCTGTTTCCACAG' + // the intron (46–76)
  'ATTGGTACTTCAACAGCCAGGTGGAGTAA' + // exon 2's coding part and its stop (76–105)
  'GCGGCCGCATAAAGT'; // 3′ UTR (105–120)

/** The coding pieces on the gene: 0-based, half-open, in genomic order. */
export const Cds = { start: [12, 76], end: [46, 105] };

/** The exons, the UTRs with them: `Exon 1` 0–46, `Exon 2` 76–120. */
export const Exons = {
  start: new Int32Array([0, 76]),
  end: new Int32Array([46, 120]),
  id: ['E1', 'E2'],
  label: ['Exon 1', 'Exon 2'],
  category: ['exon-1', 'exon-2'],
};

/** The domain: residues 5 to 14, 1-based (0-based half-open 4–14). */
export const Domain = { start: new Int32Array([4]), end: new Int32Array([14]), id: ['D'], label: ['Domain'] };

/** The phosphosite: residue 17, 1-based. */
export const Site = { start: new Int32Array([16]), end: new Int32Array([17]), id: ['pS17'], label: ['pS17'] };

/** The coding sequence, the pieces joined (the stop included). */
export function coding(): string {
  let s = '';
  for (let i = 0; i < Cds.start.length; i++) s += Gene.slice(Cds.start[i], Cds.end[i]);
  return s;
}

/** The protein translated from the coding sequence, the stop excluded. */
export function protein(): string {
  const c = coding();
  let p = '';
  for (let i = 0; i + 3 <= c.length; i += 3) p += translateCodon(c.slice(i, i + 3));
  return p.replace(/\*$/, '');
}

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const gene = b.data.sequence('gene', { alphabet: 'dna', letters: Gene });
  const prot = b.data.sequence('protein', { alphabet: 'protein', letters: protein() });
  b.data.map('protein-on-gene', { from: prot, to: gene, method: 'codon', exons: Cds, strand: '+' });
  const exons = b.data.ranges('exons', { target: gene, rows: Exons });
  const domain = b.data.ranges('domain', { target: prot, rows: Domain });
  const site = b.data.ranges('site', { target: prot, rows: Site });
  b.views.view('main', (v) => {
    const g = v.section('gene', { of: gene, direction: 'genomic', pieces: 'whole' });
    g.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    g.track('dna', { label: 'DNA' }).layer('letters', { type: 'letters' });
    g.track('protein', { label: 'Protein' }).layer('letters', { type: 'letters', data: prot });
    g.track('exons', { label: 'Exons' }).layer('ranges', {
      type: 'ranges',
      data: exons,
      color: theme.field('category'),
    });
    g.track('domain', { label: 'Domain' }).layer('ranges', { type: 'ranges', data: domain });
    g.track('site', { label: 'Phosphosite' }).layer('variants', { type: 'variants', data: site });
  });
  return b.build();
}
