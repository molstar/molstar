/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { printLink } from '@molstar/seqx-core/query';
import { compose, type ComposedOptions, Entry, structureSection, unitAddress } from './composed.js';

// Addressing, with a structure beside the view (examples.md, Showcases; addressing.md): a short text whose links act on
// addresses in Mol*'s link grammar (`!seq-focus=/refs/dbd`) and reach the view and, through the structure link, the fake
// structure viewer: a focus moves its camera frame, a highlight and a selection light residues of every chain SIFTS pairs
// them with. One residue is named in three numberings (UniProt's, the entity's `label_seq_id` and the chain's author
// numbers); "Copy link" turns what the view selects, or the pointer is over, back into a link. A link is data: its text
// and its href, never markup. The refs are bookmarks the text names.

export const title = 'Addressing with a structure';

/** A run of a passage: plain text, or a link's text and its href (`!seq-…`). */
export type Run = string | { readonly text: string; readonly href: string };

const link = (text: string, action: 'focus' | 'highlight' | 'select', address: string): Run => ({
  text,
  href: printLink({ actions: [{ action, address }] }),
});

/** The text beside the view: each link an action on an address. */
export const Passages: readonly (readonly Run[])[] = [
  [
    'p53 binds DNA through its ',
    link('DNA-binding domain', 'focus', '/refs/dbd'),
    ', which the crystal 2OCJ holds in four copies: ',
    link('hover or select', 'select', '/refs/dbd'),
    ' it and the same residues light on all four chains.',
  ],
  [
    'Most cancer mutations fall on a few ',
    link('hotspot residues', 'highlight', '/refs/hotspots'),
    ' (',
    link('select them', 'select', '/refs/hotspots'),
    ') and ',
    link('residues 245 to 249', 'select', '/data/protein/(range :from 245 :to 249)'),
    '.',
  ],
  [
    'One residue in three numberings: ',
    link('R175 of UniProt', 'focus', '/data/protein/(at :label 175)'),
    ', ',
    link('residue 82 of the entity', 'focus', `${unitAddress('A')}/(at :label 82 :numbering ordinal)`),
    ' and ',
    link('175 in chain A’s author numbers', 'focus', `${unitAddress('A')}/(at :label 175 :numbering auth)`),
    ': the same place.',
  ],
];

export function build(options: ComposedOptions = {}): Node {
  const b = createSeqSpec({ title });
  const c = compose(b, options);
  b.refs.ref('dbd', {
    to: '/data/protein/(range :from 94 :to 292)',
    label: 'DNA-binding domain',
    color: '#d62728',
  });
  b.refs.ref('hotspots', {
    to: ['/data/protein/(at :label 175)', '/data/protein/(at :label 248)', '/data/protein/(at :label 273)'],
    label: 'Hotspots R175, R248, R273',
  });
  b.views.view('main', (v) => {
    v.section('protein', { of: c.protein }, (s) => {
      s.track('ruler', { label: 'P04637' }).layer('ruler', { type: 'ruler' });
      s.track('letters', { label: 'Sequence (UniProt)' }).layer('letters', { type: 'letters' });
      s.track('dna', { label: 'DNA binding (UniProt)' }).layer('dna', {
        type: 'ranges',
        data: c.features('DNA binding'),
      });
      if (c.coverage)
        s.track('coverage', { label: `${Entry} chains (SIFTS)` }).layer('coverage', {
          type: 'ranges',
          data: c.coverage,
        });
    });
    structureSection(v, c, '/views/main/protein');
  });
  return b.build();
}
