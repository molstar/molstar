/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { printLink } from '@molstar/seqx-core/query';
import { addDomains, conservation, P53 } from './data.js';

// p53 addressed: the tracks view with refs naming its domains, and a short text whose links focus, highlight and
// select in the view. One location is linked in two numberings (ordinal and index), so the same residues are named
// twice. A link is data: its text and its href, never markup.

export const title = 'p53, addressed';

/** The view's spec: the tracks view and the refs the text's links name. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const values = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    s.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: values });
  });
  b.refs.ref('dbd', { to: `${domains}/(item :id DBD)`, label: 'DNA-binding domain', color: '#d62728' });
  b.refs.ref('tad', { to: `${domains}/(item :id TAD)`, label: 'Transactivation domain' });
  b.refs.ref('hotspots', {
    to: [`${p53}/(at :label 175)`, `${p53}/(at :label 248)`, `${p53}/(at :label 273)`],
    label: 'Hotspots R175, R248, R273',
  });
  return b.build();
}

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
    '; its N-terminus holds the ',
    link('transactivation domain', 'highlight', '/refs/tad'),
    '.',
  ],
  [
    'Most cancer mutations fall on a few ',
    link('hotspot residues', 'select', '/refs/hotspots'),
    ', such as ',
    link('residues 245 to 249', 'select', '/data/p53/(range :from 245 :to 249)'),
    '.',
  ],
  [
    'One residue in two numberings: ',
    link('R175 (1-based)', 'focus', '/data/p53/(at :label 175)'),
    ' is ',
    link('174 (0-based)', 'focus', '/data/p53/(at :label 174 :numbering index)'),
    '.',
  ],
];
