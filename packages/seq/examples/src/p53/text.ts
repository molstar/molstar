/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, addSecondary, P53 } from './data.js';

// The p53 text figure, three ways: the sequence as text (groups of ten, a number a line and a group, its secondary
// structure and domains as thin bands under the letters, a band left out of a line it is empty on); the same figure
// centered and capped, its group numbers inline and a number at each end of a line; and the text under a tracks
// overview in one view (a section's own `layout`). All data is inline.

export const title = 'p53 text';

/** The text figure: groups of ten, numbers at the start, group labels above; secondary structure and domains under. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const ss = addSecondary(b, p53);
  const domains = addDomains(b, p53);
  b.views.view('main', { layout: 'text' }, (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('secondary', { label: 'Secondary structure' }).layer('secondary', { type: 'secondary', data: ss });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}

/** The figure centered in at most 640 px, its group numbers inline, a number at each end of a line. */
export function centered(): Node {
  const b = createSeqSpec({ title: `${title}, centered` });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  b.views.view('main', { layout: 'text' }, (v) => {
    const s = v.section('protein', {
      of: p53,
      align: 'center',
      max_width: 640,
      numbers: 'both',
      group_labels: 'inline',
    });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}

/** A tracks overview of the protein and its text under it, in one view: the text section's own `layout`. */
export function underTracks(): Node {
  const b = createSeqSpec({ title: `${title} under its tracks` });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  b.views.view('main', (v) => {
    const o = v.section('overview', { of: p53 });
    o.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    const s = v.section('text', { of: p53, layout: 'text' });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}
