/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, conservation, P53 } from './data.js';

// p53 as tracks: its letters, its domains (in lanes), its conservation (two residues missing). All data is inline.

export const title = 'p53 tracks';

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const values = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer({ type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer({ type: 'ranges', data: domains });
    s.track('conservation', { label: 'Conservation' }).layer({ type: 'values', data: values });
  });
  return b.build();
}
