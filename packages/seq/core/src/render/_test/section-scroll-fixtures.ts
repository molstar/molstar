/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A tracks section past its cap: the fixture its tests and its allocation case share.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** A ruler and `n` gene tracks of one feature each; `height` the section's own. */
export function geneTracks(n = 30, height?: number | 'content'): Node {
  const b = createSeqSpec({ title: 'tracks past the cap' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, ...(height !== undefined ? { height } : {}) });
    sec.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    for (let i = 0; i < n; i++) {
      const d = b.data.ranges(`g${i}`, { target: s, rows: { start: [1 + i], end: [200 + i], id: [`gene ${i}`] } });
      sec.track(`t${i}`, { label: `Gene ${i}` }).layer('l', { type: 'ranges', data: d });
    }
  });
  return b.build();
}
