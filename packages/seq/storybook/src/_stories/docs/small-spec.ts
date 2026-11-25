/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';

// A small spec: a protein, its domains and a signal as data; one view of one section with three tracks.

export function build(): Node {
  const b = createSeqSpec({ title: 'A small spec' });
  const letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP';
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: [0, 39], end: [23, 57], label: ['TAD1', 'TAD2'], category: ['tad', 'tad'] },
  });
  const signal = b.data.values('signal', {
    target: p53,
    values: Float32Array.from({ length: letters.length }, (_, i) => 50 + 40 * Math.sin(i / 5)),
  });
  b.views.view('main', (v) => {
    const section = v.section('protein', { of: p53 });
    section.track('letters').layer('letters', { type: 'letters', color: theme.scheme('zappo') });
    section.track('domains').layer('domains', { type: 'ranges', data: domains, color: theme.field('category') });
    section.track('signal').layer('signal', { type: 'values', data: signal, domain: [0, 100], encoding: 'line' });
  });
  return b.build();
}
