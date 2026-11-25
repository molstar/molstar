/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, measure, theme } from '@molstar/seqx-core/lang';

// A layer's channels take themes: a scheme, a field of its data (categorical), rules over a measure (a value).

export function build(): Node {
  const b = createSeqSpec({ title: 'Layers and themes' });
  const letters = 'MKTAYIAKQRQISFVKSHFSRQMKTAYIAKQRQISFVKSHFSRQ';
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters });
  const regions = b.data.ranges('regions', {
    target: seq,
    rows: { start: [0, 10, 24], end: [8, 22, 40], category: ['helix', 'strand', 'helix'] },
  });
  const signal = b.data.values('signal', {
    target: seq,
    values: Float32Array.from({ length: letters.length }, (_, i) => Math.round(50 + 45 * Math.sin(i / 4))),
  });
  b.views.view('main', (v) => {
    const section = v.section('seq', { of: seq });
    section.track('letters').layer('letters', { type: 'letters', color: theme.scheme('zappo') });
    // categorical: a color for each value of the column `category`
    section.track('regions').layer('regions', { type: 'ranges', data: regions, color: theme.field('category') });
    // values: the bars in the alert color where the value is over 80, in the base color elsewhere
    section.track('signal').layer('signal', {
      type: 'values',
      data: signal,
      domain: [0, 100],
      color: theme.rules([theme.rule([theme.gt(measure.value(), 80)], '#d62728', 'high')], '#4682b4'),
    });
  });
  return b.build();
}
