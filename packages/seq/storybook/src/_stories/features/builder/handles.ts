/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, type SectionBuilder, type TrackBuilder } from '@molstar/seqx-core/lang';

// Every call returns the address of what it adds, typed by its kind: pass it where a reference is wanted.

export function build(): Node {
  const b = createSeqSpec({ title: 'Handles' });
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' }, (s) => {
    // a numbering is written in the sequence's body; the call returns its address
    s.numbering('author', { method: 'labels', labels: Array.from({ length: 22 }, (_, i) => `${101 + i}`) });
  });
  const hits = b.data.ranges('hits', {
    target: seq, // the sequence's handle, where a `ranges` node names its target
    rows: { start: [2, 12], end: [8, 18], id: ['a', 'b'] },
  });
  b.views.view('main', (v) => {
    const section = v.section('seq', { of: seq });
    section.track('ruler').layer('ruler', { type: 'ruler', numbering: 'author' });
    section.track('letters').layer('letters', { type: 'letters' });
    section.track('hits').layer('hits', { type: 'ranges', data: hits });
  });
  return b.build();
}

/** Never called: a layer's `data` takes a sequence, ranges, values…, and a section is none of them. */
export function wrongKind(section: SectionBuilder, track: TrackBuilder) {
  // @ts-expect-error Address<'section'> is not a RefValue of a data kind
  track.layer('wrong', { type: 'ranges', data: section.address });
}
