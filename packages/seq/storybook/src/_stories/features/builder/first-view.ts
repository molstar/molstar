/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// The smallest spec that draws: letters, one track, one layer.

export function build(): Node {
  const b = createSeqSpec({ title: 'A first view' });
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' });
  b.views.view('main', (v) => {
    const section = v.section('seq', { of: seq });
    section.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}
