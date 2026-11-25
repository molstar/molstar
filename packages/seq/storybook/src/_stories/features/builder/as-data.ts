/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// `build()` is the spec as data: a normalized tree, JSON by default, and what `printSexpr` prints. A host that never
// writes with the builder sends the same.

export function build(): Node {
  const b = createSeqSpec({ title: 'build() as data' });
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' });
  b.views.view('main', (v) => {
    v.section('seq', { of: seq }).track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}
