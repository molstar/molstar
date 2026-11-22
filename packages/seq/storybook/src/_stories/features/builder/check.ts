/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// A mistake the types cannot catch: an address written out as a string that names nothing. `check()` lists the issues
// of the spec so far.

function write() {
  const b = createSeqSpec({ title: 'check()' });
  const seq = b.data.sequence('hello', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' });
  b.views.view('main', (v) => {
    const section = v.section('seq', { of: seq });
    section.track('letters').layer('letters', { type: 'letters' });
    // any string is a reference, so this compiles; /data/helo names nothing
    section.track('hits').layer('hits', { type: 'ranges', data: '/data/helo' });
  });
  return b;
}

export const build = (): Node => write().build();
export const issues = () => write().check();
