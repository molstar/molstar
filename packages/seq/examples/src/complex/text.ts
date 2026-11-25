/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { crystalData } from './blocks.js';
import type { Derived } from './derived.js';
import { Data } from './data.js';

// 1FIN as text (EC-99): its first copy's chains as blocks of one text section, CDK2 then cyclin A, then the ATP as a
// chip that flows with the letters like a word; the features and the interface as thin bands under the letters.

export const title = 'Complex: 1FIN as text';

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const x = crystalData(b, d);
  b.views.view('main', { layout: 'text' }, (v) => {
    const s = v.section('complex', { copies: 'each' });
    s.block('a', { of: x.cdk2, label: 'CDK2 (A)' });
    s.block('b', { of: x.cyclin, label: 'cyclin A (B)' });
    s.block('atp-a', { of: x.atp, label: 'ATP' });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
    s.track('features', { label: 'Features' }, (t) => {
      t.layer('cdk2', { type: 'ranges', data: x.cdk2Features });
      t.layer('cyclin', { type: 'ranges', data: x.cyclinFeatures });
    });
    s.track('interface', { label: 'Interface' }, (t) => {
      t.layer('cdk2', { type: 'ranges', data: x.cdk2Interface });
      t.layer('cyclin', { type: 'ranges', data: x.cyclinInterface });
    });
  });
  return b.build();
}
