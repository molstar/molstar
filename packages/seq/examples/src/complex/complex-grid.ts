/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { Data } from './data.js';
import { PaeColor, predictionData } from './pae.js';

// The grid over the complex (the whole-complex PAE plot): the prediction's chains side by side as blocks of one tracks
// section (CDK2, then cyclin A2), and a grid section that follows it on both axes (`follow`, `follow_y: x`): its x and y
// are the blocks' display axis, the chains in the section's order, so the plot is the N×N matrix over the complex
// (298 + 432 residues, synthetic), its chain boundaries the blocks' gaps across it. The matrix is read through the matrix
// read and held as one grid node a pair of chains (a grid draws only on its own systems), all four in one track.

export const title = 'Complex: PAE over the complex (synthetic)';

export function build(d: Derived = Data): Node {
  const b = createSeqSpec({ title });
  const p = predictionData(b, d);
  const chains = ['cdk2', 'cyclin'] as const;
  b.views.view('main', (v) => {
    const s = v.section('chains', {});
    s.block('cdk2', { of: p.cdk2.protein, label: 'CDK2' });
    s.block('cyclin', { of: p.cyclin.protein, label: 'cyclin A2' });
    s.track('confidence', { label: 'pLDDT' }, (t) => {
      t.layer('cdk2', { type: 'values', data: p.cdk2.plddt });
      t.layer('cyclin', { type: 'values', data: p.cyclin.plddt });
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
    v.section('pae', { layout: 'grid', follow: s.address }).track('pae', { label: 'PAE (Å, synthetic)' }, (t) => {
      for (const x of chains)
        for (const y of chains)
          t.layer(`${x}-${y}`, { type: 'grid', data: p.grid(x, y), aggregate: 'min', color: PaeColor() });
    });
  });
  return b.build();
}
