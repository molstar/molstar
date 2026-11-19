/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// The Complex dataset, a stand-in until its dataset is recorded: CDK2 with cyclin A as in 1FIN (two copies of the
// complex, an ATP on each CDK2), its chains as blocks of one section, polymers first, then the ligands as
// one-position blocks (EC-101); the copies merged into one block labeled ×2 (`copies: 'merge'`, the default, EC-65) or
// split (`each`). The chains' lengths and features' places are CDK2's and cyclin A's construct's (PSTAIRE 45–56, the
// activation loop with T160, the cyclin box); their letters are generated from a seed (never stored): the recorded
// entry replaces them.

export const title = 'Complex blocks';

/** The chains' lengths (CDK2, the cyclin A construct of 1FIN). */
export const ComplexLengths = { cdk2: 298, cyclin: 260 } as const;

/** A chain's stand-in letters. */
export function letters(seed: number, n: number): string {
  const r = new PCG(seed);
  const aa = 'ACDEFGHIKLMNPQRSTVWY';
  let s = '';
  for (let i = 0; i < n; i++) s += aa[r.int() % aa.length];
  return s;
}

export interface BlocksOptions {
  readonly arrange?: 'row' | 'stack';
  readonly copies?: 'merge' | 'each';
  readonly layout?: 'rings';
}

/** The complex's chains and ligands as blocks, its features under them. */
export function build(options: BlocksOptions = {}): Node {
  const b = createSeqSpec({ title });
  const cdk2 = b.data.sequence('cdk2', { alphabet: 'protein', letters: letters(2, ComplexLengths.cdk2) });
  const cyclin = b.data.sequence('cyclin', { alphabet: 'protein', letters: letters(3, ComplexLengths.cyclin) });
  // ATP as one position (its per-atom system waits for the dataset)
  const atp = b.data.sequence('atp', { alphabet: 'protein', letters: 'X', name: 'ATP' });
  const cdk2f = b.data.ranges('cdk2-features', {
    target: cdk2,
    rows: {
      start: [44, 144, 159],
      end: [56, 172, 160],
      id: ['psta', 'aloop', 't160'],
      label: ['PSTAIRE', 'activation loop', 'T160'],
    },
  });
  const cyclinf = b.data.ranges('cyclin-features', {
    target: cyclin,
    rows: { start: [40, 140], end: [130, 230], id: ['box1', 'box2'], label: ['cyclin box N', 'cyclin box C'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('complex', {
      ...(options.arrange ? { arrange: options.arrange } : {}),
      ...(options.copies ? { copies: options.copies } : {}),
      ...(options.layout ? { layout: options.layout } : {}),
    });
    s.block('a', { of: cdk2, label: 'CDK2 (A)' });
    s.block('b', { of: cyclin, label: 'cyclin A (B)' });
    s.block('c', { of: cdk2, label: 'CDK2 (C)' });
    s.block('d', { of: cyclin, label: 'cyclin A (D)' });
    s.block('atp-a', { of: atp, label: 'ATP' });
    s.block('atp-c', { of: atp, label: 'ATP' });
    s.track('features', { label: 'Features' }, (t) => {
      t.layer('cdk2', { type: 'ranges', data: cdk2f });
      t.layer('cyclin', { type: 'ranges', data: cyclinf });
    });
    s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The copies split: a block a chain (`copies: 'each'`). */
export const split = () => build({ copies: 'each' });

/** The chains stacked: a block a row, aligned at their first residue (`arrange: 'stack'`). */
export const stacked = () => build({ arrange: 'stack' });

/** The chains on a ring: a block an arc. */
export const ring = () => build({ layout: 'rings' });
