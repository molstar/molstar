/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A test fixture: a heavy chain's start with Kabat and IMGT numberings (declared `numbering` children).

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

/** A heavy chain's start, numbered twice: Kabat (an expression tag unnumbered, insertion codes) and IMGT. */
export const H = 'EVQLVESGGGLVQPGGSLRL';
// biome-ignore format: a label a position, in index order
export const Kabat = [null, '1', '2', '3', '4', '5', '52', '52A', '53', '54', '82', '82A', '82B', '82C', '83', '100', '100A', '100B', '101', '102'];
// IMGT's CDR3 is numbered from both ends inward: 111, 111.1 then 112.1, 112, out of the order of labels (EC-7)
// biome-ignore format: a label a position, in index order
export const Imgt = ['1', '2', '3', '4', '5', '6', '55', '56', '57', '58', '105', '106', '107', '108', '109', '111', '111.1', '112.1', '112', '113'];

export function antibodySpec(): Node {
  const b = createSeqSpec({ title: 'numberings' });
  b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) => {
    s.numbering('kabat', { method: 'labels', labels: Kabat, scheme: 'kabat' });
    s.numbering('imgt', { method: 'labels', labels: Imgt, scheme: 'imgt' });
  });
  return b.build();
}
