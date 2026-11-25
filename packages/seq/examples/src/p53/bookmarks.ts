/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { addRefPatch, type RefInput, SeqLanguage } from '@molstar/seqx-core/lang';
import * as geneByRegion from './gene-by-region.js';
import { GeneOffset, RegionSequence, RegionSource } from './region.js';
import { GeneLength } from './data.js';

// Bookmarks on the gene's region: the gene-by-region view with refs added, each by one patch. Two name transcripts
// that arrive as dynamic nodes only when their region is read: such a ref waits (an issue, not an error) until its node
// is held, then resolves; a focus asked of it meanwhile runs when it arrives.

export const title = 'p53’s gene region, bookmarked';

/** The bookmarks the view starts with. */
export const Bookmarks: readonly RefInput[] = [
  {
    id: 'gene',
    to: `${RegionSequence}/(range :from ${GeneOffset + 1} :to ${GeneOffset + GeneLength})`,
    label: 'TP53’s gene',
  },
  { id: 't2', to: `${RegionSource}/(node :key T2)`, label: 'Transcript T2' },
  { id: 'wrap53', to: `${RegionSource}/(node :key W1)`, label: 'WRAP53 (a neighbor)' },
];

/** The view's spec: gene by region, and its bookmarks, each added as one patch. */
export function build(): Node {
  let root = geneByRegion.build();
  for (const ref of Bookmarks) {
    const r = applyPatch(root, addRefPatch(root, ref), { language: SeqLanguage });
    if (!r.ok) throw new Error(`The bookmark ${ref.id} did not apply: ${r.issues.map((i) => i.code).join(', ')}`);
    root = r.root;
  }
  return root;
}

export const registry = geneByRegion.registry;
