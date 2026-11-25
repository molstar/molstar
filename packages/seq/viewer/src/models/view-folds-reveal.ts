/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a reveal reads of the folded groups (_spec/README.md "ViewModel", Groups): loaded through `groups-door.ts` the
// first time a spec holds a section with a group, so a viewer that never shows one holds none of it.

import type { ResolvedSection } from '@molstar/seqx-core/view';
import { groupIdOf } from './view-folds.js';

/**
 * The folded groups a reveal opens (_spec/README.md "ViewModel", Groups): the groups of the view's sections whose shown
 * tracks alone draw what the data owners make (an item's key: a layer's data, or a form over it); none when a track
 * not folded draws it (the reveal goes there) or nothing does.
 */
export function foldedHolding(
  sections: readonly ResolvedSection[],
  folded: (section: string) => readonly string[],
  owners: ReadonlySet<string>,
): string[] {
  const out: string[] = [];
  for (const s of sections) {
    const ids = folded(s.address);
    for (const t of s.tracks) {
      if (!t.layers.some((l) => l.data && owners.has(ownerOf(l.data.address)))) continue;
      const g = t.group;
      // an authored summary track stands for its folded group: drawn while folded
      if (!g || t.summaryOf || !ids.includes(groupIdOf(g))) return [];
      if (!out.includes(g)) out.push(g);
    }
  }
  return out;
}

/** A data address's node part (`/data/x/(item …)` → `/data/x`). */
const ownerOf = (address: string): string => {
  const k = address.indexOf('/(');
  return k < 0 ? address : address.slice(0, k);
};
