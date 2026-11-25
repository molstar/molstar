/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The track groups' door (_spec/README.md "Render", Track groups: lazy): the folding, the summary bands and the
// headers' layout (`groups.ts`) and their drawing and hit (`groups-draw.ts`) load the first time a spec holds a
// section with a group, a model waiting for them as for the categorical axes (`withViewKind`), so a section with
// groups is never drawn without them. Why lazy: groups are rare; a view without one pays nothing for them.

import { withViewKind } from '@molstar/seqx-core/resolve';
import { holdsGroup } from '@molstar/seqx-core/view';

/** The chevron's room before a group's label, px (here: the label column's measure reads it without the chunk). */
export const ChevronPx = 14;
/** The gap between a group's label and its count, px. */
export const CountGapPx = 6;

/** A folded group's quiet count: its tracks ("1 track", "3 tracks") and what its summary leaves out. */
export function groupCount(tracks: number, left = 0): string {
  const n = `${tracks} track${tracks === 1 ? '' : 's'}`;
  return left ? `${n} · +${left} not shown` : n;
}

export const GroupsHooks: {
  /** The folding of a section's groups (`groups.ts`), once landed. */
  fold?: typeof import('./groups.js').foldGroups;
  /** The headers' layout (`groups.ts`), once landed. */
  headers?: typeof import('./groups.js').groupHeaders;
  /** The headers' drawing (`groups-draw.ts`), once landed. */
  draw?: typeof import('./groups-draw.js').drawGroupHeaders;
  /** The headers' hit (`groups-draw.ts`), once landed. */
  hit?: typeof import('./groups-draw.js').hitGroupHeader;
} = {};

/** Loads the track groups' chunk. */
export function loadGroups(): Promise<void> {
  return Promise.all([import('./groups.js'), import('./groups-draw.js')]).then(([g, d]) => {
    GroupsHooks.fold = g.foldGroups;
    GroupsHooks.headers = g.groupHeaders;
    GroupsHooks.draw = d.drawGroupHeaders;
    GroupsHooks.hit = d.hitGroupHeader;
  });
}

let registered = false;

/** Registers the chunk with resolve (the tracks layout calls it): loaded when a spec holds a section with a group. */
export function registerGroups(): void {
  if (registered) return;
  registered = true;
  withViewKind('section', loadGroups, holdsGroup);
}
