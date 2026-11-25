/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The categorical axes' door (_spec/README.md "Render", The grid layer: categorical axes): the members' orders
// (`members.ts`) and the labels' chrome (`layout/grid-members.ts`) load the first time a spec holds a section whose
// `of` or `of_y` is a rows system, a model waiting for them as for the bridge chunk (`withViewKind`), so such a grid is
// never drawn without them. Why lazy: a categorical axis is rare; a view without one pays nothing for it.

import type { Node } from '@molstar/core/repr/tree';
import { memberSystem, withViewKind } from '@molstar/seqx-core/resolve';

export const MemberHooks: {
  /** The members of a categorical axis in its order (`members.ts`), once landed. */
  membersOf?: typeof import('./members.js').membersOf;
  /** The labels' chrome (`layout/grid-members.ts`), once landed. */
  chrome?: typeof import('./layout/grid-members.js');
} = {};

/** Loads the categorical axes' chunk. */
export function loadMembers(): Promise<void> {
  return Promise.all([import('./members.js'), import('./layout/grid-members.js')]).then(([m, c]) => {
    MemberHooks.membersOf = m.membersOf;
    MemberHooks.chrome = c;
  });
}

/** Whether a section names a rows system on an axis (`of`, `of_y`). */
function categorical(n: Node): boolean {
  const of = n.params?.of;
  const ofY = n.params?.of_y;
  return (
    (typeof of === 'string' && memberSystem(of) !== undefined) ||
    (typeof ofY === 'string' && memberSystem(ofY) !== undefined)
  );
}

let registered = false;

/** Registers the chunk with resolve (the built-ins call it): loaded when a spec holds a categorical section. */
export function registerMembers(): void {
  if (registered) return;
  registered = true;
  withViewKind('section', loadMembers, categorical);
}
