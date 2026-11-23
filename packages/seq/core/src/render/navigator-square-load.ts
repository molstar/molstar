/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedView } from '@molstar/seqx-core/view';
import { isSquareNavigator } from './navigator-square-kind.js';

// The square's door (_spec/README.md "Render", Navigators, The square: lazy): `navigator-square.ts` and what it imports
// (its cells, its set, its lattice) load the first time a view shows a square (an export waits for it, a viewer's
// navigators part too), so a composition whose views have no square pays nothing for one. Why lazy: the square is
// one navigator form of several, about 16 kB.

type SquareModule = typeof import('./navigator-square.js');

export { isSquareNavigator };

/** The square's chunk once landed (`loadSquareNavigator`). */
export const SquareChunk: { module?: SquareModule } = {};

let loading: Promise<SquareModule> | undefined;

/** Loads the square's chunk (once; a failed load is tried again by the next ask). */
export function loadSquareNavigator(): Promise<SquareModule> {
  loading ??= import('./navigator-square.js').then(
    (m) => (SquareChunk.module = m),
    (e) => {
      loading = undefined;
      throw e;
    },
  );
  return loading;
}

/** Whether a view shows a square the chunk is not landed for (a draw of it waits: `loadSquareNavigator`). */
export function squareDue(view: ResolvedView): boolean {
  if (SquareChunk.module) return false;
  for (const n of view.navigators) if (!n.hidden && isSquareNavigator(n, view)) return true;
  return false;
}
