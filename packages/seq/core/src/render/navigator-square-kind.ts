/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedNavigator, ResolvedView } from '@molstar/seqx-core/view';

// Which navigators are squares (_spec/README.md "Render", Navigators, The square): asked without the square's chunk
// (`navigator-square-load.ts`), so a module that asks never reaches it (a dynamic import is a chunk of every bundle
// that reaches its module).

/** Whether a navigator is a square: two axes over blocks, bound to two sections of blocks (or one section twice). */
export function isSquareNavigator(nav: ResolvedNavigator, view: ResolvedView): boolean {
  if (!nav.blocks || nav.axes !== 'both' || nav.of.length !== 2) return false;
  return nav.of.every((a) => view.sections.find((s) => s.address === a)?.blocks?.length);
}
