/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultColumnSize, defaultLayout, normalize, RightGroup, type ShellLayout, type ViewRef } from '../layout.js';

/** The layout a model makes as views register at the right group: the main area, the right group beside it. */
export function rightLayout(views: readonly ViewRef[], active?: ViewRef): ShellLayout {
  return normalize({
    ...defaultLayout(),
    root: {
      kind: 'split',
      id: 'root',
      dir: 'row',
      children: [
        { kind: 'main', id: 'main' },
        { kind: 'group', id: RightGroup, views, active, size: DefaultColumnSize },
      ],
    },
  });
}
