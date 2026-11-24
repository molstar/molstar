/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's bridge parts' door (_spec/README.md "Bridges"): the chunk (`bridge-view.ts`: the lock's carrying, the
// band's gestures, `bridge.align`, `bridge.selectAcross`, the lock toggle, the cursor across a band, the card) loads
// when a spec holds a bridge, as render's bridge chunk does (core's `withViewKind`: the model waits for both), and sets
// `BridgeView.parts`; the view model, the input, the cursor, the describer and the viewer call it, one optional call
// each. Why: a viewer that never shows a bridge pays for none of it.

import { withViewKind } from '@molstar/seqx-core/resolve';

/** What the bridge chunk gives the eager viewer (`bridge-view.ts`'s `BridgeViewParts`). */
export const BridgeView: { parts?: typeof import('../bridge-view.js').BridgeViewParts } = {};

withViewKind('bridge', () =>
  import('../bridge-view.js').then((m) => {
    BridgeView.parts = m.BridgeViewParts;
  }),
);
