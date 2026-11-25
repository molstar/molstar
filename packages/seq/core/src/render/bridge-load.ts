/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bridge chunk's door (_spec/README.md "Render", Bridges: lazy): `bridges.ts` and what it imports (the bridge's
// layout, the ribbons) load the first time a spec holds a bridge, a model waiting for it as for the genome chunk
// (`withViewKind`), so a view with a bridge is never drawn without it. Why lazy: a view without a bridge pays nothing.

import { withViewKind } from '@molstar/seqx-core/resolve';
import { BridgePass } from '@molstar/seqx-core/view';
import { BridgeHooks } from './bridge-hooks.js';

/** Loads the bridge chunk and the view's bridge pass. */
export function loadBridges(): Promise<void> {
  return Promise.all([
    import('./bridges.js').then((m) => {
      BridgeHooks.rooms = m.bridgeRooms;
      BridgeHooks.layout = m.layoutBridges;
    }),
    import('@molstar/seqx-core/view/bridges').then((m) => {
      BridgePass.resolve = m.resolveBridges;
    }),
  ]).then(() => undefined);
}

/**
 * Registers bridges with the layout and resolve (the pairs type's module calls it): the pairs layer's region reader,
 * and the chunks loaded when a spec holds a bridge (the bands' room is the chunk's: a view with a bridge waits for it).
 */
export function registerBridges(
  region: NonNullable<(typeof BridgeHooks)['region']>,
  state: NonNullable<(typeof BridgeHooks)['state']>,
): void {
  if (BridgeHooks.region) return;
  BridgeHooks.region = region;
  BridgeHooks.state = state;
  withViewKind('bridge', loadBridges);
}
