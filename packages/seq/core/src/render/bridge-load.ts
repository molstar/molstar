/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bridge chunk's door (_spec/README.md "Render", Bridges: lazy): `bridges.ts` and what it imports (the bridge's
// layout, the ribbons) load the first time a spec holds a bridge, a model waiting for it as for the genome chunk
// (`withViewKind`), so a view with a bridge is never drawn without it. Why lazy: a view without a bridge pays nothing.

import { withViewKind } from '@molstar/seqx-core/resolve';
import { BridgePass, type ResolvedView } from '@molstar/seqx-core/view';
import { BridgeHooks } from './bridge-hooks.js';
import type { RenderConfig } from './config.js';

/** Loads the bridge chunk and the view's bridge pass. */
export function loadBridges(): Promise<void> {
  return Promise.all([
    import('./bridges.js').then((m) => {
      BridgeHooks.layout = m.layoutBridges;
    }),
    import('@molstar/seqx-core/view/bridges').then((m) => {
      BridgePass.resolve = m.resolveBridges;
    }),
  ]).then(() => undefined);
}

/**
 * Registers bridges with the layout and resolve (the pairs type's module calls it): their room, the pairs layer's region
 * reader, and the chunks loaded when a spec holds a bridge.
 */
export function registerBridges(region: NonNullable<(typeof BridgeHooks)['region']>): void {
  if (BridgeHooks.rooms) return;
  BridgeHooks.rooms = bridgeRooms;
  BridgeHooks.region = region;
  withViewKind('bridge', loadBridges);
}

/**
 * Each placed, shown bridge's room by the section below it (_spec/README.md "Render", Bridges): its height and a
 * section gap, taken above that section.
 */
export function bridgeRooms(view: ResolvedView, config: RenderConfig): Map<string, number> {
  const out = new Map<string, number>();
  const { stack = [], bridges = [], sections } = view;
  for (let i = 0; i + 1 < stack.length; i++) {
    const e = stack[i];
    if (e.kind !== 'bridge' || stack[i + 1].kind !== 'section') continue;
    const br = bridges[e.index];
    if (!br.placed || br.hidden) continue;
    out.set(sections[stack[i + 1].index].address, (br.height ?? config.bridges.height) + config.sectionGap);
  }
  return out;
}
