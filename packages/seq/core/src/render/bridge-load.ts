/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bridge chunk's door (_spec/README.md "Render", Bridges: lazy): `bridges.ts` and what it imports (the bridge's
// layout, the ribbons) load the first time a spec holds a bridge, a model waiting for it as for the genome chunk
// (`withViewKind`), so a view with a bridge is never drawn without it. Why lazy: a view without a bridge pays nothing.

import { withViewKind } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';

/** The chunk's entry once landed (`layoutBridges`), else empty: the layout reserves the bands and draws no ribbon. */
export const BridgeChunk: { layout?: typeof import('./bridges.js').layoutBridges } = {};

/**
 * The pairs layer's reads a bridge reuses (`pairs.ts` sets them where a composition has the pairs type): through this
 * record, so the chunk imports none of the pairs layer's modules (a module the entry shares would load with it).
 */
export const PairsReads: { region?: typeof import('./pairs-region.js').regionPairs } = {};

/** Loads the bridge chunk (once while in flight). */
export function loadBridges(): Promise<void> {
  return import('./bridges.js').then((m) => {
    BridgeChunk.layout = m.layoutBridges;
  });
}

withViewKind('bridge', loadBridges);

/**
 * Each placed, shown bridge's room by the section below it (_spec/README.md "Render", Bridges): its height and a
 * section gap, taken above that section.
 */
export function bridgeRooms(view: ResolvedView, config: RenderConfig): Map<string, number> {
  const out = new Map<string, number>();
  const { stack, bridges, sections } = view;
  for (let i = 0; i + 1 < stack.length; i++) {
    const e = stack[i];
    if (e.kind !== 'bridge' || stack[i + 1].kind !== 'section') continue;
    const br = bridges[e.index];
    if (!br.placed || br.hidden) continue;
    out.set(sections[stack[i + 1].index].address, (br.height ?? config.bridges.height) + config.sectionGap);
  }
  return out;
}
