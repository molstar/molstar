/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the layout reads of bridges (_spec/README.md "Render", Bridges): set by `bridge-load.ts`, which the pairs type
// registers, so a composition without the pairs type draws no bridge and holds none of their code (its sides stacked,
// as a reader without bridges draws them).

import type { ResolvedView } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';

export const BridgeHooks: {
  /** Each placed, shown bridge's room by the section below it. */
  rooms?: (view: ResolvedView, config: RenderConfig) => Map<string, number>;
  /** The chunk's entry once landed (`layoutBridges`). */
  layout?: typeof import('./bridges.js').layoutBridges;
  /** The pairs layer's region reader (a bridge reads a source's pairs as that layer does). */
  region?: typeof import('./pairs-region.js').regionPairs;
  /** A region read's kept state (its pending tiles: the pending hatch). */
  state?: typeof import('./pairs-region.js').regionState;
} = {};
