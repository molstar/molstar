/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { withTooltips } from '@molstar/seqx-core/resolve';
import type * as Chunk from './tooltip-chunk.js';

// Render's tooltips behind the door (_spec/README.md "Render", Tooltips, Loaded when used): the lookup, the tables'
// wants and sources, the export's titles from the spec's tooltips load with the tooltips (resolve's `withTooltips`),
// before a model holding a tooltip settles. Until then a view has no tooltips (its pass has not landed either).

/** Render's tooltip chunk, when it has landed. @internal */
export const TooltipRender: { chunk?: typeof Chunk } = {};

withTooltips(() =>
  import('./tooltip-chunk.js').then((m) => {
    TooltipRender.chunk = m;
  }),
);
