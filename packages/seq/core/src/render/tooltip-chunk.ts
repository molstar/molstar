/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { ItemTitle } from './accessibility.js';
import type { StoreReader } from './layer-type.js';
import { type TooltipOut, tooltipsAt } from './tooltips.js';

// Render's tooltip chunk (`tooltip-hooks.ts`): what loads with the tooltips, one dynamic import.

export { tooltipsAt } from './tooltips.js';
export { tooltipSources, tooltipWants } from './tooltip-wants.js';

/** An export's `card` titles from the spec's tooltips alone (a headless export: no viewer's providers). */
export function tooltipTitles(view: ResolvedView, model: StaticModel, data: StoreReader): ItemTitle {
  return (_, l, it, name) => {
    const system = it.system ?? l.layer.data?.system;
    if (system === undefined) return undefined;
    const out: TooltipOut = { rows: [], more: 0, failed: [] };
    tooltipsAt(
      { model, view, data },
      { layer: l.layer, item: { key: it.key, start: it.start }, parts: [{ system, start: it.start, end: it.end }] },
      out,
    );
    const rows = out.rows.filter((r) => !r.loading).map((r) => (r.value ? `${r.label}: ${r.value}` : r.label));
    return rows.length ? [name, ...rows].join('\n') : undefined;
  };
}
