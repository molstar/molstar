/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Want } from '@molstar/seqx-core/data';
import { type AccessPlan, inversePath, type StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { dataRegions, type LayerSource, layerLevel } from './layer-data.js';
import type { StoreReader } from './layer-type.js';
import type { SectionLayout } from './tracks-types.js';

// A tooltip table read with what it decorates (_spec/README.md "Render", Tooltip tables): a table no layer draws is
// wanted by no one, so its tiles are wanted beside the marks it can cover, over the same window: each section's axis
// (a strip's, a block's) on its system, or reaching it by a path, carried back onto the table's system, at the level
// a layer there reads. Its wants carry the tooltip's address (the store bounds them together, as a layer's) and are kept
// as the view's are: never evicted while in view. A view with no tooltip table wants nothing more.

/** A table's plan on an axis's system: an assembly's by its sequence (`bySystem`), else its own. */
function planOn(model: StaticModel, address: string, system: string): { plan: AccessPlan; data: string } | undefined {
  const a = model.annotations.get(address);
  if (!a) return undefined;
  if (a.bySystem) {
    const plan = a.bySystem.get(system);
    return plan ? { plan, data: system } : undefined;
  }
  return a.access ? { plan: a.access, data: a.target } : undefined;
}

/** The wants of a view's `ranges` tooltip tables beside its sections' windows, pushed onto `wants`. */
export function tooltipWants(
  model: StaticModel,
  view: ResolvedView,
  data: StoreReader,
  sections: readonly SectionLayout[],
  wants: Want[],
): void {
  const tips = view.tooltips;
  if (!tips) return;
  for (const t of tips) {
    if (t.target === undefined || t.schema !== 'ranges' || t.category) continue;
    for (const s of sections) {
      const st = s.placement.strips;
      const n = st.n || (s.axis ? 1 : 0);
      for (let i = 0; i < n; i++) {
        const axis: AxisTransform = st.n ? st.axes[i] : s.axis!;
        const on = planOn(model, t.address, axis.system);
        if (!on) continue;
        const path = on.data === axis.system ? undefined : model.path(on.data, axis.system).path;
        if (on.data !== axis.system && !path) continue;
        const { plan } = on;
        const src: LayerSource = {
          source: plan.kind === 'inline' ? plan.source.address : plan.source,
          system: plan.system,
          dataSystem: on.data,
          shape: plan.shape,
          ...(path ? { path, back: inversePath(path) } : {}),
          ratio: path?.step === 1 ? 1 / 3 : path?.step === -1 ? 3 : 1,
        };
        const regions = dataRegions(src, axis);
        const level = layerLevel(src, data, axis);
        for (let k = 0; k + 1 < regions.length; k += 2)
          wants.push({
            source: src.source,
            system: src.system,
            start: regions[k],
            end: regions[k + 1],
            level,
            priority: 'view',
            layer: t.address,
          });
      }
    }
  }
}
