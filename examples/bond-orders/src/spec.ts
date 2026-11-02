/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { BondOrders } from '@molstar/bond-orders-extension';
import { BallAndStickPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';

/**
 * The default Prism viewer with the bond-orders extension's entry: bond-order perception is one more option of a
 * structure's Bonds field (Structure panel, pick the structure). No sequence viewer, no assistant.
 */
export function BondOrdersViewerSpec(): PrismViewerSpec {
  const spec = DefaultPrismViewerSpec();
  return {
    ...spec,
    // ball-and-stick everywhere (not in the default registry), so the bond orders are what you see
    plugin: { ...spec.plugin, registry: [...(spec.plugin?.registry ?? []), BallAndStickPresetEntry, BondOrders] },
    sequence: false,
    agent: false,
  };
}
