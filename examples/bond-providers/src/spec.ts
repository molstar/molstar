/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BallAndStickPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { NaiveDistanceBondProvider, NaiveOrdersBondProvider, RandomBondProvider } from './providers.js';
import { BondOrders } from '@molstar/bond-orders-extension';

/** The example's own bond providers, registered as any registry entry is. */
export const NaiveBonds: PluginRegistryEntry = {
  structure: { bondProviders: [NaiveDistanceBondProvider, NaiveOrdersBondProvider, RandomBondProvider] },
};

/**
 * The default Prism viewer with the example's providers and the bond-orders extension's entry (`BondOrders`: bond-order
 * perception) beside the default registry's (`DefaultBonds`: 'none' and 'per-unit'). The Structure panel's inspector
 * shows a structure's Bonds field; no sequence viewer, no assistant.
 */
export function BondProvidersViewerSpec(): PrismViewerSpec {
  const spec = DefaultPrismViewerSpec();
  return {
    ...spec,
    // the ball-and-stick preset (not in the default registry), to see every bond, the polymers' included
    plugin: {
      ...spec.plugin,
      registry: [...(spec.plugin?.registry ?? []), BallAndStickPresetEntry, NaiveBonds, BondOrders],
    },
    sequence: false,
    agent: false,
  };
}
