/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BallAndStickPresetEntry } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { NaiveDistanceBondProvider, NaiveOrdersBondProvider } from './providers.js';

/** The example's own bond providers, registered as any registry entry is. */
export const NaiveBonds: PluginRegistryEntry = {
  structure: { bondProviders: [NaiveDistanceBondProvider, NaiveOrdersBondProvider] },
};

/**
 * The default Prism viewer with the example's providers beside the default registry's (`DefaultBonds`: 'none' and
 * 'per-unit'). The Structure panel's inspector shows a structure's Bonds field; no sequence viewer, no assistant.
 */
export function BondProvidersViewerSpec(): PrismViewerSpec {
  const spec = DefaultPrismViewerSpec();
  return {
    ...spec,
    // ball-and-stick everywhere (not in the default registry), so the bonds are what you see
    plugin: { ...spec.plugin, registry: [...(spec.plugin?.registry ?? []), BallAndStickPresetEntry, NaiveBonds] },
    sequence: false,
    agent: false,
  };
}
