/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BondOrderPerceptionProvider } from './provider.js';

export { type BondOrdersParams, BondOrderPerceptionProvider } from './provider.js';
export type { BondOrdersMode } from './perception/perceive.js';

/** Add to a plugin spec's `registry` to offer bond-order perception in a structure's `bonds` param (_spec §2). */
export const BondOrders: PluginRegistryEntry = {
  structure: { bondProviders: [BondOrderPerceptionProvider] },
};
