/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';
import { getInterUnitNeighbours } from './perception/inter-unit-neighbours.js';
import { type BondOrdersMode, perceiveBondOrders } from './perception/perceive.js';

const BondOrdersParams = {
  mode: PD.Select<BondOrdersMode>(
    'auto',
    [
      ['auto', 'Auto'],
      ['forceCompute', 'Force compute'],
    ],
    {
      description:
        'Auto: perceive the bond orders missing from the file and Mol* dictionary for canonical residue. Force compute: perceive every bond.',
    },
  ),
};
export type BondOrdersParams = typeof BondOrdersParams;

/**
 * The built-in bonds with perceived bond orders (_spec §2): refines `ctx.builtIn()`, intra-unit only (the bonds
 * between units are the built-in ones).
 */
export const BondOrderPerceptionProvider: BondProvider<BondOrdersParams> = {
  name: 'bond-order-perception',
  label: 'Bond order perception',
  description: 'The built-in bond graph, with the bond orders perceived from the 3D coordinates',
  getParams: () => BondOrdersParams,
  create: ({ mode }) => instance(mode),
};

/** One instance per mode: equal props, equal instance (a child unit restricts its parent's bonds of the same one). */
const Instances = new Map<BondOrdersMode, BondProvider.Instance>();

function instance(mode: BondOrdersMode): BondProvider.Instance {
  let i = Instances.get(mode);
  if (i) return i;
  i = {
    provider: BondOrderPerceptionProvider,
    cacheKey: mode,
    // perception spans rings and conjugated systems: a child gets the root unit's orders (_spec §1.8)
    childBonds: 'root',
    intra: (unit, ctx) =>
      perceiveBondOrders(
        unit,
        ctx.builtIn(),
        () => getInterUnitNeighbours(unit.model, () => ctx.modelStructure()),
        mode,
      ),
  };
  Instances.set(mode, i);
  return i;
}
