/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { StructureElement, Unit } from '@molstar/model/model/structure';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import type { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';

// Two naive bond providers, written as a third party would: against the public interface only (`BondProvider`, the
// unit, the model, a graph builder). They show the two ways a provider works; neither is meant to be right.

/**
 * Replaces the bonds: every pair of atoms in a unit closer than `cutoff` is bonded, ignoring the component dictionary,
 * struct_conn and the element thresholds. Between units it takes the built-in bonds, but pairs waters and ions too.
 */
export const NaiveDistanceBondProvider: BondProvider<{ cutoff: PD.Numeric }> = {
  name: 'naive-distance',
  label: 'Naive distance',
  description: 'Every pair of atoms closer than the cutoff, nothing else',
  getParams: () => ({ cutoff: PD.Numeric(1.9, { min: 1, max: 3, step: 0.05 }, { label: 'Cutoff (Å)' }) }),
  create: ({ cutoff }) => ({
    provider: NaiveDistanceBondProvider,
    // equal cutoffs give equal bonds: a model's units are computed once per cutoff
    cacheKey: `naive-distance:${cutoff}`,
    intra: (unit) => bondsWithin(unit, cutoff),
    inter: (_, ctx) => ctx.builtIn({ ignoreWater: false, ignoreIon: false }),
  }),
};

function bondsWithin(unit: Unit.Atomic, cutoff: number): IntraUnitBonds {
  const { elements, lookup3d } = unit;
  const { x, y, z } = unit.model.atomicConformation;
  const a: StructureElement.UnitIndex[] = [];
  const b: StructureElement.UnitIndex[] = [];
  for (let i = 0; i < elements.length; i++) {
    const e = elements[i];
    const { indices, count } = lookup3d.find(x[e], y[e], z[e], cutoff);
    for (let k = 0; k < count; k++) {
      if (indices[k] <= i) continue;
      a.push(i as StructureElement.UnitIndex);
      b.push(indices[k]);
    }
  }
  const builder = new IntAdjacencyGraph.EdgeBuilder(elements.length, a, b);
  const order = new Int8Array(builder.slotCount);
  const flags = new Uint16Array(builder.slotCount);
  const key = new Int32Array(builder.slotCount);
  for (let i = 0; i < builder.edgeCount; i++) {
    builder.addNextEdge();
    builder.assignProperty(order, 1);
    builder.assignProperty(flags, BondType.Flag.Covalent | BondType.Flag.Computed);
    builder.assignProperty(key, -1);
  }
  // the bonds depend on the coordinates: recomputed for another frame when bonds are dynamic
  return builder.createGraph({ order, flags, key }, { canRemap: false, cacheable: true });
}

/**
 * Refines the built-in bonds: single bonds between two carbons shorter than 1.42 Å become double. Only the intra-unit
 * bonds; between units the built-in ones are kept (no `inter`).
 */
export const NaiveOrdersBondProvider: BondProvider<{}> = {
  name: 'naive-orders',
  label: 'Naive bond orders',
  description: 'The built-in bonds, with short C–C single bonds made double',
  getParams: () => ({}),
  create: () => NaiveOrdersInstance,
};

const NaiveOrdersInstance: BondProvider.Instance = {
  provider: NaiveOrdersBondProvider,
  cacheKey: 'naive-orders',
  intra: (unit, ctx) => withDoubleBonds(unit, ctx.builtIn()),
};

const MaxDoubleBondLength = 1.42;

function withDoubleBonds(unit: Unit.Atomic, bonds: IntraUnitBonds): IntraUnitBonds {
  const { elements } = unit;
  const { x, y, z } = unit.model.atomicConformation;
  const { type_symbol } = unit.model.atomicHierarchy.atoms;
  const { offset, a, b, edgeCount, edgeProps } = bonds;
  // both directions of an edge are stored, each with its own slot: the same answer for both
  const order = Int8Array.from(edgeProps.order);
  for (let i = 0; i < a.length; i++) {
    if (order[i] !== 1 || !BondType.isCovalent(edgeProps.flags[i])) continue;
    const ea = elements[a[i]];
    const eb = elements[b[i]];
    if (type_symbol.value(ea) !== 'C' || type_symbol.value(eb) !== 'C') continue;
    const d = Math.hypot(x[ea] - x[eb], y[ea] - y[eb], z[ea] - z[eb]);
    if (d < MaxDoubleBondLength) order[i] = 2;
  }
  return IntAdjacencyGraph.create(offset, a, b, edgeCount, { ...edgeProps, order }, bonds.props);
}
