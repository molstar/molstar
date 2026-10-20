/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { StructureElement, Unit } from '@molstar/model/model/structure';
import { InterUnitBonds, type IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import type { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';

// Three naive bond providers, written as a third party would: against the public interface only (`BondProvider`, the
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

/**
 * Random bonds, for fun and for testing what draws them: each pair of atoms in a unit closer than `maxDistance` is
 * bonded with probability `probability`, with a random order from 1 to 3. Random but deterministic: the choice is a
 * hash of the seed and the two atoms' model indices, so the same params give the same bonds (the provider contract:
 * static bonds, a sound cache key), and a child unit picks exactly the pairs its parent did. Between units, none.
 */
export const RandomBondProvider: BondProvider<{
  probability: PD.Numeric;
  maxDistance: PD.Numeric;
  seed: PD.Numeric;
}> = {
  name: 'random',
  label: 'Random',
  description: 'Random bonds between nearby atoms, the same for the same seed',
  getParams: () => ({
    probability: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, { label: 'Probability' }),
    maxDistance: PD.Numeric(4, { min: 1, max: 8, step: 0.5 }, { label: 'Max distance (Å)' }),
    seed: PD.Numeric(42, { min: 0, max: 1000, step: 1 }, { label: 'Seed' }),
  }),
  create: ({ probability, maxDistance, seed }) => ({
    provider: RandomBondProvider,
    cacheKey: `random:${probability}:${maxDistance}:${seed}`,
    intra: (unit) => randomBonds(unit, probability, maxDistance, seed),
    inter: () => new InterUnitBonds(new Map()),
  }),
};

/**
 * An integer mixer with uniform output bits (lowbias32). Core's `hash1`…`hash4` are for bucketing: their top bits
 * barely change for small inputs, so they make poor random numbers.
 */
function mix(x: number) {
  x ^= x >>> 16;
  x = Math.imul(x, 0x7feb352d);
  x ^= x >>> 15;
  x = Math.imul(x, 0x846ca68b);
  x ^= x >>> 16;
  return x >>> 0;
}

/** A number in [0, 1) from the seed and an atom pair, the same for the same inputs. */
function pairRandom(seed: number, a: number, b: number, salt: number) {
  return mix(mix(mix(seed * 2 + salt) ^ a) ^ b) / 4294967296;
}

function randomBonds(unit: Unit.Atomic, probability: number, maxDistance: number, seed: number): IntraUnitBonds {
  const { elements, lookup3d } = unit;
  const { x, y, z } = unit.model.atomicConformation;
  const a: StructureElement.UnitIndex[] = [];
  const b: StructureElement.UnitIndex[] = [];
  const orders: number[] = [];
  for (let i = 0; i < elements.length; i++) {
    const e = elements[i];
    const { indices, count } = lookup3d.find(x[e], y[e], z[e], maxDistance);
    for (let k = 0; k < count; k++) {
      const j = indices[k];
      if (j <= i) continue;
      if (pairRandom(seed, e, elements[j], 0) >= probability) continue;
      a.push(i as StructureElement.UnitIndex);
      b.push(j);
      orders.push(1 + Math.floor(pairRandom(seed, e, elements[j], 1) * 3));
    }
  }
  const builder = new IntAdjacencyGraph.EdgeBuilder(elements.length, a, b);
  const order = new Int8Array(builder.slotCount);
  const flags = new Uint16Array(builder.slotCount);
  const key = new Int32Array(builder.slotCount);
  for (let i = 0; i < builder.edgeCount; i++) {
    builder.addNextEdge();
    builder.assignProperty(order, orders[i]);
    builder.assignProperty(flags, BondType.Flag.Covalent | BondType.Flag.Computed);
    builder.assignProperty(key, -1);
  }
  return builder.createGraph({ order, flags, key }, { canRemap: false, cacheable: true });
}
