/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import type { Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { getFlags, getOrder, isPerceivable, type State, setBond } from './common.js';
import type { InterUnitNeighbours } from './inter-unit-neighbours.js';

/** The perceived multiple and aromatic bonds of a component, by atom-name pair. */
type BondingPattern = Map<string, { order: number; flags: number }>;

/** Model-wide, keyed by `chemCompHash` (_spec §1.5). */
const PerceivedCache = new WeakMap<Model, Map<string, BondingPattern>>();

function getModelCache(model: Model) {
  let c = PerceivedCache.get(model);
  if (!c) {
    c = new Map();
    PerceivedCache.set(model, c);
  }
  return c;
}

/**
 * The key of a residue's pattern (_spec §1.5): the mode prefix, the component id, its sorted unique atom names (an
 * incomplete residue gets its own key) and its links, the residue's atoms bonded outside the residue with their
 * partner (`C12-SER.OG`), in this unit or in another unit of the model. A pattern is only reused for copies with the
 * same atoms and the same links.
 */
export function chemCompHash(
  unit: Unit.Atomic,
  bonds: IntraUnitBonds,
  interUnit: InterUnitNeighbours,
  start: UnitIndex,
  end: UnitIndex,
  prefix: string,
) {
  const { elements } = unit;
  const { label_atom_id, label_comp_id } = unit.model.atomicHierarchy.atoms;
  const { offset, b, edgeProps } = bonds;
  const names = new Set<string>();
  const linkSet = new Set<string>();
  for (let u = start; u < end; u++) {
    const name = label_atom_id.value(elements[u]);
    names.add(name);
    for (let t = offset[u]; t < offset[u + 1]; t++) {
      const v = b[t];
      if ((v >= start && v < end) || !BondType.isCovalent(edgeProps.flags[t])) continue;
      linkSet.add(`${name}-${label_comp_id.value(elements[v])}.${label_atom_id.value(elements[v])}`);
    }
    const partners = interUnit.get(elements[u]);
    if (!partners) continue;
    for (const eI of partners) {
      if (SortedArray.has(elements, eI)) continue; // an intra-unit bond, listed above
      linkSet.add(`${name}-${label_comp_id.value(eI)}.${label_atom_id.value(eI)}`);
    }
  }
  const compId = label_comp_id.value(elements[start]);
  return `${prefix}${compId}|${[...names].sort().join(',')}|${[...linkSet].sort().join(',')}`;
}

function extractPattern(state: State): BondingPattern {
  const { unit, unitIndices, n, bonds } = state;
  const { label_atom_id } = unit.model.atomicHierarchy.atoms;
  const pattern: BondingPattern = new Map();
  for (let i = 0; i < n; i++) {
    const u = unitIndices[i];
    for (const j of state.heavyNeighbours[i]) {
      if (i >= j) continue;
      const v = unitIndices[j];
      const order = getOrder(bonds, u, v);
      const flags = getFlags(bonds, u, v);
      if (order > 1 || flags & BondType.Flag.Aromatic) {
        const a = label_atom_id.value(unit.elements[u]);
        const b = label_atom_id.value(unit.elements[v]);
        const key = a < b ? `${a}|${b}` : `${b}|${a}`;
        pattern.set(key, { order, flags: flags & BondType.Flag.Aromatic });
      }
    }
  }
  return pattern;
}

/** Apply the cached pattern of `hash` onto the residue's perceivable bonds, matched by atom names; false if none. */
export function applyCachedChemCompPattern(
  unit: Unit.Atomic,
  bonds: IntraUnitBonds,
  start: UnitIndex,
  end: UnitIndex,
  hash: string,
): boolean {
  const pattern = getModelCache(unit.model).get(hash);
  if (!pattern || pattern.size === 0) return false;

  const { elements } = unit;
  const { label_atom_id } = unit.model.atomicHierarchy.atoms;
  const { offset, b, edgeProps } = bonds;
  for (let u = start; u < end; u++) {
    const nameA = label_atom_id.value(elements[u]);
    for (let t = offset[u]; t < offset[u + 1]; t++) {
      const v = b[t] as UnitIndex;
      if (u > v || !BondType.isCovalent(edgeProps.flags[t])) continue;
      const nameB = label_atom_id.value(elements[v]);
      const p = pattern.get(nameA < nameB ? `${nameA}|${nameB}` : `${nameB}|${nameA}`);
      if (!p || !isPerceivable(getFlags(bonds, u, v), getOrder(bonds, u, v))) continue;
      setBond(bonds, u, v, p.order, p.flags);
    }
  }
  return true;
}

export function cacheChemCompPattern(state: State, hash: string) {
  const cache = getModelCache(state.unit.model);
  if (cache.has(hash)) return;
  cache.set(hash, extractPattern(state));
}
