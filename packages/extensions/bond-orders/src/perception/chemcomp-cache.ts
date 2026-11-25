/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { getFlags, getOrder, isPerceivable, type State, setBond } from './common.js';

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

/** A residue already hashed: its atom names as read (file order, altloc duplicates included) and its key. */
interface HashedResidue {
  readonly names: readonly string[];
  readonly hash: string;
}

/** Model-wide, by mode prefix, component id and atom count. */
const HashedResidues = new WeakMap<Model, Map<string, HashedResidue[]>>();

/**
 * The key of a residue's pattern (_spec §1.5): the mode prefix, the component id and its sorted unique atom names (an
 * incomplete residue gets its own key). A residue whose atom names, in order, are those of one already hashed gets
 * its key without sorting: same names in the same order, same unique names.
 */
export function chemCompHash(unit: Unit.Atomic, start: UnitIndex, end: UnitIndex, prefix: string) {
  const { elements, model } = unit;
  const { label_atom_id, label_comp_id } = model.atomicHierarchy.atoms;
  let byKey = HashedResidues.get(model);
  if (!byKey) {
    byKey = new Map();
    HashedResidues.set(model, byKey);
  }
  const compId = label_comp_id.value(elements[start]);
  const key = `${prefix}${compId}|${end - start}`;
  let hashed = byKey.get(key);
  if (!hashed) {
    hashed = [];
    byKey.set(key, hashed);
  }
  search: for (const { names, hash } of hashed) {
    for (let u = start, k = 0; u < end; u++, k++) {
      if (label_atom_id.value(elements[u]) !== names[k]) continue search;
    }
    return hash;
  }

  const names: string[] = [];
  for (let u = start; u < end; u++) names.push(label_atom_id.value(elements[u]));
  const hash = `${prefix}${compId}|${[...new Set(names)].sort().join(',')}`;
  hashed.push({ names, hash });
  return hash;
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

/**
 * Apply the cached pattern of `hash` onto the residue's perceivable bonds, matched by atom names; false if none is
 * cached. An empty pattern is a result (a component without multiple or aromatic bonds): applied, it writes nothing.
 * Only bonds inside the residue, as the pattern records them: a bond leaving the residue is never written (_spec
 * §1.11), even when its atom names match a pattern's pair.
 */
export function applyCachedChemCompPattern(
  unit: Unit.Atomic,
  bonds: IntraUnitBonds,
  start: UnitIndex,
  end: UnitIndex,
  hash: string,
): boolean {
  const pattern = getModelCache(unit.model).get(hash);
  if (!pattern) return false;
  if (pattern.size === 0) return true; // only single bonds, nothing to apply

  const { elements } = unit;
  const { label_atom_id } = unit.model.atomicHierarchy.atoms;
  const { offset, b, edgeProps } = bonds;
  for (let u = start; u < end; u++) {
    const nameA = label_atom_id.value(elements[u]);
    for (let t = offset[u]; t < offset[u + 1]; t++) {
      const v = b[t];
      // `u < v` and `u >= start`: `v < end` keeps the bond inside the residue
      if (u > v || v >= end || !BondType.isCovalent(edgeProps.flags[t])) continue;
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
