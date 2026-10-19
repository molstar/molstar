/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { StructureElement } from '../../element.js';
import type { Structure } from '../../structure.js';
import { InterUnitBonds, type InterUnitEdgeProps } from './data.js';

/**
 * The parent UnitIndex of every element of `child`, ascending; `undefined` if an element is not in `parent`. Takes
 * binary searches, O(c log p), or a merged walk, O(c + p), whichever is cheaper, so that a small child of a large (for
 * example merged) unit does not pay for the whole parent; nothing here or in the callers is sized by the parent.
 */
function getParentIndices(parent: SortedArray<number>, child: ArrayLike<number>): Int32Array | undefined {
  const c = child.length;
  const p = parent.length;
  const ret = new Int32Array(c);
  if (c * Math.log2(p + 1) < p) {
    for (let i = 0, j = 0; i < c; i++) {
      const k = SortedArray.indexOfInRange(parent, child[i], j, p);
      if (k < 0) return undefined;
      ret[i] = k;
      j = k + 1;
    }
  } else {
    for (let i = 0, j = 0; i < c; i++) {
      const e = child[i];
      while (j < p && parent[j] < e) j++;
      if (j === p || parent[j] !== e) return undefined;
      ret[i] = j++;
    }
  }
  return ret;
}

/** The child index of `parentIndex` in `parentIndices` (from `getParentIndices`), or -1. */
function childIndexOf(parentIndices: Int32Array, parentIndex: number) {
  return SortedArray.indexOfInRange(
    parentIndices as ArrayLike<number> as SortedArray,
    parentIndex,
    0,
    parentIndices.length,
  );
}

/**
 * The inter-unit bonds of `child` restricted from `parent`'s, or `undefined` when `child` has a unit that is not a
 * subset of `parent`'s unit with the same id, model and operator (05-model-structure P5). The caller checks that the
 * bond props (`dynamicBonds`, `interBondsValidUnit(Pair)`, `masterModel`) are the parent's.
 *
 * `computeInterUnitBonds` decides mostly per atom pair (thresholds, struct_conn, index pairs, alt-locs), so this
 * equals a recompute on the child except where it decides from atoms the child may lack or from unit sizes:
 * - an atom of the iterated (smaller) unit with a struct_conn bond to the other gets no distance-based bonds there,
 *   so a struct_conn partner outside the child, or a child that swaps which unit is smaller, can change them;
 * - the default pair test classifies water/ion units by their first element, and pairs them only through a common
 *   struct_conn record.
 * In those cases the subset keeps what the parent shows, as it does for bonds the parent carried over from another
 * frame (`remapModel` without `dynamicBonds`). The edge order follows the parent's units, not a recompute's.
 */
export function subsetInterUnitBonds(bonds: InterUnitBonds, parent: Structure, child: Structure) {
  // the parent indices of a child unit's elements, or null when it has every element of the parent unit
  const maps = new Map<number, Int32Array | null>();
  let whole = child.units.length === parent.units.length;
  for (const u of child.units) {
    const p = parent.unitMap.get(u.id);
    if (!p || p.kind !== u.kind || p.model !== u.model) return;
    const opP = p.conformation.operator;
    const opU = u.conformation.operator;
    // as `Unit.areOperatorsEqual`, which this module cannot import (unit.ts imports it)
    if (opP !== opU && !Mat4.areEqual(opP.matrix, opU.matrix, 1e-6)) return;
    if (p.elements === u.elements || p.elements.length === u.elements.length) {
      if (p.elements !== u.elements && !SortedArray.areEqual(p.elements, u.elements)) return;
      maps.set(u.id, null);
      continue;
    }
    const parentIndices = getParentIndices(p.elements, u.elements);
    if (!parentIndices) return;
    maps.set(u.id, parentIndices);
    whole = false;
  }
  if (whole) return bonds;

  type UnitIndex = StructureElement.UnitIndex;
  const builder = new InterUnitGraph.Builder<number, UnitIndex, InterUnitEdgeProps>();
  for (const u of child.units) {
    const mapA = maps.get(u.id)!;
    for (const pair of bonds.getConnectedUnits(u.id)) {
      // every pair is listed under both of its units; take it once
      if (pair.unitB < pair.unitA) continue;
      const mapB = maps.get(pair.unitB);
      if (mapB === undefined) continue;

      builder.startUnitPair(pair.unitA, pair.unitB);
      for (const indexA of pair.connectedIndices) {
        const a = mapA ? childIndexOf(mapA, indexA) : indexA;
        if (a < 0) continue;
        for (const { indexB, props } of pair.getEdges(indexA)) {
          const b = mapB ? childIndexOf(mapB, indexB) : indexB;
          if (b >= 0) builder.add(a as UnitIndex, b as UnitIndex, props);
        }
      }
      builder.finishUnitPair();
    }
  }
  return new InterUnitBonds(builder.getMap());
}
