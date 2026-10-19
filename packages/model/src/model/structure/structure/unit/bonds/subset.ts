/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import type { StructureElement } from '../../element.js';
import type { Structure } from '../../structure.js';
import { Unit } from '../../unit.js';
import { InterUnitBonds, type InterUnitEdgeProps } from './data.js';

/**
 * Maps the UnitIndices of `parent` to those of `child`, whose elements must be a subset of the parent's:
 * `map[parentIndex]` is the child index or -1. `undefined` if an element of `child` is not in `parent`.
 */
function getUnitIndexMap(parent: ArrayLike<number>, child: ArrayLike<number>): Int32Array | undefined {
  const map = new Int32Array(parent.length).fill(-1);
  for (let i = 0, j = 0, il = child.length, jl = parent.length; i < il; i++) {
    const e = child[i];
    while (j < jl && parent[j] < e) j++;
    if (j === jl || parent[j] !== e) return undefined;
    map[j++] = i;
  }
  return map;
}

/**
 * The inter-unit bonds of `child` restricted from `parent`'s, or `undefined` when `child` has a unit that is not a
 * subset of `parent`'s unit with the same id, model and operator (05-model-structure P5). The caller checks that the
 * bond props (`dynamicBonds`, `interBondsValidUnit(Pair)`, `masterModel`) are the parent's.
 *
 * `computeInterUnitBonds` decides per atom pair (thresholds, struct_conn, index pairs, alt-locs), so this equals a
 * recompute on the child except where it decides from atoms the child may lack: an atom with a struct_conn bond to
 * another unit gets no distance-based bonds there, and water/ion units pair only through a common struct_conn record.
 * When those partners are outside the child, a recompute could differ; the subset keeps what the parent shows.
 */
export function subsetInterUnitBonds(bonds: InterUnitBonds, parent: Structure, child: Structure) {
  // null: the child unit has every element of the parent unit
  const maps = new Map<number, Int32Array | null>();
  let whole = child.units.length === parent.units.length;
  for (const u of child.units) {
    const p = parent.unitMap.get(u.id);
    if (!p || p.kind !== u.kind || p.model !== u.model) return;
    if (p.conformation.operator !== u.conformation.operator && !Unit.areOperatorsEqual(p, u)) return;
    if (p.elements === u.elements || p.elements.length === u.elements.length) {
      if (p.elements !== u.elements && !SortedArray.areEqual(p.elements, u.elements)) return;
      maps.set(u.id, null);
      continue;
    }
    const map = getUnitIndexMap(p.elements, u.elements);
    if (!map) return;
    maps.set(u.id, map);
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
        const a = mapA ? mapA[indexA] : indexA;
        if (a < 0) continue;
        for (const { indexB, props } of pair.getEdges(indexA)) {
          const b = mapB ? mapB[indexB] : indexB;
          if (b >= 0) builder.add(a as UnitIndex, b as UnitIndex, props);
        }
      }
      builder.finishUnitPair();
    }
  }
  return new InterUnitBonds(builder.getMap());
}
