/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { StructureElement } from '../../element.js';
import type { Structure } from '../../structure.js';
import { InterUnitBonds, type InterUnitEdgeProps, IntraUnitBonds } from './data.js';

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
 * The intra-unit bonds of a child unit with `elements`, restricted from `bonds` of its parent unit with
 * `parentElements` (05-model-structure P6); `undefined` if `elements` is not a subset of `parentElements`.
 *
 * The edges equal a recompute on the child: `computeIntraUnitBonds` decides per atom pair (struct_conn, the component
 * dictionary, the hydrogen fallback, alt-locs, thresholds, index pairs). Their order is the parent's. The props
 * (`canRemap`, `cacheable`) are the parent's: a parent that may remap its bonds has children that may too, while a
 * child of one that may not could have been allowed (it then recomputes on a remap). As in the parent, bonds carried
 * over from another frame (`remapModel`) are carried over to the child.
 */
export function subsetIntraUnitBonds(
  bonds: IntraUnitBonds,
  parentElements: ArrayLike<number>,
  elements: ArrayLike<number>,
): IntraUnitBonds | undefined {
  if (elements.length <= 1) return IntraUnitBonds.Empty;
  // the parent UnitIndex of every child element
  const vertices = new Int32Array(elements.length);
  for (let i = 0, j = 0, il = elements.length, jl = parentElements.length; i < il; i++) {
    const e = elements[i];
    while (j < jl && parentElements[j] < e) j++;
    if (j === jl || parentElements[j] !== e) return undefined;
    vertices[i] = j++;
  }
  return IntAdjacencyGraph.induceByVertices(bonds, vertices, bonds.props) as IntraUnitBonds;
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
  // null: the child unit has every element of the parent unit
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
