/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { arrayPickIndices } from '@molstar/core/data/util';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import type { StructureElement } from '../../element.js';
import type { Structure } from '../../structure.js';
import { InterUnitBonds, type InterUnitEdgeProps, IntraUnitBonds } from './data.js';

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
 * The intra-unit bonds of a child unit with `elements`, restricted from `bonds` of its parent unit with
 * `parentElements` (05-model-structure P6, unit.md §4); `undefined` if `elements` is not a subset of `parentElements`.
 *
 * The edges equal a recompute on the child: `computeIntraUnitBonds` decides per atom pair (struct_conn, the component
 * dictionary, the hydrogen fallback, alt-locs, thresholds, index pairs), and from the unit's operator only for index
 * pairs with operators, where the caller restricts a parent of the same operator only. Their order is the parent's. The props
 * (`canRemap`, `cacheable`) are the parent's: a parent that may remap its bonds has children that may too, while a
 * child of one that may not could have been allowed (it then recomputes on a remap). As in the parent, bonds carried
 * over from another frame (`remapModel`) are carried over to the child.
 */
export function subsetIntraUnitBonds(
  bonds: IntraUnitBonds,
  parentElements: SortedArray<number>,
  elements: ArrayLike<number>,
): IntraUnitBonds | undefined {
  const n = elements.length;
  if (n <= 1) return IntraUnitBonds.Empty;
  const vertices = getParentIndices(parentElements, elements);
  if (!vertices) return undefined;
  const { order, flags, key } = bonds.edgeProps;
  // no edges: an empty graph of the child's size with the parent's props, as a recompute gives (the parent's graph
  // may be the shared `IntraUnitBonds.Empty`, which has no offsets)
  if (bonds.edgeCount === 0) {
    const none = new Int32Array(0);
    return IntAdjacencyGraph.create(
      new Int32Array(n + 1),
      none as ArrayLike<number> as ArrayLike<StructureElement.UnitIndex>,
      none as ArrayLike<number> as ArrayLike<StructureElement.UnitIndex>,
      0,
      { order: arrayPickIndices(order, none), flags: arrayPickIndices(flags, none), key: arrayPickIndices(key, none) },
      bonds.props,
    );
  }

  // the induced subgraph, visiting only the child's vertices and their parent edges (IntAdjacencyGraph
  // .induceByVertices visits every parent vertex)
  const { offset, b } = bonds;
  // a child that is a large share of its parent looks its neighbours up in a parent-sized map, which costs no more
  // than the walk of `getParentIndices` did; a small one searches its own vertices
  const p = parentElements.length;
  let map: Int32Array | undefined;
  if (n * Math.log2(p + 1) >= p) {
    map = new Int32Array(p).fill(-1);
    for (let i = 0; i < n; i++) map[vertices[i]] = i;
  }
  let slots = 0;
  for (let i = 0; i < n; i++) slots += offset[vertices[i] + 1] - offset[vertices[i]];
  const newOffset = new Int32Array(n + 1);
  const newA = new Int32Array(slots);
  const newB = new Int32Array(slots);
  const edgeIndices = new Int32Array(slots);
  let eo = 0;
  for (let i = 0; i < n; i++) {
    const v = vertices[i];
    for (let t = offset[v], tl = offset[v + 1]; t < tl; t++) {
      const j = map ? map[b[t]] : childIndexOf(vertices, b[t]);
      if (j < 0) continue;
      newA[eo] = i;
      newB[eo] = j;
      edgeIndices[eo++] = t;
    }
    newOffset[i + 1] = eo;
  }
  const picked = edgeIndices.subarray(0, eo);
  return IntAdjacencyGraph.create(
    newOffset,
    newA.slice(0, eo) as ArrayLike<number> as ArrayLike<StructureElement.UnitIndex>,
    newB.slice(0, eo) as ArrayLike<number> as ArrayLike<StructureElement.UnitIndex>,
    eo / 2,
    {
      order: arrayPickIndices(order, picked),
      flags: arrayPickIndices(flags, picked),
      key: arrayPickIndices(key, picked),
    },
    bonds.props,
  );
}

/**
 * Whether a child unit operator is its parent unit's under the common rigid motion `sameMotion` tests
 * (`SymmetryOperator.rigidMotionTest`): `Structure.transform` and its components (perf-audit-2 05/P14). It keeps every
 * distance between atoms of the units, so their inter-unit bonds are the parent's. A moved operator must also keep its
 * key, which selects index-pair bonds; `Structure.transform` keeps it.
 */
function sameRigidMotion(
  sameMotion: ReturnType<typeof SymmetryOperator.rigidMotionTest>,
  parent: SymmetryOperator,
  child: SymmetryOperator,
) {
  return sameMotion(parent, child) && (!sameMotion.moved || parent.key === child.key);
}

/**
 * Whether every unit of `child` is its `parent` unit (by id, kind and model) after one rigid motion other than the
 * identity: a moved structure (`Structure.transform`) or a component of one. The parent stays while such children are
 * made again at every move, so their inter-unit bonds are best restricted from the parent's (perf-audit-2 05/P14).
 */
export function isMovedImage(parent: Structure, child: Structure) {
  const sameMotion = SymmetryOperator.rigidMotionTest();
  for (const u of child.units) {
    const p = parent.unitMap.get(u.id);
    if (!p || p.kind !== u.kind || p.model !== u.model) return false;
    if (!sameRigidMotion(sameMotion, p.conformation.operator, u.conformation.operator)) return false;
  }
  return sameMotion.moved;
}

/**
 * The inter-unit bonds of `child` restricted from `parent`'s, or `undefined` when `child` has a unit that is not a
 * subset of `parent`'s unit with the same id and model, or the units' operators are not the parent's after one rigid
 * motion (05-model-structure P5, perf-audit-2 05/P14). The caller checks that the
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
  let wholeUnits = true;
  const sameMotion = SymmetryOperator.rigidMotionTest();
  for (const u of child.units) {
    const p = parent.unitMap.get(u.id);
    if (!p || p.kind !== u.kind || p.model !== u.model) return;
    if (!sameRigidMotion(sameMotion, p.conformation.operator, u.conformation.operator)) return;
    if (p.elements === u.elements || p.elements.length === u.elements.length) {
      if (p.elements !== u.elements && !SortedArray.areEqual(p.elements, u.elements)) return;
      maps.set(u.id, null);
      continue;
    }
    const parentIndices = getParentIndices(p.elements, u.elements);
    if (!parentIndices) return;
    maps.set(u.id, parentIndices);
    whole = false;
    wholeUnits = false;
  }
  if (whole) return bonds;
  if (wholeUnits) return wholeUnitBonds(bonds, child, maps);

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

/**
 * The bonds of a child whose units are whole parent units (a polymer or ligand component of a structure whose units
 * are chains): the parent's unit-pair lists between the child's units, shared, not rebuilt edge by edge
 * (perf-audit-2 05/P20)
 */
function wholeUnitBonds(bonds: InterUnitBonds, child: Structure, units: ReadonlyMap<number, unknown>) {
  const map = new Map<number, InterUnitGraph.UnitPairEdges<number, StructureElement.UnitIndex, InterUnitEdgeProps>[]>();
  for (const u of child.units) {
    const pairs = bonds.getConnectedUnits(u.id).filter((pair) => units.has(pair.unitB));
    if (pairs.length > 0) map.set(u.id, pairs);
  }
  return new InterUnitBonds(map);
}
