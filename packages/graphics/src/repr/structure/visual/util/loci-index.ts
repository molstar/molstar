/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IntMap } from '@molstar/core/data/int/map';
import { StructureElement, type Unit } from '@molstar/model/model/structure';

// What every visual of a representation asks of a loci it marks or applies as a layer, computed once per loci and
// shared by the visuals (loci are immutable, the same object reaches every visual), so applying a loci costs
// O(loci + visuals) instead of O(loci × visuals) (perf-audit 07 P3).

type Element = StructureElement.Loci['elements'][number];

/** below this many elements, scanning the loci beats building and keeping an index */
const MinIndexedElements = 16;

const WholeStructure = new WeakMap<StructureElement.Loci, boolean>();
/** position in `loci.elements` by unit id, `undefined` when a unit is listed twice */
const Positions = new WeakMap<StructureElement.Loci, Map<number, number> | undefined>();

/** `StructureElement.Loci.isWholeStructure`, computed once per loci */
export function isWholeStructureLoci(loci: StructureElement.Loci) {
  if (loci.elements.length < MinIndexedElements) return StructureElement.Loci.isWholeStructure(loci);
  let whole = WholeStructure.get(loci);
  if (whole === undefined) {
    whole = StructureElement.Loci.isWholeStructure(loci);
    WholeStructure.set(loci, whole);
  }
  return whole;
}

function getPositions(loci: StructureElement.Loci) {
  if (Positions.has(loci)) return Positions.get(loci);
  const { elements } = loci;
  let positions: Map<number, number> | undefined = new Map<number, number>();
  for (let i = 0, il = elements.length; i < il; ++i) positions.set(elements[i].unit.id, i);
  if (positions.size !== elements.length) positions = undefined;
  Positions.set(loci, positions);
  return positions;
}

/**
 * The elements of `loci` on `units` (of a symmetry group, `unitIndexMap` its index by unit id), in the order of
 * `loci`: `loci` itself when all of its elements are, `undefined` when none is. Takes the group's fields rather than
 * the group: symmetry groups are object literals with an accessor each, so reading them per call is megamorphic.
 */
export function getGroupLoci(
  loci: StructureElement.Loci,
  units: ReadonlyArray<Unit>,
  unitIndexMap: IntMap<number>,
): StructureElement.Loci | undefined {
  const { elements } = loci;
  const positions =
    elements.length >= MinIndexedElements && units.length < elements.length ? getPositions(loci) : undefined;

  let picked: Element[];
  if (positions) {
    const indices: number[] = [];
    for (let i = 0, il = units.length; i < il; ++i) {
      const p = positions.get(units[i].id);
      if (p !== undefined) indices.push(p);
    }
    if (indices.length === 0) return;
    if (indices.length > 1) indices.sort((a, b) => a - b);
    picked = indices.map((p) => elements[p]);
  } else {
    // a pass without allocation first: most visuals have none of a small loci (a hover), and with one unit per group
    // a comparison of ids does
    let count = 0;
    if (units.length === 1) {
      const { id } = units[0];
      for (let i = 0, il = elements.length; i < il; ++i) if (elements[i].unit.id === id) ++count;
    } else {
      for (let i = 0, il = elements.length; i < il; ++i) if (unitIndexMap.has(elements[i].unit.id)) ++count;
    }
    if (count === 0) return;
    if (count === elements.length) return loci;
    picked = [];
    for (let i = 0, il = elements.length; i < il; ++i) {
      if (unitIndexMap.has(elements[i].unit.id)) picked.push(elements[i]);
    }
  }
  return StructureElement.Loci(loci.structure, picked);
}
