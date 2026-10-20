/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, type Unit, type StructureElement } from '../../structure.js';
import { SortedArray } from '@molstar/core/data/int';
import type { StructureSubsetBuilder } from '../../structure/util/subset-builder.js';
import type { ElementIndex } from '../../model.js';

/** The element sets of one unit id across the inputs of a union. */
interface UnionParts {
  sets: StructureElement.Set[];
  /** A set as large as the source unit, which is then the union. */
  full: StructureElement.Set | undefined;
  size: number;
}

/**
 * The union of parts of `source`. Each unit's sets are merged once, after all inputs are read, not input by input:
 * merging into the accumulated set copied it per input, quadratic in the unit's size for singletons
 * (perf audit 05/P1).
 */
export function structureUnion(source: Structure, structures: Structure[]) {
  if (structures.length === 0) return Structure.Empty;
  if (structures.length === 1) return structures[0];

  const parts = new Map<number, UnionParts>();

  for (const { units } of structures) {
    for (let i = 0, _i = units.length; i < _i; i++) {
      const u = units[i];
      const isFull = u.elements.length === source.unitMap.get(u.id).elements.length;
      const p = parts.get(u.id);
      if (p) {
        if (p.full) continue;
        if (isFull) p.full = u.elements;
        else {
          p.sets.push(u.elements);
          p.size += u.elements.length;
        }
      } else {
        parts.set(u.id, { sets: [u.elements], full: isFull ? u.elements : undefined, size: u.elements.length });
      }
    }
  }

  const builder = source.subsetBuilder(true);
  parts.forEach(buildUnion, builder);
  return builder.getStructure();
}

function buildUnion(this: StructureSubsetBuilder, p: UnionParts, id: number) {
  this.setUnit(id, p.full ?? mergeSets(p.sets, p.size));
}

/** The union of sorted sets: one linear merge for two, else a concatenation, sorted and deduplicated if needed. */
function mergeSets(sets: StructureElement.Set[], size: number): StructureElement.Set {
  if (sets.length === 1) return sets[0];
  if (sets.length === 2) return SortedArray.union(sets[0], sets[1]);

  const merged = new Int32Array(size);
  let offset = 0,
    inOrder = true,
    largest = sets[0];
  for (const s of sets) {
    // sets that follow one another (residues of a unit in order) concatenate to a sorted, duplicate-free set
    if (offset > 0 && merged[offset - 1] >= s[0]) inOrder = false;
    merged.set(s, offset);
    offset += s.length;
    if (s.length > largest.length) largest = s;
  }
  if (inOrder) return SortedArray.ofSortedArray(merged);

  merged.sort();
  const union = SortedArray.deduplicate(SortedArray.ofSortedArray<ElementIndex>(merged));
  // a set that holds all the others is their union: keep it rather than a copy
  return union.length === largest.length ? largest : union;
}

export function structureAreEqual(sA: Structure, sB: Structure): boolean {
  if (sA === sB) return true;

  if (sA.units.length !== sB.units.length) return false;

  const aU = sA.units,
    bU = sB.unitMap;
  for (let i = 0, _i = aU.length; i < _i; i++) {
    const u = aU[i];
    if (!bU.has(u.id)) return false;
    const v = bU.get(u.id);
    if (!SortedArray.areEqual(u.elements, v.elements)) return false;
  }

  return true;
}

export function structureAreIntersecting(sA: Structure, sB: Structure): boolean {
  if (sA === sB) return true;

  let a, b;
  if (sA.units.length < sB.units.length) {
    a = sA;
    b = sB;
  } else {
    a = sB;
    b = sA;
  }

  const aU = a.units,
    bU = b.unitMap;

  for (let i = 0, _i = aU.length; i < _i; i++) {
    const u = aU[i];
    if (!bU.has(u.id)) continue;
    const v = bU.get(u.id);
    if (SortedArray.areIntersecting(u.elements, v.elements)) return true;
  }

  return false;
}

export function structureIntersect(sA: Structure, sB: Structure): Structure {
  if (sA === sB) return sA;
  if (!structureAreIntersecting(sA, sB)) return Structure.Empty;

  let a, b;
  if (sA.units.length < sB.units.length) {
    a = sA;
    b = sB;
  } else {
    a = sB;
    b = sA;
  }

  const aU = a.units,
    bU = b.unitMap;
  const units: Unit[] = [];

  for (let i = 0, _i = aU.length; i < _i; i++) {
    const u = aU[i];
    if (!bU.has(u.id)) continue;
    const v = bU.get(u.id);
    if (SortedArray.areIntersecting(u.elements, v.elements)) {
      const int = SortedArray.intersect(u.elements, v.elements);
      units[units.length] = u.getChild(int);
    }
  }

  return Structure.create(units, { parent: sA.parent || sB.parent });
}

export function structureSubtract(a: Structure, b: Structure): Structure {
  if (a === b) return Structure.Empty;
  if (!structureAreIntersecting(a, b)) return a;

  const aU = a.units,
    bU = b.unitMap;
  const units: Unit[] = [];

  for (let i = 0, _i = aU.length; i < _i; i++) {
    const u = aU[i];
    if (!bU.has(u.id)) {
      units[units.length] = u;
      continue;
    }
    const v = bU.get(u.id);
    const sub = SortedArray.subtract(u.elements, v.elements);
    if (sub.length > 0) {
      units[units.length] = u.getChild(sub);
    }
  }

  return Structure.create(units, { parent: a.parent || b.parent });
}
