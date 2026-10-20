/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, type Unit, type StructureElement } from '../../structure.js';
import { SortedArray } from '@molstar/core/data/int';
import type { StructureSubsetBuilder } from '../../structure/util/subset-builder.js';
import type { ElementIndex } from '../../model.js';

/** The element sets of one unit id across the inputs of a union. */
interface UnionParts {
  /** The source unit's elements. */
  all: StructureElement.Set;
  sets: StructureElement.Set[];
  /** A set as large as the source unit, which is then the union. */
  full: StructureElement.Set | undefined;
  size: number;
}

/**
 * The union of parts of `source`: every input unit is a subset of the source unit with its id. Each unit's sets are
 * merged once, after all inputs are read, not input by input: merging into the accumulated set copied it per input,
 * quadratic in the unit's size for singletons (perf audit 05/P1). An input as large as the source unit is the union
 * of its unit, whatever else the other inputs hold.
 */
export function structureUnion(source: Structure, structures: Structure[]) {
  if (structures.length === 0) return Structure.Empty;
  if (structures.length === 1) return structures[0];

  const parts = new Map<number, UnionParts>();

  for (const { units } of structures) {
    for (let i = 0, _i = units.length; i < _i; i++) {
      const u = units[i];
      const p = parts.get(u.id);
      if (p) {
        if (p.full) continue;
        if (u.elements.length === p.all.length) p.full = u.elements;
        else {
          p.sets.push(u.elements);
          p.size += u.elements.length;
        }
      } else {
        const all = source.unitMap.get(u.id).elements;
        const full = u.elements.length === all.length ? u.elements : undefined;
        parts.set(u.id, { all, sets: [u.elements], full, size: u.elements.length });
      }
    }
  }

  const builder = source.subsetBuilder(true);
  parts.forEach(buildUnion, builder);
  return builder.getStructure();
}

function buildUnion(this: StructureSubsetBuilder, p: UnionParts, id: number) {
  this.setUnit(id, p.full ?? mergeSets(p));
}

/**
 * The union of a unit's sorted sets: one linear merge for two; for more, a concatenation if they follow one another
 * (residues of a unit in order), else, by their total size against the span of the source unit's elements, a sort
 * (few elements in a large unit) or marks over the span (sets that cover much of the unit, overlapping or not).
 */
function mergeSets({ all, sets, size }: UnionParts): StructureElement.Set {
  if (sets.length === 1) return sets[0];
  if (sets.length === 2) return SortedArray.union(sets[0], sets[1]);

  let inOrder = true,
    last = -1,
    largest = sets[0];
  for (const s of sets) {
    if (s.length === 0) continue;
    if (s[0] <= last) inOrder = false;
    last = s[s.length - 1];
    if (s.length > largest.length) largest = s;
  }

  if (inOrder) return SortedArray.ofSortedArray(concat(sets, size));

  const min = all[0];
  const span = all[all.length - 1] - min + 1;
  let union: StructureElement.Set;
  if (size * 16 < span) {
    union = SortedArray.deduplicate(SortedArray.ofSortedArray<ElementIndex>(concat(sets, size).sort()));
  } else {
    const marks = new Uint8Array(span);
    let count = 0;
    for (const s of sets) {
      for (let i = 0, il = s.length; i < il; i++) {
        const m = s[i] - min;
        if (marks[m]) continue;
        marks[m] = 1;
        // the whole unit: the rest adds nothing (as the fold's full-unit shortcut)
        if (++count === all.length) return all;
      }
    }
    const merged = new Int32Array(count);
    for (let i = 0, o = 0; o < count; i++) {
      if (marks[i]) merged[o++] = min + i;
    }
    union = SortedArray.ofSortedArray<ElementIndex>(merged);
  }
  // a set that holds all the others is their union: keep it rather than a copy
  return union.length === largest.length ? largest : union;
}

function concat(sets: StructureElement.Set[], size: number) {
  const ret = new Int32Array(size);
  let offset = 0;
  for (const s of sets) {
    ret.set(s, offset);
    offset += s.length;
  }
  return ret;
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
