/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import type { ElementIndex } from '../../../model/indexing.js';
import type { Structure } from '../../structure.js';
import type { UnitIndex } from '../element.js';
import { Loci } from '../loci.js';
import { type Schema, type SchemaItem, Schema as SchemaFns } from '../schema.js';
import type { ElementRanges } from './element-ranges.js';
import { ElementRangesCache } from './row-structure.js';

/** The positions in `elements` (sorted) of the elements which lie in `ranges` (sorted and disjoint), as an `OrderedSet`;
 * `undefined` when there are none. */
export function unitIndicesInRanges(
  elements: SortedArray<ElementIndex>,
  ranges: ElementRanges,
): OrderedSet<UnitIndex> | undefined {
  const nElements = elements.length;
  const nRanges = ranges.from.length;
  if (nElements === 0 || nRanges === 0) return undefined;

  // The elements in a range are `elements[i0..i1)`. Look the ends of each range up in the elements when there are
  // fewer ranges than elements, otherwise walk the elements once.
  const begins: number[] = [];
  const ends: number[] = [];
  if (nRanges <= nElements) {
    for (let r = 0; r < nRanges; r++) {
      const i0 = SortedArray.findPredecessorIndex(elements, ranges.from[r]);
      const i1 = SortedArray.findPredecessorIndex(elements, ranges.to[r]);
      if (i1 > i0) {
        if (ends.length > 0 && ends[ends.length - 1] === i0) ends[ends.length - 1] = i1;
        else {
          begins.push(i0);
          ends.push(i1);
        }
      }
    }
  } else {
    let r = 0;
    let open = -1;
    for (let i = 0; i < nElements; i++) {
      const e = elements[i];
      while (r < nRanges && ranges.to[r] <= e) r++;
      if (r === nRanges) break;
      if (ranges.from[r] <= e) {
        if (open < 0) open = i;
      } else if (open >= 0) {
        begins.push(open);
        ends.push(i);
        open = -1;
      }
    }
    if (open >= 0) {
      begins.push(open);
      ends.push(nElements);
    }
  }

  if (begins.length === 0) return undefined;
  if (begins.length === 1) return OrderedSet.ofBounds(begins[0] as UnitIndex, ends[0] as UnitIndex);
  let n = 0;
  for (let k = 0; k < begins.length; k++) n += ends[k] - begins[k];
  const indices = new Int32Array(n);
  let o = 0;
  for (let k = 0; k < begins.length; k++) {
    for (let i = begins[k]; i < ends[k]; i++) indices[o++] = i;
  }
  return OrderedSet.ofSortedArray(indices as unknown as SortedArray<UnitIndex>);
}

/** The items of `schema` (copies: the items of an object of arrays are produced by reusing one object) */
function itemsOf(schema: Schema): SchemaItem[] {
  const items: SchemaItem[] = [];
  SchemaFns.forEachItem(schema, (item) => items.push({ ...item }));
  return items;
}

/** The `Loci` of the elements of `structure` which satisfy at least one item of `schema`, resolved through the indices of the models
 * (see `Schema.toLoci` and its `useIndex` option). */
export function schemaToLociViaIndex(structure: Structure, schema: Schema): Loci {
  const items = itemsOf(schema);
  if (items.length === 0 || structure.isEmpty) return Loci.none(structure);
  const cache = new ElementRangesCache(items);
  const elements: Loci['elements'][number][] = [];
  for (const unit of structure.units) {
    const indices = unitIndicesInRanges(unit.elements, cache.get(unit));
    if (indices) elements.push({ unit, indices });
  }
  // (`structure.units` is sorted by unit id, like the elements of the loci of a query)
  return Loci(structure, elements);
}
