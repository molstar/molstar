/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { type Structure, Unit } from '@molstar/model/model/structure';
import { ElementRanges } from './element-ranges.js';
import { IndicesAndSortings } from './indexing.js';
import type { MVSAnnotationRow } from './schemas.js';
import { getAtomRangesForRows, getGaussianRangesForRows, getSphereRangesForRows } from './selections.js';
import { isDefined } from './utils.js';

/** Helper for caching element ranges qualifying to a group of annotation rows, per `Unit`
 * (units of the same model and kind share the ranges, unless some row refers to a specific `instance_id`). */
export class ElementRangesCache {
  private readonly cache: { [key: string]: ElementRanges } = {};
  private readonly hasOperators: boolean;

  constructor(private readonly rows: readonly MVSAnnotationRow[]) {
    this.hasOperators = rows.some((row) => isDefined(row.instance_id));
  }

  get(unit: Unit): ElementRanges {
    const instanceId = unit.conformation.operator.instanceId;
    const key = `${unit.model.id}:${unit.kind}:${this.hasOperators ? instanceId : '*'}`;
    return (this.cache[key] ??= this.compute(unit));
  }
  private compute(unit: Unit): ElementRanges {
    const instanceId = unit.conformation.operator.instanceId;
    const indices = IndicesAndSortings.get(unit.model);
    switch (unit.kind) {
      case Unit.Kind.Atomic:
        return getAtomRangesForRows(this.rows, unit.model, instanceId, indices);
      case Unit.Kind.Spheres:
        return getSphereRangesForRows(this.rows, unit.model, instanceId, indices);
      case Unit.Kind.Gaussians:
        return getGaussianRangesForRows(this.rows, unit.model, instanceId, indices);
    }
  }
}

/** Return the substructure of `structure` consisting of all elements which satisfy at least one of `rows`
 * (the same result as running `rowsToExpression(rows)` on `structure`, but found through the indices of each model,
 * so the cost does not grow with the number of atoms times the number of rows). */
export function substructureFromRows(structure: Structure, rows: readonly MVSAnnotationRow[]): Structure {
  if (rows.length === 0 || structure.isEmpty) return structure.subsetBuilder(true).getStructure(); // (empty)
  const cache = new ElementRangesCache(rows);
  const builder = structure.subsetBuilder(true);
  for (const unit of structure.units) {
    // (a new array per unit: the builder keeps it)
    const elements = ElementRanges.selectElementsInRanges(unit.elements, cache.get(unit));
    if (elements.length > 0) builder.setUnit(unit.id, elements);
  }
  if (builder.elementCount === structure.elementCount) return structure; // everything selected: keep the identity (and the root-only shortcuts and caches)
  return builder.getStructure();
}
