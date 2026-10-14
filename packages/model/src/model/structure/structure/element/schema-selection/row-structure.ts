/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import type { Structure } from '../../structure.js';
import { Unit } from '../../unit.js';
import { ElementRanges } from './element-ranges.js';
import { IndicesAndSortings } from './indexing.js';
import type { SchemaItem } from '../schema.js';
import {
  getAtomRangesForRows,
  getGaussianRangesForRows,
  getSphereRangesForRows,
  rowsUseOperators,
} from './selections.js';

/** Helper for caching element ranges qualifying to a group of annotation rows, per `Unit`
 * (units of the same model and kind share the ranges, unless some row refers to a specific `instance_id` or `operator_name`). */
export class ElementRangesCache {
  private readonly cache: { [key: string]: ElementRanges } = {};
  private readonly hasOperators: boolean;

  constructor(private readonly rows: readonly SchemaItem[]) {
    this.hasOperators = rowsUseOperators(rows);
  }

  get(unit: Unit): ElementRanges {
    const { instanceId, name } = unit.conformation.operator;
    const key = `${unit.model.id}:${unit.kind}:${this.hasOperators ? `${instanceId}|${name}` : '*'}`;
    return (this.cache[key] ??= this.compute(unit));
  }
  private compute(unit: Unit): ElementRanges {
    const operator = unit.conformation.operator;
    const indices = IndicesAndSortings.get(unit.model);
    switch (unit.kind) {
      case Unit.Kind.Atomic:
        return getAtomRangesForRows(this.rows, unit.model, operator, indices);
      case Unit.Kind.Spheres:
        return getSphereRangesForRows(this.rows, unit.model, operator, indices);
      case Unit.Kind.Gaussians:
        return getGaussianRangesForRows(this.rows, unit.model, operator, indices);
    }
  }
}

/** Return the substructure of `structure` consisting of all elements which satisfy at least one of `rows`,
 * found through the indices of each model, so the cost does not grow with the number of atoms times the number of rows.
 *
 * This is what the generator `element-schema` (`Schema.toExpression({ items: rows })`) runs. Rows select what their criteria say: a residue
 * without `label_seq_id` never satisfies a `label_seq_id` criterion (not even 0), criteria which only exist for atoms select nothing on coarse
 * elements (with a warning), a range with `beg > end` selects nothing, and a range with insertion codes is a range of residues in file order.
 * When every element is selected, it returns `structure` itself. */
export function substructureFromRows(structure: Structure, rows: readonly SchemaItem[]): Structure {
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
