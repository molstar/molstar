/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { Expression } from '@molstar/query-language/language/expression';
import type { MVSAnnotationRow } from './schemas.js';

/** Convert an annotation row into a MolScript expression */
export function rowToExpression(row: MVSAnnotationRow): Expression {
  return StructureElement.Schema.toExpression(row);
}

/** Convert multiple annotation rows into a MolScript expression.
 * (with union semantics, i.e. an atom qualifies if it qualifies for at least one of the rows) */
export function rowsToExpression(rows: readonly MVSAnnotationRow[]): Expression {
  return StructureElement.Schema.toExpression({
    items: rows as StructureElement.SchemaItem[],
  });
}
