/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SchemaItem } from '../schema.js';

/** Every field of a `SchemaItem` and the type of its values (the compiler checks that none is missing) */
const FieldTypes: { [K in keyof Required<SchemaItem>]: 'string' | 'number' } = {
  operator_name: 'string',
  instance_id: 'string',
  label_entity_id: 'string',
  label_asym_id: 'string',
  auth_asym_id: 'string',
  label_seq_id: 'number',
  auth_seq_id: 'number',
  residue_index: 'number',
  label_comp_id: 'string',
  auth_comp_id: 'string',
  pdbx_PDB_ins_code: 'string',
  beg_label_seq_id: 'number',
  end_label_seq_id: 'number',
  beg_auth_seq_id: 'number',
  end_auth_seq_id: 'number',
  beg_auth_ins_code: 'string',
  end_auth_ins_code: 'string',
  label_atom_id: 'string',
  auth_atom_id: 'string',
  type_symbol: 'string',
  atom_id: 'number',
  atom_index: 'number',
};

export const SchemaItemFields = Object.keys(FieldTypes) as readonly (keyof SchemaItem)[];

/** Is `value` there? `undefined`, `null` and `NaN` (what a number which is missing is in an array of numbers) are not: they are not a criterion. */
export function hasValue<T>(value: T | undefined | null): value is T {
  return value !== undefined && value !== null && !Number.isNaN(value as number);
}

/** Does `row` have a value of another type than its field (a number given as text, an identifier given as a number, ...)? Values are matched
 * as they are, so such a row cannot select anything. */
export function hasWrongTypes(row: SchemaItem): boolean {
  for (const name of SchemaItemFields) {
    const value = row[name];
    if (hasValue(value) && typeof value !== FieldTypes[name]) return true;
  }
  return false;
}

export function isSchemaItemField(name: string): name is keyof SchemaItem {
  return Object.hasOwn(FieldTypes, name);
}

/** A value which an item can have: a missing value (`undefined`, `null`, `NaN`) is not a criterion. */
export function isFieldValue(value: unknown): value is string | number | boolean {
  return hasValue(value);
}

/** A copy of the item with the values it has (without the fields which are missing) */
export function itemWithValues(item: SchemaItem): SchemaItem {
  const result: Record<string, unknown> = {};
  for (const name of SchemaItemFields) {
    const value = item[name];
    if (isFieldValue(value)) result[name] = value;
  }
  return result as SchemaItem;
}

const RowsOfItems = new WeakSet<object>();
/** Are these the rows which `schemaItemsFromFields` made (an argument of `element-schema` has to be)? */
export function isRowsOfSchemaItem(value: unknown): value is SchemaItem[] {
  return typeof value === 'object' && value !== null && RowsOfItems.has(value);
}

/** The items of a table given by its fields: a field which is an array is a column of the table (a missing value is a hole in it),
 * any other field has the same value in every row. Without any column, there is one row. All the columns have the same length. */
export function schemaItemsFromFields(fields: { readonly [name: string]: unknown } = {}): SchemaItem[] {
  let nRows = -1;
  for (const name of Object.keys(fields)) {
    if (!isSchemaItemField(name)) throw new Error(`ArgumentError: '${name}' is not a field of a schema item.`);
    const value = fields[name];
    if (!Array.isArray(value)) continue;
    if (nRows >= 0 && value.length !== nRows) {
      throw new Error(
        `ArgumentError: the columns of a schema item have different lengths ('${name}' has ${value.length}, others ${nRows}).`,
      );
    }
    nRows = value.length;
  }

  const names = Object.keys(fields);
  const items: SchemaItem[] = [];
  for (let i = 0, n = nRows < 0 ? 1 : nRows; i < n; i++) {
    const item: Record<string, unknown> = {};
    for (const name of names) {
      const value = fields[name];
      const v = Array.isArray(value) ? value[i] : value;
      if (isFieldValue(v)) item[name] = v;
    }
    items.push(item as SchemaItem);
  }
  RowsOfItems.add(items);
  return items;
}
