/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SchemaItem } from '../schema.js';

/** Every field of a `SchemaItem` (the compiler checks that none is missing) */
const Fields: { [K in keyof Required<SchemaItem>]: true } = {
  operator_name: true,
  instance_id: true,
  label_entity_id: true,
  label_asym_id: true,
  auth_asym_id: true,
  label_seq_id: true,
  auth_seq_id: true,
  residue_index: true,
  label_comp_id: true,
  auth_comp_id: true,
  pdbx_PDB_ins_code: true,
  beg_label_seq_id: true,
  end_label_seq_id: true,
  beg_auth_seq_id: true,
  end_auth_seq_id: true,
  beg_auth_ins_code: true,
  end_auth_ins_code: true,
  label_atom_id: true,
  auth_atom_id: true,
  type_symbol: true,
  atom_id: true,
  atom_index: true,
};

export const SchemaItemFields = Object.keys(Fields) as readonly (keyof SchemaItem)[];

export function isSchemaItemField(name: string): name is keyof SchemaItem {
  return Object.hasOwn(Fields, name);
}

/** A value which an item can have: a missing value (`undefined`, `null`, `NaN`) is not a criterion. */
export function isFieldValue(value: unknown): value is string | number | boolean {
  return value !== undefined && value !== null && !(typeof value === 'number' && Number.isNaN(value));
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

/** The items of a table given by its fields: a field which is an array is a column of the table (a missing value is a hole in it),
 * any other field has the same value in every row. Without any column, there is one row. All the columns have the same length. */
export function schemaItemsFromFields(fields: { readonly [name: string]: unknown }): SchemaItem[] {
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
  return items;
}
