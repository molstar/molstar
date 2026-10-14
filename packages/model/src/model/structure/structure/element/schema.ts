/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { Expression } from '@molstar/query-language/language/expression';
import type { Structure } from '../structure.js';
import { Bundle } from './bundle.js';
import { Loci } from './loci.js';
import { SchemaItemFields, isFieldValue, itemWithValues } from './schema-selection/fields.js';
import { substructureFromRows } from './schema-selection/row-structure.js';

export interface SchemaItem {
  /** Corresponds to SymmetryOperator.name, e.g. 1_555, ASM_5 */
  operator_name?: string;
  /** Corresponds to SymmetryOperator.instanceId, e.g. 1_555, ASM-X0-5 */
  instance_id?: string;
  label_entity_id?: string;
  label_asym_id?: string;
  auth_asym_id?: string;
  label_seq_id?: number;
  auth_seq_id?: number;
  residue_index?: number;
  label_comp_id?: string;
  auth_comp_id?: string;
  pdbx_PDB_ins_code?: string;
  beg_label_seq_id?: number;
  end_label_seq_id?: number;
  beg_auth_seq_id?: number;
  end_auth_seq_id?: number;
  /** Insertion code of the first residue of the range which starts at `beg_auth_seq_id` (`''`: the residue without an insertion code).
   *
   * If `beg_auth_ins_code` or `end_auth_ins_code` is given, the range is not a range of `auth_seq_id` values but of residues
   * in the order of the input file: from the first residue with `beg_auth_seq_id` (and this insertion code) to the last residue with
   * `end_auth_seq_id` (and `end_auth_ins_code`), in each chain selected by the other criteria. A bound without an insertion code is the
   * whole `auth_seq_id` (the first residue of it for `beg`, the last for `end`), a bound which is not in the chain selects nothing from it,
   * and so does a range which ends before it begins. An insertion code without its `auth_seq_id` selects nothing. */
  beg_auth_ins_code?: string;
  /** Insertion code of the last residue of the range which ends at `end_auth_seq_id`, see `beg_auth_ins_code`. */
  end_auth_ins_code?: string;
  label_atom_id?: string;
  auth_atom_id?: string;
  type_symbol?: string;
  atom_id?: number;
  atom_index?: number;
}

export interface SchemaItems {
  // The prefix is applied to each item in the list.
  // This is useful for example for referencing multiple atoms in a single
  // residue or multiple residues in a chain, e.g.:
  //    `{ prefix: { label_asym_id: 'A' }, items: { label_seq_id: [1, 3, 5] } }`
  prefix?: SchemaItem;

  // Depending on the use-case, items can be specified either
  // as a list of objects or using a more memory-efficient object of lists.
  items: SchemaItem[] | { [K in keyof SchemaItem]: SchemaItem[K][] };
}

export type Schema = SchemaItem | SchemaItems;

function isItems(schema: Schema): schema is SchemaItems {
  return !!(schema as any).items;
}

/** The expression of one row (the fields which it has) */
function itemExpression(item: SchemaItem): Expression {
  const args: Record<string, Expression> = {};
  for (const name of SchemaItemFields) {
    const value = item[name];
    if (isFieldValue(value)) args[name] = value;
  }
  return MS.struct.type.schemaItem(args);
}

/** From this many rows on, the rows are given as columns: a list of values is passed through by the compiler as it is (so the cost of
 * the expression does not grow with the number of rows) and does not repeat the names of the fields. */
const MinRowsForColumns = 10;

/** The column of a field: a missing value is `null` (what JSON makes of `undefined`) */
function columnOf(values: ArrayLike<unknown>): Expression {
  const column: Expression.ListValue[] = new Array(values.length);
  for (let i = 0; i < values.length; i++) {
    const value = values[i];
    column[i] = isFieldValue(value) ? value : null;
  }
  return column;
}

/** The expression of a table given by its columns: the prefix and the lists are the fields of one item (the prefix is the same in every row). */
function columnsExpression(
  prefix: SchemaItem | undefined,
  columns: { readonly [name: string]: ArrayLike<unknown> | undefined },
): Expression {
  const args: Record<string, Expression> = {};
  for (const name of SchemaItemFields) {
    const value = prefix?.[name];
    if (isFieldValue(value)) args[name] = value;
  }
  for (const name of SchemaItemFields) {
    const column = columns[name];
    if (column) args[name] = columnOf(column);
  }
  return MS.struct.type.schemaItem(args);
}

/** The expression of a schema is the generator `element-schema`: it finds the elements through the indices of the models (instead of testing
 * every element for every item), as a single set of the elements which satisfy at least one item. A table of many rows is given by columns. */
function toExpression(schema: Schema): Expression {
  if (isItems(schema) && !Array.isArray(schema.items)) {
    return MS.struct.generator.elementSchema([columnsExpression(schema.prefix, schema.items)]);
  }
  const items: SchemaItem[] = [];
  forEachItem(schema, (item) => items.push(itemWithValues(item)));
  if (items.length < MinRowsForColumns) return MS.struct.generator.elementSchema(items.map(itemExpression));
  const columns: { [name: string]: unknown[] } = {};
  for (const name of SchemaItemFields) {
    if (items.some((item) => isFieldValue(item[name]))) columns[name] = items.map((item) => item[name]);
  }
  return MS.struct.generator.elementSchema([columnsExpression(undefined, columns)]);
}

/**
 * Iterate over all items in a structure element schema.
 * @param schema Schema to iterate over
 * @param f Function called for each item in the schema.
 *          The value passed to the function can be mutable and should not be
 *          modified => make a copy if the value is used outside the callback.
 */
function forEachItem(schema: Schema, f: (item: SchemaItem) => void) {
  if (isItems(schema)) {
    if (Array.isArray(schema.items)) {
      if (schema.prefix) {
        for (const item of schema.items) {
          f(Object.assign({}, schema.prefix, item));
        }
      } else {
        for (const item of schema.items) {
          f(item);
        }
      }
    } else {
      const current: any = { ...schema.prefix };
      const keys = Object.keys(schema.items);
      const n = (schema.items as any)[keys[0]].length;

      for (let i = 0; i < n; i++) {
        for (const k of keys) {
          current[k] = (schema.items as any)[k][i];
        }
        f(current);
      }
    }
  } else {
    f(schema);
  }
}

function toLoci(structure: Structure, schema: Schema): Loci {
  return Loci.fromExpression(structure, toExpression(schema));
}

function toBundle(structure: Structure, schema: Schema): Bundle {
  return Bundle.fromLoci(toLoci(structure, schema));
}

/** The substructure of `structure` of the elements which satisfy at least one item of `schema`
 * (`structure` itself if every element satisfies one). */
function toStructure(structure: Structure, schema: Schema): Structure {
  const items: SchemaItem[] = [];
  // (copies: the items of an object of arrays are produced by reusing one object; a missing value is not a criterion)
  forEachItem(schema, (item) => items.push(itemWithValues(item)));
  return substructureFromRows(structure, items);
}

export const Schema = {
  forEachItem,
  toExpression,
  toLoci,
  toBundle,
  toStructure,
};
