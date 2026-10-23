/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRangeArray, makeBuckets } from '@molstar/core/data/util';
import { Column, ColumnHelpers } from '@molstar/core/data/db';
import type { RuntimeContext } from '@molstar/core/task';
import type { AtomSite } from './schema.js';

export type SortedAtomSite = {
  atom_site: AtomSite;
  sourceIndex: Column<number>;
};

export async function sortAtomSite(
  ctx: RuntimeContext,
  atom_site: AtomSite,
  start: number,
  end: number,
): Promise<SortedAtomSite> {
  const indices = createRangeArray(start, end - 1);

  const { label_entity_id, label_asym_id, label_seq_id } = atom_site;
  const entityBuckets = makeBuckets(indices, label_entity_id.value);
  if (ctx.shouldUpdate) await ctx.update();
  for (let ei = 0, _eI = entityBuckets.length - 1; ei < _eI; ei++) {
    const chainBuckets = makeBuckets(indices, label_asym_id.value, {
      start: entityBuckets[ei],
      end: entityBuckets[ei + 1],
    });
    for (let cI = 0, _cI = chainBuckets.length - 1; cI < _cI; cI++) {
      const aI = chainBuckets[cI];
      // are we in HETATM territory?
      if (label_seq_id.valueKind(aI) !== Column.ValueKinds.Present) continue;

      makeBuckets(indices, label_seq_id.value, { sort: true, start: aI, end: chainBuckets[cI + 1] });
      if (ctx.shouldUpdate) await ctx.update();
    }
    if (ctx.shouldUpdate) await ctx.update();
  }

  // Already in order (every model of most files): typed-array columns are read in place, where `Table.view` copied
  // them (perf-audit 04/P5); the others are views as before (a plain array keeps the array fast paths of the
  // hierarchy comparison). The whole table when the range is all of it.
  if (isRange(indices, start)) {
    if (start === 0 && end === atom_site._rowCount) return { atom_site, sourceIndex: Column.ofIntArray(indices) };
    return {
      atom_site: lazyTable(atom_site, indices.length, (c) =>
        floatRows(
          c,
          indices.length,
          () => {
            // a typed array read as its values (`Column.view`'s own test for reading the array)
            const inPlace =
              !!c.__array && ColumnHelpers.isTypedArray(c.__array) && typeof c.value(0) === typeof c.__array[0];
            return inPlace ? Column.window(c, start, end) : Column.view(c, indices);
          },
          (params) => {
            const { start: s, end: e } = ColumnHelpers.getArrayBounds(indices.length, params);
            const rows = c.toArray({ array: Float32Array, start: start + s, end: start + e });
            return rows instanceof Float32Array ? rows : Float32Array.from(rows);
          },
        ),
      ),
      sourceIndex: Column.ofIntArray(indices),
    };
  }

  return {
    atom_site: lazyTable(atom_site, indices.length, (c) =>
      floatRows(
        c,
        indices.length,
        () => Column.view(c, indices),
        (params) => {
          const source = c.toArray({ array: Float32Array });
          const { start, end } = ColumnHelpers.getArrayBounds(indices.length, params);
          const out = new Float32Array(end - start);
          for (let i = start; i < end; i++) out[i - start] = source[indices[i]];
          return out;
        },
      ),
    ),
    sourceIndex: Column.ofIntArray(indices),
  };
}

/**
 * The rows of `atom_site` as a table whose columns are made by `column` on their first read (perf-audit-2 04/P5):
 * the columns a model never reads are not copied, and a BinaryCIF source column is not decoded (`binary/field.ts`).
 */
function lazyTable(atom_site: AtomSite, rowCount: number, column: (c: Column<any>) => Column<any>): AtomSite {
  const ret = Object.create(null);
  ret._rowCount = rowCount;
  ret._columns = Object.keys(atom_site._schema);
  ret._schema = atom_site._schema;
  for (const k of ret._columns) {
    let c: Column<any> | undefined;
    Object.defineProperty(ret, k, {
      get: () => c ?? (c = column((atom_site as any)[k])),
      set: (v: Column<any>) => {
        c = v;
      },
      enumerable: true,
      configurable: true,
    });
  }
  return ret;
}

/**
 * A float column's rows: `toArray({ array: Float32Array })` reads them through `float32` from the source's own Float32
 * array (a BinaryCIF column decodes straight to Float32: the conformation's coordinates, perf-audit-2 04/P5 part 2),
 * everything else through the `rows` column, made on first use. Other columns are `rows()` as they are.
 */
function floatRows(
  c: Column<any>,
  rowCount: number,
  rows: () => Column<any>,
  float32: (params?: Column.ToArrayParams<number>) => Float32Array,
): Column<any> {
  if (c.schema.valueType !== 'float' || !c.isDefined || c.rowCount === 0) return rows();
  let column: Column<any> | undefined;
  const get = () => column ?? (column = rows());
  return {
    schema: c.schema,
    get __array() {
      return get().__array;
    },
    isDefined: true,
    rowCount,
    value: (row) => get().value(row),
    valueKind: (row) => get().valueKind(row),
    areValuesEqual: (rowA, rowB) => get().areValuesEqual(rowA, rowB),
    toArray: (params) => (params?.array === Float32Array ? float32(params) : get().toArray(params)),
  };
}

/** Whether `indices[i] === start + i` for every i. */
function isRange(indices: ArrayLike<number>, start: number) {
  for (let i = 0, il = indices.length; i < il; ++i) {
    if (indices[i] !== start + i) return false;
  }
  return true;
}
