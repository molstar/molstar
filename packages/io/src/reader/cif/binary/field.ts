/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import type * as Data from '../data-model.js';
import { type EncodedColumn, decode } from '@molstar/io/common/binary-cif';
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '@molstar/core/util/number-parser';

/**
 * A BinaryCIF field. With `rowCount` (the category's), the column decodes on the first access to its values, so the
 * columns a reader never touches are never decoded (perf-audit-2 04/P5); without it, it decodes here.
 */
export function Field(column: EncodedColumn, rowCount?: number): Data.CifField {
  let impl: Data.CifField | undefined;
  const get = () => impl ?? (impl = decodedField(column));
  const field = rowCount === void 0 ? get() : void 0;

  return {
    get __array() {
      return get().__array;
    },
    binaryEncoding: column.data.encoding,
    isDefined: true,
    rowCount: field ? field.rowCount : rowCount!,
    str: (row) => get().str(row),
    int: (row) => get().int(row),
    float: (row) => get().float(row),
    valueKind: (row) => get().valueKind(row),
    areValuesEqual: (rowA, rowB) => get().areValuesEqual(rowA, rowB),
    toStringArray: (params) => get().toStringArray(params),
    toIntArray: (params) => get().toIntArray(params),
    toFloatArray: (params) => get().toFloatArray(params),
  };
}

function decodedField(column: EncodedColumn): Data.CifField {
  const mask = column.mask ? (decode(column.mask) as number[]) : void 0;
  const data = decode(column.data);
  const isNumeric = ColumnHelpers.isTypedArray(data);
  /** Masked values are missing data, hence `int`/`float` return the default value (0) for them,
   * consistent with the text and mmCIF parsers. If the mask is trivial, i.e. nothing is
   * actually masked, we can skip the per-value check altogether. */
  const hasMissingValues = !!mask && mask.some((v) => v !== Column.ValueKinds.Present);
  /** only accessed when `hasMissingValues` is true, which implies `mask` is defined */
  const presentMask = mask as number[];

  const str: Data.CifField['str'] = isNumeric
    ? mask
      ? (row) => (mask[row] === Column.ValueKinds.Present ? '' + data[row] : '')
      : (row) => '' + data[row]
    : mask
      ? (row) => (mask[row] === Column.ValueKinds.Present ? data[row] : '')
      : (row) => data[row];

  const int: Data.CifField['int'] = isNumeric
    ? hasMissingValues
      ? (row) => (presentMask[row] === Column.ValueKinds.Present ? data[row] : 0)
      : (row) => data[row]
    : (row) => {
        const v = data[row];
        return fastParseInt(v, 0, v.length);
      };

  const float: Data.CifField['float'] = isNumeric
    ? hasMissingValues
      ? (row) => (presentMask[row] === Column.ValueKinds.Present ? data[row] : 0)
      : (row) => data[row]
    : (row) => {
        const v = data[row];
        return fastParseFloat(v, 0, v.length);
      };

  const valueKind: Data.CifField['valueKind'] = mask
    ? (row) => mask[row] as Column.ValueKind
    : (row) => Column.ValueKinds.Present;

  const rowCount = data.length;

  return {
    __array: data,
    binaryEncoding: column.data.encoding,
    isDefined: true,
    rowCount,
    str,
    int,
    float,
    valueKind,
    areValuesEqual: (rowA, rowB) => data[rowA] === data[rowB],
    toStringArray: (params) => ColumnHelpers.createAndFillArray(rowCount, str, params),
    toIntArray: isNumeric
      ? hasMissingValues
        ? (params) => ColumnHelpers.typedArrayWindowMasked(data, presentMask, 0, params)
        : (params) => ColumnHelpers.typedArrayWindow(data, params)
      : (params) => ColumnHelpers.createAndFillArray(rowCount, int, params),
    toFloatArray: isNumeric
      ? hasMissingValues
        ? (params) => ColumnHelpers.typedArrayWindowMasked(data, presentMask, 0, params)
        : (params) => ColumnHelpers.typedArrayWindow(data, params)
      : (params) => ColumnHelpers.createAndFillArray(rowCount, float, params),
  };
}
