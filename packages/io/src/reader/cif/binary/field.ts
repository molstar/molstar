/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import type * as Data from '../data-model.js';
import { type EncodedColumn, Encoding, decode, decodeStringDictionary } from '@molstar/io/common/binary-cif';
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '@molstar/core/util/number-parser';

/**
 * A BinaryCIF field. With `rowCount` (the category's), the column decodes on the first access to its values, so the
 * columns a reader never touches are never decoded (perf-audit-2 04/P5); without it, it decodes here.
 */
export function Field(column: EncodedColumn, rowCount?: number): Data.CifField {
  let impl: Data.CifField | undefined;
  const get = () => impl ?? (impl = decodedField(column));
  const field = rowCount === void 0 ? get() : void 0;

  // A quantized float column read as Float32 (the conformation's coordinates) decodes straight to Float32: each value
  // computed in double as by the source type and rounded once, so equal to a Float32 copy of the source-type array.
  // The other accessors decode to the source type, so writers see exact values; both kinds of access decode twice.
  const quantized =
    column.data.encoding[0]?.kind === 'FixedPoint' || column.data.encoding[0]?.kind === 'IntervalQuantization';
  let float32: Data.CifField | undefined;
  const toFloatArray: Data.CifField['toFloatArray'] = quantized
    ? (params) => {
        if (params?.array !== Float32Array) return get().toFloatArray(params);
        if (!float32) {
          const [scale, ...rest] = column.data.encoding;
          const data = { ...column.data, encoding: [{ ...scale, srcType: Encoding.FloatDataType.Float32 }, ...rest] };
          float32 = decodedField({ ...column, data } as EncodedColumn);
        }
        return float32.toFloatArray(params);
      }
    : (params) => get().toFloatArray(params);

  return {
    storage: () => get().storage!(),
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
    toFloatArray,
  };
}

function decodedField(column: EncodedColumn): Data.CifField {
  const mask = column.mask ? (decode(column.mask) as number[]) : void 0;
  const dictionary = column.data.encoding[0]?.kind === 'StringArray' ? decodeStringDictionary(column.data) : void 0;
  if (dictionary) {
    const field = dictionaryField(column, dictionary.indices, dictionary.values, mask);
    if (field) return field;
  }
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
  const storage: Column.Storage = { kind: 'array', array: data };

  return {
    storage: () => storage,
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

/**
 * A string column kept as its dictionary (perf-audit-2 04/P16): `values[indices[row]]`, no array of a string per row,
 * and views and windows of it gather the indices (`Column.storage`). The same values as the array of strings; a mask
 * whose rows do not all read as `''` keeps that array (`undefined`).
 */
function dictionaryField(
  column: EncodedColumn,
  indices: ArrayLike<number>,
  values: ArrayLike<string>,
  mask: ArrayLike<number> | undefined,
): Data.CifField | undefined {
  const rowCount = indices.length;
  if (mask) {
    for (let i = 0; i < rowCount; i++) {
      if (mask[i] !== Column.ValueKinds.Present && values[indices[i]] !== '') return;
    }
  }
  const str: Data.CifField['str'] = (row) => values[indices[row]];
  const int: Data.CifField['int'] = (row) => {
    const v = values[indices[row]];
    return fastParseInt(v, 0, v.length);
  };
  const float: Data.CifField['float'] = (row) => {
    const v = values[indices[row]];
    return fastParseFloat(v, 0, v.length);
  };
  const storage: Column.Storage = { kind: 'dictionary', indices, values };
  return {
    storage: () => storage,
    binaryEncoding: column.data.encoding,
    isDefined: true,
    rowCount,
    str,
    int,
    float,
    valueKind: mask ? (row) => mask[row] as Column.ValueKind : (row) => Column.ValueKinds.Present,
    // `values` may repeat a string
    areValuesEqual: (rowA, rowB) => {
      const a = indices[rowA];
      const b = indices[rowB];
      return a === b || values[a] === values[b];
    },
    toStringArray: (params) => ColumnHelpers.createAndFillArray(rowCount, str, params),
    toIntArray: (params) => ColumnHelpers.createAndFillArray(rowCount, int, params),
    toFloatArray: (params) => ColumnHelpers.createAndFillArray(rowCount, float, params),
  };
}
