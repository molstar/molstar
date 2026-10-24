/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import { parseFloatBytes, parseIntBytes } from '@molstar/core/util/number-parser';
import { trimBytesStr, type TextByteTokens } from '../text-byte-tokenizer.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

/**
 * `FixedColumnProvider` (`fixed.ts`) over byte tokens of lines. Offsets and widths count bytes, which are the chars of
 * the fixed-column formats (ASCII).
 */
export function TextByteFixedColumnProvider(lines: TextByteTokens) {
  return function <T extends Column.Schema>(offset: number, width: number, type: T) {
    return TextByteFixedColumn(lines, offset, width, type);
  };
}

function skipSpaces(data: ReadonlyUint8Array, start: number, end: number) {
  while (start < end && data[start] === 32) start++;
  return start;
}

export function TextByteFixedColumn<T extends Column.Schema>(
  lines: TextByteTokens,
  offset: number,
  width: number,
  schema: T,
): Column<T['T']> {
  const { data, indices, count: rowCount } = lines;
  const { valueType: type } = schema;

  const value: Column<T['T']>['value'] =
    type === 'str'
      ? (row) => {
          const s = indices[2 * row] + offset,
            le = indices[2 * row + 1];
          if (s >= le) return '';
          let e = s + width;
          if (e > le) e = le;
          return trimBytesStr(data, s, e);
        }
      : type === 'int'
        ? (row) => {
            const s = indices[2 * row] + offset;
            const le = indices[2 * row + 1];
            if (s > le) return 0;
            // to the line's end: past the data's end a number read 0 (`fixed.ts` still does)
            const e = Math.min(s + width, le);
            return parseIntBytes(data, skipSpaces(data, s, e), e);
          }
        : (row) => {
            const s = indices[2 * row] + offset;
            const le = indices[2 * row + 1];
            if (s > le) return 0;
            const e = Math.min(s + width, le);
            return parseFloatBytes(data, skipSpaces(data, s, e), e);
          };
  return {
    schema: schema,
    isDefined: true,
    rowCount,
    value,
    valueKind: (row) => Column.ValueKinds.Present,
    toArray: (params) => ColumnHelpers.createAndFillArray(rowCount, value, params),
    areValuesEqual: (rowA, rowB) => value(rowA) === value(rowB),
  };
}
