/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import { parseFloatBytes, parseIntBytes } from '@molstar/core/util/number-parser';
import { utf8ReadValue } from '@molstar/core/util/utf8';
import { areTextByteValuesEqualProvider } from '../common/text/column/text-byte-token.js';
import type { TextByteTokens } from '../common/text/text-byte-tokenizer.js';
import type { CsvColumn } from './data-model.js';
import { ByteStringCache } from '../common/text/byte-string-cache.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

function skipBlanks(data: ReadonlyUint8Array, start: number, end: number) {
  while (start < end && (data[start] === 32 || data[start] === 9)) start++;
  return start;
}

/**
 * A column of CSV values (`parser.ts`): a value is its bytes as they are; an empty value is not present (`.` and `?`
 * are values, unlike CIF). Short values decode through a string cache (`ByteStringCache`); numbers parse past the
 * spaces and tabs around them.
 */
export function CsvField(tokens: TextByteTokens, quote: number): CsvColumn {
  const { data, indices, count: rowCount } = tokens;

  let cache: ByteStringCache | undefined;
  const cached = rowCount >= ByteStringCache.MinRows;
  const str: CsvColumn['str'] = (row) => {
    const s = indices[2 * row],
      e = indices[2 * row + 1];
    // a quoted value with a doubled quote (its range starts with the quote) is unescaped, not cached
    if (e - s > ByteStringCache.MaxLength || !cached || data[s] === quote) return CsvField.decode(data, s, e, quote);
    return ByteStringCache.read((cache ??= ByteStringCache()), data, s, e - s);
  };
  const rawEqual = areTextByteValuesEqualProvider(tokens);
  // a value with a quote has two forms (`x"y` and `"x""y"`): compared by its string then
  const areValuesEqual: CsvColumn['areValuesEqual'] = (rowA, rowB) =>
    data[indices[2 * rowA]] === quote || data[indices[2 * rowB]] === quote
      ? str(rowA) === str(rowB)
      : rawEqual(rowA, rowB);
  const int: CsvColumn['int'] = (row) => {
    const e = indices[2 * row + 1];
    return parseIntBytes(data, skipBlanks(data, indices[2 * row], e), e) || 0;
  };
  const float: CsvColumn['float'] = (row) => {
    const e = indices[2 * row + 1];
    return parseFloatBytes(data, skipBlanks(data, indices[2 * row], e), e) || 0;
  };
  const valueKind: CsvColumn['valueKind'] = (row) =>
    indices[2 * row] === indices[2 * row + 1] ? Column.ValueKinds.NotPresent : Column.ValueKinds.Present;

  return {
    binaryEncoding: void 0,
    isDefined: true,
    rowCount,
    str,
    int,
    float,
    valueKind,
    areValuesEqual,
    toStringArray: (params) => ColumnHelpers.createAndFillArray(rowCount, str, params),
    toIntArray: (params) => ColumnHelpers.createAndFillArray(rowCount, int, params),
    toFloatArray: (params) => ColumnHelpers.createAndFillArray(rowCount, float, params),
  };
}

export namespace CsvField {
  /** The value at `start` to `end`: a range starting with the quote is a quoted value with a doubled quote. */
  export function decode(data: ReadonlyUint8Array, start: number, end: number, quote: number) {
    if (start < end && data[start] === quote) {
      const q = String.fromCharCode(quote);
      return utf8ReadValue(data, start + 1, end - start - 2).replaceAll(q + q, q);
    }
    return utf8ReadValue(data, start, end - start);
  }
}
