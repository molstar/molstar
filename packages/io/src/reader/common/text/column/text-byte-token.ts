/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import { parseFloatBytes, parseIntBytes } from '@molstar/core/util/number-parser';
import { utf8ReadValue } from '@molstar/core/util/utf8';
import type { TextByteTokens } from '../text-byte-tokenizer.js';
import { ByteStringCache } from '../byte-string-cache.js';

/** `TokenColumnProvider` (`token.ts`) over byte tokens. */
export function TextByteTokenColumnProvider(tokens: TextByteTokens) {
  return function <T extends Column.Schema>(type: T) {
    return TextByteTokenColumn(tokens, type);
  };
}

/**
 * `TokenColumn` over byte tokens: values decode as UTF-8, short ones through a string cache (`ByteStringCache`);
 * numbers parse from the bytes.
 */
export function TextByteTokenColumn<T extends Column.Schema>(tokens: TextByteTokens, schema: T): Column<T['T']> {
  const { data, indices, count: rowCount } = tokens;
  const { valueType: type } = schema;
  let cache: ByteStringCache | undefined;
  const cached = rowCount >= ByteStringCache.MinRows;

  const value: Column<T['T']>['value'] =
    type === 'str'
      ? (row) => {
          const s = indices[2 * row];
          const l = indices[2 * row + 1] - s;
          if (l > ByteStringCache.MaxLength || !cached) return utf8ReadValue(data, s, l);
          return ByteStringCache.read((cache ??= ByteStringCache()), data, s, l);
        }
      : type === 'int'
        ? (row) => parseIntBytes(data, indices[2 * row], indices[2 * row + 1]) || 0
        : (row) => parseFloatBytes(data, indices[2 * row], indices[2 * row + 1]) || 0;

  return {
    schema: schema,
    isDefined: true,
    rowCount,
    value,
    valueKind: (row) => Column.ValueKinds.Present,
    toArray: (params) => ColumnHelpers.createAndFillArray(rowCount, value, params),
    areValuesEqual: areTextByteValuesEqualProvider(tokens),
  };
}

export function areTextByteValuesEqualProvider(tokens: TextByteTokens) {
  const { data, indices } = tokens;
  return function (rowA: number, rowB: number) {
    const aS = indices[2 * rowA],
      bS = indices[2 * rowB];
    const len = indices[2 * rowA + 1] - aS;
    if (len !== indices[2 * rowB + 1] - bS) return false;
    for (let i = 0; i < len; i++) {
      if (data[i + aS] !== data[i + bS]) {
        return false;
      }
    }
    return true;
  };
}

export function areTextByteTokensEmpty(tokens: TextByteTokens) {
  const { count, indices } = tokens;
  for (let i = 0; i < count; ++i) {
    if (indices[2 * i] !== indices[2 * i + 1]) return false;
  }
  return true;
}
