/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '@molstar/core/data/db';
import { Task } from '@molstar/core/task';
import { parseFloatBytes, parseIntBytes } from '@molstar/core/util/number-parser';
import type { StringLike } from '@molstar/core/util/string-like';
import { TextByteTokenColumnProvider as TokenColumn } from '../common/text/column/text-byte-token.js';
import { TextByteTokenBuilder, TextByteTokenizer } from '../common/text/text-byte-tokenizer.js';
import { ReaderResult as Result } from '../result.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

export interface XyzFile {
  readonly molecules: {
    readonly comment: string;
    readonly count: number;
    readonly x: Column<number>;
    readonly y: Column<number>;
    readonly z: Column<number>;
    readonly type_symbol: Column<string>;
  }[];
}

/** The next space- or tab-separated field of the bytes `position` to `end`: its start, its end in `fieldEnd`. */
let fieldEnd = 0;
function nextField(data: ReadonlyUint8Array, position: number, end: number) {
  while (position < end && (data[position] === 32 || data[position] === 9)) ++position;
  let e = position;
  while (e < end && data[e] !== 32 && data[e] !== 9) ++e;
  fieldEnd = e;
  return position;
}

/** A coordinate field as a number; a missing field is NaN (as `+undefined` was). */
function coordinate(data: ReadonlyUint8Array, start: number, end: number) {
  return start === end ? NaN : parseFloatBytes(data, start, end);
}

/** `+line` of the atom count line: an integer from the bytes, the decoded line only when it is not a plain one. */
function readCount(tokenizer: TextByteTokenizer) {
  TextByteTokenizer.markLine(tokenizer);
  const { data, tokenEnd } = tokenizer;
  let s = tokenizer.tokenStart;
  while (s < tokenEnd && (data[s] === 32 || data[s] === 9)) s++;
  let e = tokenEnd;
  while (e > s && (data[e - 1] === 32 || data[e - 1] === 9)) e--;
  let i = s;
  while (i < e && data[i] >= 48 && data[i] <= 57) i++;
  if (i === e && e - s < 10) return s === e ? 0 : parseIntBytes(data, s, e);
  return +TextByteTokenizer.getTokenString(tokenizer);
}

function handleMolecule(tokenizer: TextByteTokenizer): XyzFile['molecules'][number] {
  let count = tokenizer.position >= tokenizer.data.length - 1 ? 0 : readCount(tokenizer);
  if (isNaN(count)) count = 0;

  const comment = TextByteTokenizer.readLine(tokenizer);

  const x = new Float64Array(count);
  const y = new Float64Array(count);
  const z = new Float64Array(count);
  const type_symbol = TextByteTokenBuilder.create(tokenizer.data, count * 2);

  const { data } = tokenizer;
  for (let i = 0; i < count; ++i) {
    TextByteTokenizer.markLine(tokenizer);
    const end = tokenizer.tokenEnd;
    let s = nextField(data, tokenizer.tokenStart, end);
    TextByteTokenBuilder.addUnchecked(type_symbol, s, fieldEnd);
    s = nextField(data, fieldEnd, end);
    x[i] = coordinate(data, s, fieldEnd);
    s = nextField(data, fieldEnd, end);
    y[i] = coordinate(data, s, fieldEnd);
    s = nextField(data, fieldEnd, end);
    z[i] = coordinate(data, s, fieldEnd);
  }

  return {
    count,
    comment,
    x: Column.ofFloatArray(x),
    y: Column.ofFloatArray(y),
    z: Column.ofFloatArray(z),
    type_symbol: TokenColumn(type_symbol)(Column.Schema.str),
  };
}

function parseInternal(data: ReadonlyUint8Array): Result<XyzFile> {
  const tokenizer = TextByteTokenizer(data);

  const molecules: XyzFile['molecules'] = [];
  while (true) {
    const mol = handleMolecule(tokenizer);
    if (mol.count === 0) break;
    molecules.push(mol);
  }

  const result: XyzFile = { molecules };
  return Result.success(result);
}

/** Parses XYZ from its text or its UTF-8 bytes (read, not modified). */
export function parseXyz(data: StringLike | ReadonlyUint8Array) {
  return Task.create<Result<XyzFile>>('Parse Mol', async () => {
    return parseInternal(TextByteTokenizer.bytesOf(data));
  });
}
