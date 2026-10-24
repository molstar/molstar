/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, type RuntimeContext } from '@molstar/core/task';
import { CifCategory, CifField } from '../cif/data-model.js';
import { toTable } from '../cif/schema.js';
import type { StringLike } from '@molstar/core/util/string-like';
import { TextByteTokenBuilder, TextByteTokenizer, type TextByteTokens } from '../common/text/text-byte-tokenizer.js';
import { ReaderResult as Result } from '../result.js';
import { DynamoTblSchema, type DynamoTblTable } from './schema.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

export type DynamoTblFieldName = keyof DynamoTblSchema;

export interface DynamoTblFile {
  readonly rowCount: number;
  readonly fieldNames: ReadonlyArray<DynamoTblFieldName>;
  readonly fields: DynamoTblTable;
}

const RequiredColumnCount = 26;

/**
 * Dynamo TBL column layout: 1-based column number per Dynamo spec maps to the
 * 0-based source-token index used by the parser. Column numbers 33, 38, 39,
 * and 40 are intentionally absent from the Dynamo spec.
 */
const DynamoTblColumns: ReadonlyArray<{ readonly index: number; readonly name: DynamoTblFieldName }> = [
  { index: 0, name: 'tag' },
  { index: 1, name: 'aligned' },
  { index: 2, name: 'averaged' },
  { index: 3, name: 'dx' },
  { index: 4, name: 'dy' },
  { index: 5, name: 'dz' },
  { index: 6, name: 'tdrot' },
  { index: 7, name: 'tilt' },
  { index: 8, name: 'narot' },
  { index: 9, name: 'cc' },
  { index: 10, name: 'cc2' },
  { index: 11, name: 'cpu' },
  { index: 12, name: 'ftype' },
  { index: 13, name: 'ymintilt' },
  { index: 14, name: 'ymaxtilt' },
  { index: 15, name: 'xmintilt' },
  { index: 16, name: 'xmaxtilt' },
  { index: 17, name: 'fs1' },
  { index: 18, name: 'fs2' },
  { index: 19, name: 'tomo' },
  { index: 20, name: 'reg' },
  { index: 21, name: 'class' },
  { index: 22, name: 'annotation' },
  { index: 23, name: 'x' },
  { index: 24, name: 'y' },
  { index: 25, name: 'z' },
  { index: 26, name: 'dshift' },
  { index: 27, name: 'daxis' },
  { index: 28, name: 'dnarot' },
  { index: 29, name: 'dcc' },
  { index: 30, name: 'otag' },
  { index: 31, name: 'npar' },
  { index: 33, name: 'ref' },
  { index: 34, name: 'sref' },
  { index: 35, name: 'apix' },
  { index: 36, name: 'def' },
  { index: 40, name: 'eig1' },
  { index: 41, name: 'eig2' },
];

async function parseInternal(data: ReadonlyUint8Array, ctx: RuntimeContext) {
  const tokenizer = TextByteTokenizer(data);
  const warnings: string[] = [];
  let skippedRows = 0;
  let prevPosition = 0;
  let rowCount = 0;
  let expectedTokenCount = -1;

  // One Tokens builder per spec column. Each parsed row appends one
  // (start, end) pair to every column; columns absent (i.e. column index
  // beyond the row's token count) get an empty (0, 0) range so that
  // CifField.ofTextByteTokens reports them as ValueKinds.NotPresent.
  const sizeHint = Math.max(10, (data.length / 80) | 0);
  const columnTokens: TextByteTokens[] = DynamoTblColumns.map(() => TextByteTokenBuilder.create(data, sizeHint));

  // Reverse map: token index in a row -> index into DynamoTblColumns (or -1).
  const maxColumnIndex = DynamoTblColumns[DynamoTblColumns.length - 1].index;
  const tokenIndexToColumn = new Int32Array(maxColumnIndex + 1);
  for (let i = 0; i < tokenIndexToColumn.length; ++i) tokenIndexToColumn[i] = -1;
  for (let i = 0; i < DynamoTblColumns.length; ++i) tokenIndexToColumn[DynamoTblColumns[i].index] = i;

  // The (start, end) of the row's tokens up to the last spec column: added to the columns only once the row is
  // accepted, so a skipped (incomplete) row adds nothing.
  const rowTokens = new Uint32Array(2 * (maxColumnIndex + 1));

  while (tokenizer.position < tokenizer.length) {
    if (tokenizer.position - prevPosition > 100000 && ctx.shouldUpdate) {
      prevPosition = tokenizer.position;
      await ctx.update({ current: tokenizer.position, max: tokenizer.length });
    }

    // Read the next line.
    TextByteTokenizer.markStart(tokenizer);
    TextByteTokenizer.eatLine(tokenizer);
    const lineStart = tokenizer.tokenStart;
    const lineEnd = tokenizer.tokenEnd;
    if (lineStart >= lineEnd) continue;

    // Skip leading whitespace and detect empty/comment lines.
    let i = lineStart;
    while (i < lineEnd) {
      const c = data[i];
      if (c !== 9 && c !== 32) break;
      i++;
    }
    if (i >= lineEnd) continue;
    const firstChar = data[i];
    if (firstChar === 35 /* # */ || firstChar === 37 /* % */ || firstChar === 59 /* ; */) continue;

    // Scan whitespace-separated tokens within the line.
    let tokenCount = 0;
    while (i < lineEnd) {
      let c = data[i];
      while (i < lineEnd && (c === 9 || c === 32)) {
        i++;
        c = data[i];
      }
      if (i >= lineEnd) break;
      const ts = i;
      while (i < lineEnd) {
        const cc = data[i];
        if (cc === 9 || cc === 32) break;
        i++;
      }
      if (tokenCount <= maxColumnIndex) {
        rowTokens[2 * tokenCount] = ts;
        rowTokens[2 * tokenCount + 1] = i;
      }
      tokenCount++;
    }

    if (expectedTokenCount === -1) {
      if (tokenCount < RequiredColumnCount) {
        skippedRows += 1;
        continue;
      }
      expectedTokenCount = tokenCount;
    } else if (tokenCount !== expectedTokenCount) {
      throw new Error(`Inconsistent Dynamo TBL row width: expected ${expectedTokenCount} columns, got ${tokenCount}.`);
    }

    // Each spec column gets the row's token, or an empty range when its index lies beyond the row width.
    for (let c = 0; c < DynamoTblColumns.length; ++c) {
      const index = DynamoTblColumns[c].index;
      if (index < tokenCount) TextByteTokenBuilder.add(columnTokens[c], rowTokens[2 * index], rowTokens[2 * index + 1]);
      else TextByteTokenBuilder.add(columnTokens[c], 0, 0);
    }
    rowCount++;
  }

  if (rowCount === 0) {
    return Result.error<DynamoTblFile>('No readable Dynamo table rows were found.');
  }

  if (skippedRows > 0) {
    warnings.push(`Skipped ${skippedRows} Dynamo row${skippedRows === 1 ? '' : 's'} that were incomplete.`);
  }

  const cifFields: { [name: string]: CifField } = {};
  const fieldNames: DynamoTblFieldName[] = [];
  for (let c = 0; c < DynamoTblColumns.length; ++c) {
    const name = DynamoTblColumns[c].name;
    cifFields[name] = CifField.ofTextByteTokens(columnTokens[c]);
    fieldNames.push(name);
  }
  const category = CifCategory.ofFields('dynamo_tbl', cifFields);
  const fields = toTable(DynamoTblSchema, category) as DynamoTblTable;

  return Result.success<DynamoTblFile>(
    {
      rowCount,
      fieldNames,
      fields,
    },
    warnings,
  );
}

/** Parses a Dynamo table from its text or its UTF-8 bytes (a string is encoded once). */
export function parseDynamoTbl(data: StringLike | ReadonlyUint8Array) {
  return Task.create<Result<DynamoTblFile>>('Parse Dynamo TBL', async (ctx) => {
    return await parseInternal(TextByteTokenizer.bytesOf(data), ctx);
  });
}
