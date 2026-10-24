/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * CSV parser, RFC 4180 (https://www.rfc-editor.org/rfc/rfc4180), over the UTF-8 bytes of the text:
 * - a field is the bytes between delimiters, spaces included (no trimming; numbers parse past spaces around them);
 * - a field starting with the quote is quoted: it may contain delimiters, line breaks and the quote doubled (`""`),
 *   read as one quote; only spaces and tabs may follow its closing quote;
 * - records end at CRLF, LF or CR; the last line break is optional; empty lines are skipped (so a record is never
 *   empty: a one-column file reads an empty value as `""`);
 * - the first record is the column names (unless `noColumnNames`) and sets the column count: a shorter record reads
 *   the missing fields as empty (not present), the fields of a longer one past the count are dropped;
 * - a leading UTF-8 byte order mark is skipped;
 * - not in the RFC: lines starting with `comment`, when set, are skipped (off by default).
 * An unterminated quoted field, or text after a closing quote, is an error with its line.
 */

import { Task, type RuntimeContext, chunkedSubtask } from '@molstar/core/task';
import type { StringLike } from '@molstar/core/util/string-like';
import { TextByteTokenBuilder, TextByteTokenizer, type TextByteTokens } from '../common/text/text-byte-tokenizer.js';
import { ReaderResult as Result } from '../result.js';
import * as Data from './data-model.js';
import { CsvField } from './field.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

interface CsvOptions {
  /** One ASCII char; `"` by default. */
  quote: string;
  /** One ASCII char starting a comment line, `''` (the default) for none: comments are not in RFC 4180. */
  comment: string;
  /** One ASCII char; `,` by default. */
  delimiter: string;
  /** The first record is data: the columns are named by their index (`'0'`, `'1'`, ...). */
  noColumnNames: boolean;
}

const DefaultOptions: CsvOptions = { quote: '"', comment: '', delimiter: ',', noColumnNames: false };

class CsvError {
  constructor(
    public message: string,
    public line: number,
  ) {}
}

interface State {
  data: ReadonlyUint8Array;
  length: number;
  position: number;
  lineNumber: number;

  quote: number;
  delimiter: number;
  /** -1 for none */
  comment: number;

  /** ranges of the fields of the record read last, a start and an end each */
  fields: number[];
  fieldCount: number;
}

function asciiCode(char: string, name: string) {
  const c = char.charCodeAt(0);
  if (char.length !== 1 || c > 127) throw new Error(`CSV ${name} must be one ASCII char, got '${char}'.`);
  return c;
}

function State(data: ReadonlyUint8Array, opts: CsvOptions): State {
  const quote = asciiCode(opts.quote, 'quote');
  const delimiter = asciiCode(opts.delimiter, 'delimiter');
  const comment = opts.comment ? asciiCode(opts.comment, 'comment') : -1;
  if (quote === delimiter || quote === 10 || quote === 13 || delimiter === 10 || delimiter === 13) {
    throw new Error('CSV quote and delimiter must differ and not be a line break.');
  }
  if (comment === quote || comment === delimiter || comment === 10 || comment === 13) {
    throw new Error('CSV comment must differ from the quote and the delimiter and not be a line break.');
  }
  const bom = data[0] === 0xef && data[1] === 0xbb && data[2] === 0xbf;
  return {
    data,
    length: data.length,
    position: bom ? 3 : 0,
    lineNumber: 1,
    quote,
    delimiter,
    comment,
    fields: [],
    fieldCount: 0,
  };
}

/** Moves past the line break at the position, if any. */
function eatLineBreak(state: State) {
  const c = state.data[state.position];
  if (c === 13) {
    ++state.position;
    if (state.data[state.position] === 10) ++state.position;
    ++state.lineNumber;
  } else if (c === 10) {
    ++state.position;
    ++state.lineNumber;
  }
}

/** Skips empty lines and comment lines; false at the end of the data. */
function skipToRecord(state: State) {
  const { data, length, comment } = state;
  while (state.position < length) {
    const c = data[state.position];
    if (c === 10 || c === 13) {
      eatLineBreak(state);
    } else if (c === comment) {
      while (state.position < length && data[state.position] !== 10 && data[state.position] !== 13) ++state.position;
    } else {
      return true;
    }
  }
  return false;
}

function addField(state: State, start: number, end: number) {
  const i = 2 * state.fieldCount++;
  state.fields[i] = start;
  state.fields[i + 1] = end;
}

/**
 * Reads the quoted field at the position. Its range is the value between the quotes, or the field with its quotes
 * when it has a doubled quote (`CsvField` unescapes such values: no unquoted value starts with the quote).
 */
function readQuoted(state: State) {
  const { data, length, quote } = state;
  const open = state.position;
  const line = state.lineNumber;
  let position = open + 1;
  let escaped = false;
  while (true) {
    if (position >= length) throw new CsvError('Unterminated quoted field.', line);
    const c = data[position];
    if (c === quote) {
      if (data[position + 1] === quote) {
        escaped = true;
        position += 2;
        continue;
      }
      break;
    }
    if (c === 10) {
      // \n, not of a \r\n
      if (data[position - 1] !== 13) ++state.lineNumber;
    } else if (c === 13) {
      ++state.lineNumber;
    }
    ++position;
  }
  const close = position;
  ++position;
  // spaces and tabs, unless one is the delimiter
  while (position < length) {
    const c = data[position];
    if ((c !== 32 && c !== 9) || c === state.delimiter) break;
    ++position;
  }
  if (position < length) {
    const c = data[position];
    if (c !== state.delimiter && c !== 10 && c !== 13) {
      throw new CsvError('Unexpected text after a closing quote.', state.lineNumber);
    }
  }
  state.position = position;
  if (escaped) addField(state, open, close + 1);
  else addField(state, open + 1, close);
}

/** Reads the record at the position (not an empty or comment line) into `state.fields`. */
function readRecord(state: State) {
  const { data, length, quote, delimiter } = state;
  state.fieldCount = 0;
  while (true) {
    if (data[state.position] === quote) {
      readQuoted(state);
    } else {
      // locals, not the state's fields, per byte: the hot loop of a large file
      const start = state.position;
      let position = start;
      while (position < length) {
        const c = data[position];
        if (c === delimiter || c === 10 || c === 13) break;
        ++position;
      }
      state.position = position;
      addField(state, start, position);
    }
    if (state.position < length && data[state.position] === delimiter) {
      ++state.position;
      continue;
    }
    eatLineBreak(state);
    return;
  }
}

function readRecordsChunk(chunkSize: number, ctx: { state: State; columns: TextByteTokens[]; rowCount: number }) {
  const { state, columns } = ctx;
  const columnCount = columns.length;
  let read = 0;
  while (read < chunkSize && skipToRecord(state)) {
    readRecord(state);
    const { fields, fieldCount } = state;
    for (let i = 0; i < columnCount; i++) {
      if (i < fieldCount) TextByteTokenBuilder.add(columns[i], fields[2 * i], fields[2 * i + 1]);
      else TextByteTokenBuilder.add(columns[i], 0, 0);
    }
    ++read;
  }
  ctx.rowCount += read;
  return read;
}

async function parseInternal(
  data: ReadonlyUint8Array,
  ctx: RuntimeContext,
  opts: CsvOptions,
): Promise<Result<Data.CsvFile>> {
  const state = State(data, opts);
  ctx.update({ message: 'Parsing...', current: 0, max: data.length });

  if (!skipToRecord(state)) return Result.success(Data.CsvFile(Data.CsvTable(0, [], Object.create(null))));

  const firstStart = state.position;
  const firstLine = state.lineNumber;
  readRecord(state);
  const columnCount = state.fieldCount;
  const columnNames: string[] = [];
  for (let i = 0; i < columnCount; i++) {
    columnNames.push(
      opts.noColumnNames ? `${i}` : CsvField.decode(data, state.fields[2 * i], state.fields[2 * i + 1], state.quote),
    );
  }
  // the rows estimated by the bytes of the first record, at most a million to start with (the tokens grow)
  const rowEstimate = Math.min(Math.ceil(data.length / Math.max(state.position - firstStart, 8)), 1 << 20);
  if (opts.noColumnNames) {
    state.position = firstStart;
    state.lineNumber = firstLine;
  }

  const columns: TextByteTokens[] = [];
  for (let i = 0; i < columnCount; i++) columns.push(TextByteTokenBuilder.create(data, 2 * rowEstimate));

  const chunkState = { state, columns, rowCount: 0 };
  await chunkedSubtask(ctx, 100000, chunkState, readRecordsChunk, (ctx, { state }) =>
    ctx.update({ message: 'Parsing...', current: state.position, max: state.length }),
  );

  const table: Data.CsvColumns = Object.create(null);
  for (let i = 0; i < columnCount; ++i) {
    table[columnNames[i]] = CsvField(columns[i], state.quote);
  }
  return Result.success(Data.CsvFile(Data.CsvTable(chunkState.rowCount, columnNames, table)));
}

/** Parses CSV from its UTF-8 bytes, or from a string, encoded to them first. */
export function parseCsv(data: StringLike | ReadonlyUint8Array, opts?: Partial<CsvOptions>) {
  const completeOpts = { ...DefaultOptions, ...opts };
  return Task.create<Result<Data.CsvFile>>('Parse CSV', async (ctx) => {
    try {
      return await parseInternal(TextByteTokenizer.bytesOf(data), ctx, completeOpts);
    } catch (e) {
      if (e instanceof CsvError) return Result.error(e.message, e.line);
      throw e;
    }
  });
}
