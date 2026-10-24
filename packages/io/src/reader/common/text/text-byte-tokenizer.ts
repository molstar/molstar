/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * `Tokenizer` (`tokenizer.ts`) over the UTF-8 bytes of a text: the same functions, reading a `Uint8Array` instead of
 * `charCodeAt` on a string. The delimiters of the text formats are ASCII, so scanning bytes is safe, and values decode
 * as UTF-8 when read (`getTokenString`, `TextByteTokenColumn`). Why: a downloaded or read file is its bytes
 * (`Data.Response`); a parser over them neither decodes the whole text nor holds a string of it, and the scan is about
 * 2x faster than over a string (`.v6/designs/cif-parsing.md`; the CIF tokenizer is the byte tokenizer of its own).
 *
 * Offsets are byte offsets (Uint32 in the tokens): a text of 4 GiB or more is refused (`TextByteTokenizer.bytesOf`).
 * The bytes are `ReadonlyUint8Array` (core `type-helpers`): a parser only reads them, so the plugin hands it a response's
 * shared bytes without a copy (`SO.Data.readBytes`), and the types catch a write.
 */

import type { RuntimeContext } from '@molstar/core/task/execution/runtime-context';
import { chunkedSubtask } from '@molstar/core/task/util/chunked';
import type { StringLike } from '@molstar/core/util/string-like';
import { utf8Encode, utf8ReadValue } from '@molstar/core/util/utf8';
import { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

export { TextByteTokenizer };

interface TextByteTokenizer {
  data: ReadonlyUint8Array;

  position: number;
  length: number;

  lineNumber: number;
  tokenStart: number;
  tokenEnd: number;
}

/** Ranges of the bytes of a text, a start and an end per token: `Tokens` over bytes. */
export interface TextByteTokens {
  data: ReadonlyUint8Array;
  count: number;
  indices: ArrayLike<number>;
}

/** The largest text the byte tokenizer reads: its offsets are Uint32. */
const MaxLength = 0xffffffff;

function TextByteTokenizer(data: ReadonlyUint8Array): TextByteTokenizer {
  TextByteTokenizer.checkLength(data);
  return {
    data,
    position: 0,
    length: data.length,
    lineNumber: 1,
    tokenStart: 0,
    tokenEnd: 0,
  };
}

namespace TextByteTokenizer {
  /** Throws for a text of 4 GiB or more: the offsets are Uint32 and would wrap around. */
  export function checkLength(data: ReadonlyUint8Array) {
    if (data.length > MaxLength) throw new Error('Text of 4 GiB or more is not supported.');
  }

  /**
   * The bytes a parser reads: the bytes themselves, or a string encoded to UTF-8 (one pass, a copy), without a leading
   * UTF-8 byte order mark (a view past it), as decoding to a string dropped it. A parser that takes either calls it once
   * at the start.
   */
  export function bytesOf(data: StringLike | ReadonlyUint8Array): ReadonlyUint8Array {
    let bytes = ReadonlyUint8Array.is(data) ? data : utf8Encode(data);
    if (bytes[0] === 0xef && bytes[1] === 0xbb && bytes[2] === 0xbf) bytes = bytes.subarray(3);
    checkLength(bytes);
    return bytes;
  }

  export function getTokenString(state: TextByteTokenizer) {
    return utf8ReadValue(state.data, state.tokenStart, state.tokenEnd - state.tokenStart);
  }

  /** Resets the state */
  export function reset(state: TextByteTokenizer) {
    state.position = 0;
    state.lineNumber = 1;
    state.tokenStart = 0;
    state.tokenEnd = 0;
  }

  /**
   * Eat everything until a newline occurs.
   */
  export function eatLine(state: TextByteTokenizer): boolean {
    // locals, not the state's fields, per byte: the hot loop of a large file
    const { data, length } = state;
    let position = state.position;
    while (position < length) {
      const c = data[position];
      if (c === 10) {
        // \n
        state.tokenEnd = position;
        state.position = position + 1;
        ++state.lineNumber;
        return true;
      }
      if (c === 13) {
        // \r
        state.tokenEnd = position;
        ++position;
        ++state.lineNumber;
        if (data[position] === 10) ++position;
        state.position = position;
        return true;
      }
      ++position;
    }
    state.position = position;
    state.tokenEnd = position;
    return state.tokenStart !== state.tokenEnd;
  }

  /** Sets the current token start to the current position */
  export function markStart(state: TextByteTokenizer) {
    state.tokenStart = state.position;
  }

  /** Sets the current token start to current position and moves to the next line. */
  export function markLine(state: TextByteTokenizer) {
    state.tokenStart = state.position;
    return eatLine(state);
  }

  /** Advance the state and return line as string. */
  export function readLine(state: TextByteTokenizer): string {
    markLine(state);
    return getTokenString(state);
  }

  /** Advance the state and return trimmed line as string. */
  export function readLineTrim(state: TextByteTokenizer): string {
    markLine(state);
    const position = state.position;
    trim(state, state.tokenStart, state.tokenEnd);
    state.position = position;
    return getTokenString(state);
  }

  function readLinesChunk(state: TextByteTokenizer, count: number, tokens: TextByteTokens) {
    let read = 0;
    for (let i = 0; i < count; i++) {
      if (!markLine(state)) return read;
      TextByteTokenBuilder.addUnchecked(tokens, state.tokenStart, state.tokenEnd);
      read++;
    }
    return read;
  }

  /** Advance the state by the given number of lines and return them*/
  export function markLines(state: TextByteTokenizer, count: number): TextByteTokens {
    const lineTokens = TextByteTokenBuilder.create(state.data, count * 2);
    readLinesChunk(state, count, lineTokens);
    return lineTokens;
  }

  /** Advance the state by the given number of lines and return them */
  export function readLines(state: TextByteTokenizer, count: number): string[] {
    const ret: string[] = [];
    for (let i = 0; i < count; i++) {
      ret.push(readLine(state));
    }
    return ret;
  }

  /** Advance the state by the given number of lines and return line starts/ends as tokens. */
  export async function readLinesAsync(
    state: TextByteTokenizer,
    count: number,
    ctx: RuntimeContext,
    initialLineCount = 100000,
  ): Promise<TextByteTokens> {
    const lineTokens = TextByteTokenBuilder.create(state.data, count * 2);

    let linesAlreadyRead = 0;
    await chunkedSubtask(
      ctx,
      initialLineCount,
      state,
      (chunkSize, state) => {
        const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
        readLinesChunk(state, linesToRead, lineTokens);
        linesAlreadyRead += linesToRead;
        return linesToRead;
      },
      (ctx, state) => ctx.update({ message: 'Parsing...', current: state.position, max: state.length }),
    );

    return lineTokens;
  }

  export function readAllLines(data: ReadonlyUint8Array) {
    const state = TextByteTokenizer(data);
    const tokens = TextByteTokenBuilder.create(state.data, Math.max(data.length / 80, 2));
    while (markLine(state)) {
      TextByteTokenBuilder.add(tokens, state.tokenStart, state.tokenEnd);
    }
    return tokens;
  }

  function readLinesChunkChecked(state: TextByteTokenizer, count: number, tokens: TextByteTokens) {
    let read = 0;
    for (let i = 0; i < count; i++) {
      if (!markLine(state)) return read;
      TextByteTokenBuilder.add(tokens, state.tokenStart, state.tokenEnd);
      read++;
    }
    return read;
  }

  export async function readAllLinesAsync(data: ReadonlyUint8Array, ctx: RuntimeContext, chunkSize = 100000) {
    const state = TextByteTokenizer(data);
    const tokens = TextByteTokenBuilder.create(state.data, Math.max(data.length / 80, 2));

    await chunkedSubtask(
      ctx,
      chunkSize,
      state,
      (chunkSize, state) => {
        readLinesChunkChecked(state, chunkSize, tokens);
        return state.position < state.length ? chunkSize : 0;
      },
      (ctx, state) => ctx.update({ message: 'Parsing...', current: state.position, max: state.length }),
    );

    return tokens;
  }

  /**
   * Eat everything until a whitespace/newline occurs.
   */
  export function eatValue(state: TextByteTokenizer) {
    const { data, length } = state;
    let position = state.position;
    while (position < length) {
      const c = data[position];
      if (c <= 32 && (c === 32 || c === 10 || c === 13 || c === 9)) break;
      ++position;
    }
    state.position = position;
    state.tokenEnd = position;
  }

  /**
   * Skips all the whitespace - space, tab, newline, CR
   * Handles incrementing line count.
   */
  export function skipWhitespace(state: TextByteTokenizer): number {
    const { data, length } = state;
    let prev = -1;
    while (state.position < length) {
      const c = data[state.position];
      switch (c) {
        case 9: // '\t'
        case 32: // ' '
          prev = c;
          ++state.position;
          break;
        case 10: // \n
          // handle \r\n
          if (prev !== 13) {
            ++state.lineNumber;
          }
          prev = c;
          ++state.position;
          break;
        case 13: // \r
          prev = c;
          ++state.position;
          ++state.lineNumber;
          break;
        default:
          return prev;
      }
    }
    return prev;
  }

  /** Skips all the whitespace */
  export function skipStrictWhitespace(state: TextByteTokenizer): number {
    const { data, length } = state;
    let prev = -1;
    while (state.position < length) {
      const c = data[state.position];
      switch (c) {
        case 9: // '\t'
        case 32: // ' '
          prev = c;
          ++state.position;
          break;
        default:
          return prev;
      }
    }
    return prev;
  }

  /** Trims spaces and tabs */
  export function trim(state: TextByteTokenizer, start: number, end: number) {
    const { data } = state;
    let s = start,
      e = end - 1;

    let c = data[s];
    while ((c === 9 || c === 32) && s <= e) c = data[++s];
    c = data[e];
    while ((c === 9 || c === 32) && e >= s) c = data[--e];

    state.tokenStart = s;
    state.tokenEnd = e + 1;
    state.position = end;
    return state;
  }
}

/** The bytes `start` to `end` without the spaces and tabs around them, decoded: `trimStr` over bytes. */
export function trimBytesStr(data: ReadonlyUint8Array, start: number, end: number) {
  let s = start,
    e = end - 1;
  let c = data[s];
  while ((c === 9 || c === 32) && s <= e) c = data[++s];
  c = data[e];
  while ((c === 9 || c === 32) && e >= s) c = data[--e];
  return utf8ReadValue(data, s, e + 1 - s);
}

export namespace TextByteTokenBuilder {
  interface Builder extends TextByteTokens {
    offset: number;
    indices: Uint32Array;
    indicesLenMinus2: number;
  }

  function resize(builder: Builder) {
    // scale the size using golden ratio, because why not.
    const newBuffer = new Uint32Array((1.61 * builder.indices.length) | 0);
    newBuffer.set(builder.indices);
    builder.indices = newBuffer;
    builder.indicesLenMinus2 = (newBuffer.length - 2) | 0;
  }

  export function add(tokens: TextByteTokens, start: number, end: number) {
    const builder = tokens as Builder;
    if (builder.offset > builder.indicesLenMinus2) {
      resize(builder);
    }
    builder.indices[builder.offset++] = start;
    builder.indices[builder.offset++] = end;
    tokens.count++;
  }

  export function addToken(tokens: TextByteTokens, tokenizer: TextByteTokenizer) {
    add(tokens, tokenizer.tokenStart, tokenizer.tokenEnd);
  }

  export function addUnchecked(tokens: TextByteTokens, start: number, end: number) {
    (tokens as Builder).indices[(tokens as Builder).offset++] = start;
    (tokens as Builder).indices[(tokens as Builder).offset++] = end;
    tokens.count++;
  }

  export function create(data: ReadonlyUint8Array, size: number): TextByteTokens {
    size = Math.max(10, size);
    return <Builder>{
      data,
      indicesLenMinus2: (size - 2) | 0,
      count: 0,
      offset: 0,
      indices: new Uint32Array(size),
    };
  }
}
