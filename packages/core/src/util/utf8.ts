/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ChunkedBigString, MAX_STRING_LENGTH, type StringLike } from './string-like.js';
import type { ReadonlyUint8Array } from './type-helpers.js';

export function utf8Write(data: Uint8Array, offset: number, str: string) {
  for (let i = 0, l = str.length; i < l; i++) {
    // A surrogate pair is one code point
    const codePoint = str.codePointAt(i)!;
    if (codePoint > 0xffff) i++;

    // One byte of UTF-8
    if (codePoint < 0x80) {
      data[offset++] = ((codePoint >>> 0) & 0x7f) | 0x00;
      continue;
    }

    // Two bytes of UTF-8
    if (codePoint < 0x800) {
      data[offset++] = ((codePoint >>> 6) & 0x1f) | 0xc0;
      data[offset++] = ((codePoint >>> 0) & 0x3f) | 0x80;
      continue;
    }

    // Three bytes of UTF-8.
    if (codePoint < 0x10000) {
      data[offset++] = ((codePoint >>> 12) & 0x0f) | 0xe0;
      data[offset++] = ((codePoint >>> 6) & 0x3f) | 0x80;
      data[offset++] = ((codePoint >>> 0) & 0x3f) | 0x80;
      continue;
    }

    // Four bytes of UTF-8
    if (codePoint < 0x110000) {
      data[offset++] = ((codePoint >>> 18) & 0x07) | 0xf0;
      data[offset++] = ((codePoint >>> 12) & 0x3f) | 0x80;
      data[offset++] = ((codePoint >>> 6) & 0x3f) | 0x80;
      data[offset++] = ((codePoint >>> 0) & 0x3f) | 0x80;
      continue;
    }
    throw new Error('bad codepoint ' + codePoint);
  }
}

const __chars = (function () {
  const data: string[] = [];
  for (let i = 0; i < 1024; i++) data[i] = String.fromCharCode(i);
  return data;
})();

function throwError(err: string) {
  throw new Error(err);
}

function _utf8Read(data: ReadonlyUint8Array, offset: number, length: number) {
  const chars = __chars;
  let str: string[] | undefined = void 0,
    chunkOffset = 0;
  const chunk: string[] = [],
    chunkSize = 512;

  for (let i = offset, end = offset + length; i < end; i++) {
    const byte = data[i];
    if ((byte & 0x80) === 0x00) {
      // One byte character
      chunk[chunkOffset++] = chars[byte];
    } else if ((byte & 0xe0) === 0xc0) {
      // Two byte character
      chunk[chunkOffset++] = chars[((byte & 0x0f) << 6) | (data[++i] & 0x3f)];
    } else if ((byte & 0xf0) === 0xe0) {
      // Three byte character
      chunk[chunkOffset++] = String.fromCharCode(
        ((byte & 0x0f) << 12) | ((data[++i] & 0x3f) << 6) | ((data[++i] & 0x3f) << 0),
      );
    } else if ((byte & 0xf8) === 0xf0) {
      // Four byte character
      chunk[chunkOffset++] = String.fromCharCode(
        ((byte & 0x07) << 18) | ((data[++i] & 0x3f) << 12) | ((data[++i] & 0x3f) << 6) | ((data[++i] & 0x3f) << 0),
      );
    } else throwError('Invalid byte ' + byte.toString(16));

    if (chunkOffset === chunkSize) {
      str = str || [];
      str[str.length] = chunk.join('');
      chunkOffset = 0;
    }
  }
  if (!str) return chunk.slice(0, chunkOffset).join('');
  if (chunkOffset > 0) {
    str[str.length] = chunk.slice(0, chunkOffset).join('');
  }
  return str.join('');
}

/**
 * Up to this many bytes, ASCII is made with `String.fromCharCode`: a `TextDecoder` call costs ~0.4 µs, about 10x a short
 * string, and a parser of text values (CIF fields) reads millions of them.
 */
const ShortAsciiLength = 24;

/** The string of `length` ASCII bytes at `offset`, or `undefined` when one is not ASCII or the range is out of bounds. */
function shortAscii(d: ReadonlyUint8Array, o: number, length: number): string | undefined {
  if (o + length > d.length) return;
  let mask = 0;
  for (let i = o, e = o + length; i < e; i++) mask |= d[i];
  if (mask & 0x80) return;
  switch (length) {
    case 0:
      return '';
    case 1:
      return String.fromCharCode(d[o]);
    case 2:
      return String.fromCharCode(d[o], d[o + 1]);
    case 3:
      return String.fromCharCode(d[o], d[o + 1], d[o + 2]);
    case 4:
      return String.fromCharCode(d[o], d[o + 1], d[o + 2], d[o + 3]);
  }
  let ret = '';
  for (; length >= 8; o += 8, length -= 8) {
    ret += String.fromCharCode(d[o], d[o + 1], d[o + 2], d[o + 3], d[o + 4], d[o + 5], d[o + 6], d[o + 7]);
  }
  for (; length >= 4; o += 4, length -= 4) ret += String.fromCharCode(d[o], d[o + 1], d[o + 2], d[o + 3]);
  for (; length > 0; o++, length--) ret += String.fromCharCode(d[o]);
  return ret;
}

const utf8Decoder = typeof TextDecoder !== 'undefined' ? new TextDecoder() : undefined;
/** Decode UTF8 data. Return as primitive `string` type, or fail if the result is longer than MAX_STRING_LENGTH. */
export function utf8Read(data: ReadonlyUint8Array, offset: number = 0, length: number = data.length): string {
  if (length <= ShortAsciiLength) {
    const ascii = shortAscii(data, offset, length);
    if (ascii !== void 0) return ascii;
  }
  if (utf8Decoder) {
    const input = offset || length !== data.length ? data.subarray(offset, offset + length) : data;
    return utf8Decoder.decode(input);
  } else {
    return _utf8Read(data, offset, length);
  }
}

const utf8ValueDecoder = typeof TextDecoder !== 'undefined' ? new TextDecoder('utf-8', { ignoreBOM: true }) : undefined;
/**
 * Decode UTF8 data of a value in a text, keeping a leading U+FEFF: `utf8Read` (a `TextDecoder`) drops it, as the byte
 * order mark of a file.
 */
export function utf8ReadValue(data: ReadonlyUint8Array, offset: number = 0, length: number = data.length): string {
  if (length <= ShortAsciiLength) {
    const ascii = shortAscii(data, offset, length);
    if (ascii !== void 0) return ascii;
  }
  if (utf8ValueDecoder) return utf8ValueDecoder.decode(data.subarray(offset, offset + length));
  return _utf8Read(data, offset, length);
}

/** Decode UTF8 data, potentially exceeding MAX_STRING_LENGTH. Return as primitive `string` if possible; or as `ChunkedBigString` if the result is longer than MAX_STRING_LENGTH. */
export function utf8ReadLong(data: ReadonlyUint8Array, offset: number = 0, length: number = data.length): StringLike {
  if (length <= MAX_STRING_LENGTH) {
    return utf8Read(data, offset, length);
  }
  const out = ChunkedBigString.fromUtf8Data(data, offset, offset + length);
  return out.length <= MAX_STRING_LENGTH ? out.toString() : out;
}

const utf8Encoder = typeof TextEncoder !== 'undefined' ? new TextEncoder() : undefined;
/** Chars per `TextEncoder.encode` call past an ASCII prefix: a piece is a string, and a `ChunkedBigString` is longer */
const EncodePieceLength = 1 << 24;

/**
 * The UTF-8 bytes of `str` in one array, also of a `StringLike` longer than a string can be. ASCII, the usual case of
 * text formats, is one `encodeInto` pass into an array of the exact size.
 */
export function utf8Encode(str: StringLike): Uint8Array {
  const length = str.length;
  const pieces: Uint8Array[] = [];
  let total = 0;
  let start = 0;
  if (typeof str === 'string' && utf8Encoder) {
    const out = new Uint8Array(length);
    const { read, written } = utf8Encoder.encodeInto(str, out);
    if (read === length) return out;
    pieces.push(out.subarray(0, written));
    total = written;
    start = read;
  }
  while (start < length) {
    let end = Math.min(start + EncodePieceLength, length);
    // not between the halves of a surrogate pair
    const last = str.charCodeAt(end - 1);
    if (end < length && last >= 0xd800 && last <= 0xdbff) end++;
    const piece = str.substring(start, end);
    let bytes: Uint8Array;
    if (utf8Encoder) {
      bytes = utf8Encoder.encode(piece);
    } else {
      bytes = new Uint8Array(utf8ByteCount(piece));
      utf8Write(bytes, 0, piece);
    }
    pieces.push(bytes);
    total += bytes.length;
    start = end;
  }
  if (pieces.length === 1) return pieces[0];
  const out = new Uint8Array(total);
  let offset = 0;
  for (const p of pieces) {
    out.set(p, offset);
    offset += p.length;
  }
  return out;
}

export function utf8ByteCount(str: string) {
  let count = 0;
  for (let i = 0, l = str.length; i < l; i++) {
    // A surrogate pair is one code point
    const codePoint = str.codePointAt(i)!;
    if (codePoint > 0xffff) i++;
    if (codePoint < 0x80) {
      count += 1;
      continue;
    }
    if (codePoint < 0x800) {
      count += 2;
      continue;
    }
    if (codePoint < 0x10000) {
      count += 3;
      continue;
    }
    if (codePoint < 0x110000) {
      count += 4;
      continue;
    }
    throwError('bad codepoint ' + codePoint);
  }
  return count;
}

/**
 * The index of ASCII `needle` in UTF-8 `data` at or after `from`, or -1: `String.indexOf` on the text the bytes decode
 * to, without decoding them (an ASCII byte never occurs inside a multi-byte sequence). The index is a byte index.
 *
 * The native `indexOf` finds the candidates for the needle's first lowercase letter (its first byte when it has
 * none): in text formats lowercase is rare outside names and keys (in a 124 MB mmCIF, 327 `a` against 2.2M `_`), so
 * the search runs at about the speed of a single native scan.
 */
export function asciiIndexOf(data: ReadonlyUint8Array, needle: string, from = 0): number {
  const n = needle.length;
  if (n === 0) return Math.min(from, data.length);
  let k = 0;
  for (let j = 0; j < n; j++) {
    const c = needle.charCodeAt(j);
    if (c >= 97 && c <= 122) {
      k = j;
      break;
    }
  }
  const anchor = needle.charCodeAt(k);
  const last = data.length - n;
  for (let i = data.indexOf(anchor, from + k); i >= 0 && i - k <= last; i = data.indexOf(anchor, i + 1)) {
    const start = i - k;
    let j = 0;
    while (j < n && data[start + j] === needle.charCodeAt(j)) j++;
    if (j === n) return start;
  }
  return -1;
}

/** Whether UTF-8 `data` starts with ASCII `needle` (at byte `position`). */
export function asciiStartsWith(data: ReadonlyUint8Array, needle: string, position = 0): boolean {
  if (position + needle.length > data.length) return false;
  for (let j = 0; j < needle.length; j++) if (data[position + j] !== needle.charCodeAt(j)) return false;
  return true;
}
