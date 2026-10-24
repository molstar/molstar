/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8Read } from '@molstar/core/util/utf8';

/**
 * Decode an extension value from its type code and payload (a view into the buffer: copy what is kept).
 * Returning undefined gives `{ type, data }`.
 */
export type MsgPackExtDecoder = (type: number, data: Uint8Array) => unknown;

export function decodeMsgPack(buffer: Uint8Array, ext?: MsgPackExtDecoder) {
  const state = createMsgPackState(buffer);
  state.ext = ext;
  return parse(state);
}

// Loosely based on
// The MIT License (MIT)
// Copyright (c) 2013 Tim Caswell <tim@creationix.com>
// https://github.com/creationix/msgpack-js

export interface MsgPackState {
  buffer: Uint8Array;
  offset: number;
  dataView: DataView;
  ext?: MsgPackExtDecoder;
}

export function createMsgPackState(buffer: Uint8Array): MsgPackState {
  return { buffer, offset: 0, dataView: new DataView(buffer.buffer, buffer.byteOffset, buffer.byteLength) };
}

/** Decode the value at the current offset and advance past it. */
export { parse as parseMsgPackValue };

/** Advance past a map header and return the number of key-value pairs it contains. */
export function readMsgPackMapLength(state: MsgPackState) {
  const type = state.buffer[state.offset];
  // FixMap
  if ((type & 0xf0) === 0x80) {
    state.offset++;
    return type & 0x0f;
  }
  // map 16
  if (type === 0xde) {
    const length = state.dataView.getUint16(state.offset + 1);
    state.offset += 3;
    return length;
  }
  // map 32
  if (type === 0xdf) {
    const length = state.dataView.getUint32(state.offset + 1);
    state.offset += 5;
    return length;
  }
  throw new Error('Expected a map, got type 0x' + type.toString(16));
}

/** Advance past an array header and return the number of elements it contains. */
export function readMsgPackArrayLength(state: MsgPackState) {
  const type = state.buffer[state.offset];
  // FixArray
  if ((type & 0xf0) === 0x90) {
    state.offset++;
    return type & 0x0f;
  }
  // array 16
  if (type === 0xdc) {
    const length = state.dataView.getUint16(state.offset + 1);
    state.offset += 3;
    return length;
  }
  // array 32
  if (type === 0xdd) {
    const length = state.dataView.getUint32(state.offset + 1);
    state.offset += 5;
    return length;
  }
  throw new Error('Expected an array, got type 0x' + type.toString(16));
}

/**
 * A length header is trusted only as far as the bytes left can hold it: `length` items of at least `bytes` each.
 * A header without its payload would otherwise allocate gigabytes or loop for minutes.
 */
function room(state: MsgPackState, length: number, bytes: number) {
  if (length * bytes > state.buffer.length - state.offset) {
    throw new Error(
      `truncated: ${length} item(s) at byte ${state.offset}, ${state.buffer.length - state.offset} byte(s) left`,
    );
  }
}

/**
 * decode all key-value pairs of a map into an object
 */
function map(state: MsgPackState, length: number) {
  room(state, length, 2);
  const value: { [k: string]: any } = {};
  for (let i = 0; i < length; i++) {
    const key = parse(state);
    // `__proto__` is an own key, as JSON.parse makes it, never the object's prototype
    if (key === '__proto__') {
      Object.defineProperty(value, key, { value: parse(state), enumerable: true, writable: true, configurable: true });
    } else {
      value[key] = parse(state);
    }
  }
  return value;
}

/**
 * decode binary array
 */
function bin(state: MsgPackState, length: number) {
  room(state, length, 1);
  // a copy, not a view: element access on the copy is faster, and its buffer is exactly the bytes, which
  // `view(data.buffer)` callers rely on. `set` copies natively (the byte loop was 10× slower, perf-audit-3 13/P1);
  // not `slice`, which makes a `Buffer` from a `Buffer` input (species), whose own `slice` is a view.
  const o = state.offset;
  const value = new Uint8Array(length);
  value.set(state.buffer.subarray(o, o + length));
  state.offset += length;
  return value;
}

/**
 * decode an extension value, the offset at its type code
 */
function extension(state: MsgPackState, length: number) {
  room(state, 1 + length, 1);
  const type = state.dataView.getInt8(state.offset);
  const data = state.buffer.subarray(state.offset + 1, state.offset + 1 + length);
  state.offset += 1 + length;
  const value = state.ext?.(type, data);
  return value !== undefined ? value : { type, data: data.slice() };
}

/**
 * decode string
 */
function str(state: MsgPackState, length: number) {
  room(state, length, 1);
  const value = utf8Read(state.buffer, state.offset, length);
  state.offset += length;
  return value;
}

/**
 * decode array
 */
function array(state: MsgPackState, length: number) {
  room(state, length, 1);
  // Grown as items are read, never sized by the header: nested headers may each claim every byte left,
  // and sizing each before reading it would multiply the claim by the depth.
  const value: any[] = [];
  for (let i = 0; i < length; i++) {
    value.push(parse(state));
  }
  return value;
}

/**
 * recursively parse the MessagePack data and return decoded MessagePack data
 */
function parse(state: MsgPackState) {
  // Past the end the byte is undefined, which would read as a positive fixint
  if (state.offset >= state.buffer.length) throw new Error(`truncated: no value at byte ${state.offset}`);
  const type = state.buffer[state.offset];
  let value: any, length: number;
  // Positive FixInt
  if ((type & 0x80) === 0x00) {
    state.offset++;
    return type;
  }
  // FixMap
  if ((type & 0xf0) === 0x80) {
    length = type & 0x0f;
    state.offset++;
    return map(state, length);
  }
  // FixArray
  if ((type & 0xf0) === 0x90) {
    length = type & 0x0f;
    state.offset++;
    return array(state, length);
  }
  // FixStr
  if ((type & 0xe0) === 0xa0) {
    length = type & 0x1f;
    state.offset++;
    return str(state, length);
  }
  // Negative FixInt
  if ((type & 0xe0) === 0xe0) {
    value = state.dataView.getInt8(state.offset);
    state.offset++;
    return value;
  }
  switch (type) {
    // nil
    case 0xc0:
      state.offset++;
      return null;
    // false
    case 0xc2:
      state.offset++;
      return false;
    // true
    case 0xc3:
      state.offset++;
      return true;
    // bin 8
    case 0xc4:
      length = state.dataView.getUint8(state.offset + 1);
      state.offset += 2;
      return bin(state, length);
    // bin 16
    case 0xc5:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return bin(state, length);
    // bin 32
    case 0xc6:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return bin(state, length);
    // float 32
    case 0xca:
      value = state.dataView.getFloat32(state.offset + 1);
      state.offset += 5;
      return value;
    // float 64
    case 0xcb:
      value = state.dataView.getFloat64(state.offset + 1);
      state.offset += 9;
      return value;
    // uint8
    case 0xcc:
      value = state.buffer[state.offset + 1];
      state.offset += 2;
      return value;
    // uint 16
    case 0xcd:
      value = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return value;
    // uint 32
    case 0xce:
      value = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return value;
    // int 8
    case 0xd0:
      value = state.dataView.getInt8(state.offset + 1);
      state.offset += 2;
      return value;
    // int 16
    case 0xd1:
      value = state.dataView.getInt16(state.offset + 1);
      state.offset += 3;
      return value;
    // int 32
    case 0xd2:
      value = state.dataView.getInt32(state.offset + 1);
      state.offset += 5;
      return value;
    // str 8
    case 0xd9:
      length = state.dataView.getUint8(state.offset + 1);
      state.offset += 2;
      return str(state, length);
    // str 16
    case 0xda:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return str(state, length);
    // str 32
    case 0xdb:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return str(state, length);
    // array 16
    case 0xdc:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return array(state, length);
    // array 32
    case 0xdd:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return array(state, length);
    // map 16:
    case 0xde:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return map(state, length);
    // map 32
    case 0xdf:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return map(state, length);
    // fixext 1, 2, 4, 8, 16
    case 0xd4:
    case 0xd5:
    case 0xd6:
    case 0xd7:
    case 0xd8:
      state.offset += 1;
      return extension(state, 1 << (type - 0xd4));
    // ext 8, 16, 32
    case 0xc7:
      length = state.dataView.getUint8(state.offset + 1);
      state.offset += 2;
      return extension(state, length);
    case 0xc8:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      return extension(state, length);
    case 0xc9:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      return extension(state, length);
  }

  throw new Error('Unknown type 0x' + type.toString(16));
}

/**
 * Advance past the value at the current offset without decoding it.
 * Mirrors the type dispatch of `parse` but allocates nothing, which makes it suitable for
 * cheaply skipping large payloads (e.g. BinaryCIF column data) when only some fields are needed.
 */
export function skipMsgPackValue(state: MsgPackState) {
  const type = state.buffer[state.offset];
  let length: number;
  // Positive FixInt, Negative FixInt
  if ((type & 0x80) === 0x00 || (type & 0xe0) === 0xe0) {
    state.offset++;
    return;
  }
  // FixMap
  if ((type & 0xf0) === 0x80) {
    length = type & 0x0f;
    state.offset++;
    skipValues(state, 2 * length);
    return;
  }
  // FixArray
  if ((type & 0xf0) === 0x90) {
    length = type & 0x0f;
    state.offset++;
    skipValues(state, length);
    return;
  }
  // FixStr
  if ((type & 0xe0) === 0xa0) {
    length = type & 0x1f;
    state.offset += 1 + length;
    return;
  }
  switch (type) {
    // nil, false, true
    case 0xc0:
    case 0xc2:
    case 0xc3:
      state.offset += 1;
      return;
    // bin 8, str 8
    case 0xc4:
    case 0xd9:
      length = state.dataView.getUint8(state.offset + 1);
      state.offset += 2 + length;
      return;
    // bin 16, str 16
    case 0xc5:
    case 0xda:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3 + length;
      return;
    // bin 32, str 32
    case 0xc6:
    case 0xdb:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5 + length;
      return;
    // uint 8, int 8
    case 0xcc:
    case 0xd0:
      state.offset += 2;
      return;
    // uint 16, int 16
    case 0xcd:
    case 0xd1:
      state.offset += 3;
      return;
    // float 32, uint 32, int 32
    case 0xca:
    case 0xce:
    case 0xd2:
      state.offset += 5;
      return;
    // float 64
    case 0xcb:
      state.offset += 9;
      return;
    // array 16
    case 0xdc:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      skipValues(state, length);
      return;
    // array 32
    case 0xdd:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      skipValues(state, length);
      return;
    // map 16
    case 0xde:
      length = state.dataView.getUint16(state.offset + 1);
      state.offset += 3;
      skipValues(state, 2 * length);
      return;
    // map 32
    case 0xdf:
      length = state.dataView.getUint32(state.offset + 1);
      state.offset += 5;
      skipValues(state, 2 * length);
      return;
    // fixext 1, 2, 4, 8, 16; ext 8, 16, 32
    case 0xd4:
    case 0xd5:
    case 0xd6:
    case 0xd7:
    case 0xd8:
      state.offset += 2 + (1 << (type - 0xd4));
      return;
    case 0xc7:
      state.offset += 3 + state.dataView.getUint8(state.offset + 1);
      return;
    case 0xc8:
      state.offset += 4 + state.dataView.getUint16(state.offset + 1);
      return;
    case 0xc9:
      state.offset += 6 + state.dataView.getUint32(state.offset + 1);
      return;
  }

  throw new Error('Unknown type 0x' + type.toString(16));
}

function skipValues(state: MsgPackState, count: number) {
  for (let i = 0; i < count; i++) skipMsgPackValue(state);
}
