/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8Read } from '../utf8';

export function decodeMsgPack(buffer: Uint8Array) {
    return parse(createMsgPackState(buffer));
}

// Loosely based on
// The MIT License (MIT)
// Copyright (c) 2013 Tim Caswell <tim@creationix.com>
// https://github.com/creationix/msgpack-js

export interface MsgPackState {
    buffer: Uint8Array,
    offset: number,
    dataView: DataView
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
 * decode all key-value pairs of a map into an object
 */
function map(state: MsgPackState, length: number) {
    const value: { [k: string]: any } = {};
    for (let i = 0; i < length; i++) {
        const key = parse(state);
        value[key] = parse(state);
    }
    return value;
}

/**
 * decode binary array
 */
function bin(state: MsgPackState, length: number) {
    // This approach to binary parsing wastes a bit of memory to trade for speed compared to:
    //
    //   let value = buffer.subarray(offset, offset + length); //new Uint8Array(buffer.buffer, offset, length);
    //
    // It turns out that using the view created by subarray probably uses DataView
    // in the background, which causes the element access to be several times slower
    // than creating the new byte array.

    const value = new Uint8Array(length);
    const o = state.offset;
    for (let i = 0; i < length; i++) value[i] = state.buffer[i + o];
    state.offset += length;
    return value;
}

/**
 * decode string
 */
function str(state: MsgPackState, length: number) {
    const value = utf8Read(state.buffer, state.offset, length);
    state.offset += length;
    return value;
}

/**
 * decode array
 */
function array(state: MsgPackState, length: number) {
    const value: any[] = new Array(length);
    for (let i = 0; i < length; i++) {
        value[i] = parse(state);
    }
    return value;
}

/**
 * recursively parse the MessagePack data and return decoded MessagePack data
 */
function parse(state: MsgPackState) {
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
    }

    throw new Error('Unknown type 0x' + type.toString(16));
}

function skipValues(state: MsgPackState, count: number) {
    for (let i = 0; i < count; i++) skipMsgPackValue(state);
}
