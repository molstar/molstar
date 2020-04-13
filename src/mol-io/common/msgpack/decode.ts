/*
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8Read } from '../utf8';

export default function decode(buffer: Uint8Array) {
    return parse({ buffer, offset: 0, dataView: new DataView(buffer.buffer) });
}

// Loosely based on
// The MIT License (MIT)
// Copyright (c) 2013 Tim Caswell <tim@creationix.com>
// https://github.com/creationix/msgpack-js

interface State {
    buffer: Uint8Array,
    offset: number,
    dataView: DataView
}

/**
 * decode all key-value pairs of a map into an object
 */
function map(state: State, length: number) {
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
function bin(state: State, length: number) {
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
function str(state: State, length: number) {
    const value = utf8Read(state.buffer, state.offset, length);
    state.offset += length;
    return value;
}

/**
 * decode array
 */
function array(state: State, length: number) {
    const value: any[] = new Array(length);
    for (let i = 0; i < length; i++) {
        value[i] = parse(state);
    }
    return value;
}

/**
 * recursively parse the MessagePack data and return decoded MessagePack data
 */
function parse(state: State) {
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