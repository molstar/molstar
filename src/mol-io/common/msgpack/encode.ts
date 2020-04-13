/*
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8ByteCount, utf8Write } from '../utf8';

export default function encode(value: any) {
    const buffer = new ArrayBuffer(encodedSize(value));
    const view = new DataView(buffer);
    const bytes = new Uint8Array(buffer);
    encodeInternal(value, view, bytes, 0);
    return bytes;
}

function encodedSize(value: any) {
    const type = typeof value;

    // Raw Bytes
    if (type === 'string') {
        let length = utf8ByteCount(value);
        if (length < 0x20) {
            return 1 + length;
        }
        if (length < 0x100) {
            return 2 + length;
        }
        if (length < 0x10000) {
            return 3 + length;
        }
        if (length < 0x100000000) {
            return 5 + length;
        }
    }

    if (value instanceof Uint8Array) {
        let length = value.byteLength;
        if (length < 0x100) {
            return 2 + length;
        }
        if (length < 0x10000) {
            return 3 + length;
        }
        if (length < 0x100000000) {
            return 5 + length;
        }
    }

    if (type === 'number') {
        // Floating Point
        // double
        if (Math.floor(value) !== value) return 9;

        // Integers
        if (value >= 0) {
            // positive fixnum
            if (value < 0x80) return 1;
            // uint 8
            if (value < 0x100) return 2;
            // uint 16
            if (value < 0x10000) return 3;
            // uint 32
            if (value < 0x100000000) return 5;
            throw new Error('Number too big 0x' + value.toString(16));
        }
        // negative fixnum
        if (value >= -0x20) return 1;
        // int 8
        if (value >= -0x80) return 2;
        // int 16
        if (value >= -0x8000) return 3;
        // int 32
        if (value >= -0x80000000) return 5;
        throw new Error('Number too small -0x' + value.toString(16).substr(1));
    }

    // Boolean, null
    if (type === 'boolean' || value === null || value === void 0) return 1;

    // Container Types
    if (type === 'object') {
        let length: number, size = 0;
        if (Array.isArray(value)) {
            length = value.length;
            for (let i = 0; i < length; i++) {
                size += encodedSize(value[i]);
            }
        } else {
            let keys = Object.keys(value);
            length = keys.length;
            for (let i = 0; i < length; i++) {
                let key = keys[i];
                size += encodedSize(key) + encodedSize(value[key]);
            }
        }
        if (length < 0x10) {
            return 1 + size;
        }
        if (length < 0x10000) {
            return 3 + size;
        }
        if (length < 0x100000000) {
            return 5 + size;
        }
        throw new Error('Array or object too long 0x' + length.toString(16));
    }
    throw new Error('Unknown type ' + type);
}

function encodeInternal(value: any, view: DataView, bytes: Uint8Array, offset: number) {
    let type = typeof value;

    // Strings Bytes
    if (type === 'string') {
        let length = utf8ByteCount(value);
        // fix str
        if (length < 0x20) {
            view.setUint8(offset, length | 0xa0);
            utf8Write(bytes, offset + 1, value);
            return 1 + length;
        }
        // str 8
        if (length < 0x100) {
            view.setUint8(offset, 0xd9);
            view.setUint8(offset + 1, length);
            utf8Write(bytes, offset + 2, value);
            return 2 + length;
        }
        // str 16
        if (length < 0x10000) {
            view.setUint8(offset, 0xda);
            view.setUint16(offset + 1, length);
            utf8Write(bytes, offset + 3, value);
            return 3 + length;
        }
        // str 32
        if (length < 0x100000000) {
            view.setUint8(offset, 0xdb);
            view.setUint32(offset + 1, length);
            utf8Write(bytes, offset + 5, value);
            return 5 + length;
        }
    }

    if (value instanceof Uint8Array) {
        let length = value.byteLength;
        let bytes = new Uint8Array(view.buffer);
        // bin 8
        if (length < 0x100) {
            view.setUint8(offset, 0xc4);
            view.setUint8(offset + 1, length);
            bytes.set(value, offset + 2);
            return 2 + length;
        }
        // bin 16
        if (length < 0x10000) {
            view.setUint8(offset, 0xc5);
            view.setUint16(offset + 1, length);
            bytes.set(value, offset + 3);
            return 3 + length;
        }
        // bin 32
        if (length < 0x100000000) {
            view.setUint8(offset, 0xc6);
            view.setUint32(offset + 1, length);
            bytes.set(value, offset + 5);
            return 5 + length;
        }
    }

    if (type === 'number') {
        if (!isFinite(value)) {
            throw new Error('Number not finite: ' + value);
        }

        // Floating point
        if (Math.floor(value) !== value) {
            view.setUint8(offset, 0xcb);
            view.setFloat64(offset + 1, value);
            return 9;
        }

        // Integers
        if (value >= 0) {
            // positive fixnum
            if (value < 0x80) {
                view.setUint8(offset, value);
                return 1;
            }
            // uint 8
            if (value < 0x100) {
                view.setUint8(offset, 0xcc);
                view.setUint8(offset + 1, value);
                return 2;
            }
            // uint 16
            if (value < 0x10000) {
                view.setUint8(offset, 0xcd);
                view.setUint16(offset + 1, value);
                return 3;
            }
            // uint 32
            if (value < 0x100000000) {
                view.setUint8(offset, 0xce);
                view.setUint32(offset + 1, value);
                return 5;
            }
            throw new Error('Number too big 0x' + value.toString(16));
        }
        // negative fixnum
        if (value >= -0x20) {
            view.setInt8(offset, value);
            return 1;
        }
        // int 8
        if (value >= -0x80) {
            view.setUint8(offset, 0xd0);
            view.setInt8(offset + 1, value);
            return 2;
        }
        // int 16
        if (value >= -0x8000) {
            view.setUint8(offset, 0xd1);
            view.setInt16(offset + 1, value);
            return 3;
        }
        // int 32
        if (value >= -0x80000000) {
            view.setUint8(offset, 0xd2);
            view.setInt32(offset + 1, value);
            return 5;
        }
        throw new Error('Number too small -0x' + (-value).toString(16).substr(1));
    }

    // null
    if (value === null || value === undefined) {
        view.setUint8(offset, 0xc0);
        return 1;
    }

    // Boolean
    if (type === 'boolean') {
        view.setUint8(offset, value ? 0xc3 : 0xc2);
        return 1;
    }

    // Container Types
    if (type === 'object') {
        let length: number, size = 0;
        let isArray = Array.isArray(value);
        let keys: string[] | undefined;

        if (isArray) {
            length = value.length;
        } else {
            keys = Object.keys(value);
            length = keys.length;
        }

        if (length < 0x10) {
            view.setUint8(offset, length | (isArray ? 0x90 : 0x80));
            size = 1;
        } else if (length < 0x10000) {
            view.setUint8(offset, isArray ? 0xdc : 0xde);
            view.setUint16(offset + 1, length);
            size = 3;
        } else if (length < 0x100000000) {
            view.setUint8(offset, isArray ? 0xdd : 0xdf);
            view.setUint32(offset + 1, length);
            size = 5;
        }

        if (isArray) {
            for (let i = 0; i < length; i++) {
                size += encodeInternal(value[i], view, bytes, offset + size);
            }
        } else {
            for (let i = 0, _i = keys!.length; i < _i; i++) {
                const key = keys![i];
                size += encodeInternal(key, view, bytes, offset + size);
                size += encodeInternal(value[key], view, bytes, offset + size);
            }
        }

        return size;
    }
    throw new Error('Unknown type ' + type);
}