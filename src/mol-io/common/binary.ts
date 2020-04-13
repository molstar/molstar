/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export const IsNativeEndianLittle = new Uint16Array(new Uint8Array([0x12, 0x34]).buffer)[0] === 0x3412;

export function flipByteOrder(data: Uint8Array, bytes: number) {
    let buffer = new ArrayBuffer(data.length);
    let ret = new Uint8Array(buffer);
    for (let i = 0, n = data.length; i < n; i += bytes) {
        for (let j = 0; j < bytes; j++) {
            ret[i + bytes - j - 1] = data[i + j];
        }
    }
    return buffer;
}

const ChunkSize = 0x7000;
export function uint8ToString(array: Uint8Array) {
    if (array.length > ChunkSize) {
        const c = [];
        for (let i = 0; i < array.length; i += ChunkSize) {
            c.push(String.fromCharCode.apply(null, array.subarray(i, i + ChunkSize)));
        }
        return c.join('');
    } else {
        return String.fromCharCode.apply(null, array);
    }
}