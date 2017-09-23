/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function utf8Write(data: Uint8Array, offset: number, str: string) {
    for (let i = 0, l = str.length; i < l; i++) {
        let codePoint = str.charCodeAt(i);

        // One byte of UTF-8
        if (codePoint < 0x80) {
            data[offset++] = codePoint >>> 0 & 0x7f | 0x00;
            continue;
        }

        // Two bytes of UTF-8
        if (codePoint < 0x800) {
            data[offset++] = codePoint >>> 6 & 0x1f | 0xc0;
            data[offset++] = codePoint >>> 0 & 0x3f | 0x80;
            continue;
        }

        // Three bytes of UTF-8.
        if (codePoint < 0x10000) {
            data[offset++] = codePoint >>> 12 & 0x0f | 0xe0;
            data[offset++] = codePoint >>> 6 & 0x3f | 0x80;
            data[offset++] = codePoint >>> 0 & 0x3f | 0x80;
            continue;
        }

        // Four bytes of UTF-8
        if (codePoint < 0x110000) {
            data[offset++] = codePoint >>> 18 & 0x07 | 0xf0;
            data[offset++] = codePoint >>> 12 & 0x3f | 0x80;
            data[offset++] = codePoint >>> 6 & 0x3f | 0x80;
            data[offset++] = codePoint >>> 0 & 0x3f | 0x80;
            continue;
        }
        throw new Error('bad codepoint ' + codePoint);
    }
}

const __chars = function () {
    let data: string[] = [];
    for (let i = 0; i < 1024; i++) data[i] = String.fromCharCode(i);
    return data;
}();

function throwError(err: string) {
    throw new Error(err);
}

export function utf8Read(data: Uint8Array, offset: number, length: number) {
    let chars = __chars;
    let str: string[] | undefined = void 0, chunk: string[] = [], chunkSize = 512, chunkOffset = 0;

    for (let i = offset, end = offset + length; i < end; i++) {
        let byte = data[i];
        // One byte character
        if ((byte & 0x80) === 0x00) {
            chunk[chunkOffset++] = chars[byte];
        }
        // Two byte character
        else if ((byte & 0xe0) === 0xc0) {
            chunk[chunkOffset++] = chars[((byte & 0x0f) << 6) | (data[++i] & 0x3f)];
        }
        // Three byte character
        else if ((byte & 0xf0) === 0xe0) {
            chunk[chunkOffset++] = String.fromCharCode(
                ((byte & 0x0f) << 12) |
                ((data[++i] & 0x3f) << 6) |
                ((data[++i] & 0x3f) << 0)
            );
        }
        // Four byte character
        else if ((byte & 0xf8) === 0xf0) {
            chunk[chunkOffset++] = String.fromCharCode(
                ((byte & 0x07) << 18) |
                ((data[++i] & 0x3f) << 12) |
                ((data[++i] & 0x3f) << 6) |
                ((data[++i] & 0x3f) << 0)
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

export function utf8ByteCount(str: string) {
    let count = 0;
    for (let i = 0, l = str.length; i < l; i++) {
        let codePoint = str.charCodeAt(i);
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