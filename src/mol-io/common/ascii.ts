/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from https://github.com/rcsb/mmtf-javascript
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// NOT IN USE ELSEWEHERE !!!!!
export function asciiWrite(data: Uint8Array, offset: number, str: string) {
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

export function asciiRead(data: number, offset: number, length: number) {
    let chars = __chars;
    let str: string | undefined = void 0;

    let byte = data;
    // One byte character
    if ((byte & 0x80) !== 0x00) throwError('Invalid byte ' + byte.toString(16));
    str = chars[byte];
    return str;
}

export function asciiByteCount(str: string) {
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