/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function asciiWrite(data: Uint8Array, str: string) {
    for (let i = 0, il = str.length; i < il; ++i) {
        data[i] = str.charCodeAt(i);
    }
}

/**
 * Decode a region of a Uint8Array as an ASCII string.
 */
export function asciiSlice(buf: Uint8Array, start: number, end: number): string {
    let s = '';
    for (let i = start; i < end; i++) s += String.fromCharCode(buf[i]);
    return s;
}
