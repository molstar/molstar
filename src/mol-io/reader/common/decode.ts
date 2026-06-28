/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { StringLike } from '../../common/string-like';

// --- Base64 ---

export const B64_DECODE = new Int8Array(256).fill(-1);
'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'.split('').forEach((c, i) => (B64_DECODE[c.charCodeAt(0)] = i));

export function decodeB64Bytes(src: Uint8Array, srcStart: number, nOut: number): Uint8Array {
    if (nOut === 0) return new Uint8Array(0);
    const out = new Uint8Array(nOut);
    let si = srcStart, di = 0;
    while (di < nOut) {
        const a = B64_DECODE[src[si++]], b = B64_DECODE[src[si++]];
        const c = B64_DECODE[src[si++]], d = B64_DECODE[src[si++]];
        if (a < 0 || b < 0) break;
        if (di < nOut) out[di++] = (a << 2) | (b >> 4);
        if (c >= 0 && di < nOut) out[di++] = ((b & 0xF) << 4) | (c >> 2);
        if (d >= 0 && di < nOut) out[di++] = ((c & 0x3) << 6) | d;
    }
    return out;
}

export function decodeB64Str(src: StringLike, srcStart: number, nOut: number): Uint8Array {
    if (nOut === 0) return new Uint8Array(0);
    const out = new Uint8Array(nOut);
    let si = srcStart, di = 0;
    while (di < nOut) {
        const a = B64_DECODE[src.charCodeAt(si++)], b = B64_DECODE[src.charCodeAt(si++)];
        const c = B64_DECODE[src.charCodeAt(si++)], d = B64_DECODE[src.charCodeAt(si++)];
        if (a < 0 || b < 0) break;
        if (di < nOut) out[di++] = (a << 2) | (b >> 4);
        if (c >= 0 && di < nOut) out[di++] = ((b & 0xF) << 4) | (c >> 2);
        if (d >= 0 && di < nOut) out[di++] = ((c & 0x3) << 6) | d;
    }
    return out;
}
