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

// --- Typed-array decoders ---

export function decodeFloat32(raw: Uint8Array, count: number): Float32Array {
    if (raw.byteOffset % 4 === 0) return new Float32Array(raw.buffer, raw.byteOffset, count);
    const aligned = raw.slice(0, count * 4);
    return new Float32Array(aligned.buffer, 0, count);
}

export function decodeInt64AsInt32(raw: Uint8Array, count: number): Int32Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Int32Array(count);
    for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 8, true);
    return out;
}

export function decodeToFloat64(raw: Uint8Array, type: string, count: number): Float64Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Float64Array(count);
    switch (type) {
        case 'Float32': { const f32 = decodeFloat32(raw, count); for (let i = 0; i < count; i++) out[i] = f32[i]; break; }
        case 'Float64': for (let i = 0; i < count; i++) out[i] = dv.getFloat64(i * 8, true); break;
        case 'Int8': for (let i = 0; i < count; i++) out[i] = dv.getInt8(i); break;
        case 'UInt8': for (let i = 0; i < count; i++) out[i] = dv.getUint8(i); break;
        case 'Int16': for (let i = 0; i < count; i++) out[i] = dv.getInt16(i * 2, true); break;
        case 'UInt16': for (let i = 0; i < count; i++) out[i] = dv.getUint16(i * 2, true); break;
        case 'Int32': for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 4, true); break;
        case 'UInt32': for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 4, true); break;
        case 'Int64': case 'UInt64':
            for (let i = 0; i < count; i++) {
                const lo = dv.getUint32(i * 8, true);
                const hi = dv.getInt32(i * 8 + 4, true);
                out[i] = hi * 0x100000000 + lo;
            }
            break;
        default: throw new Error(`Unsupported scalar type: "${type}".`);
    }
    return out;
}

export function decodePositions(raw: Uint8Array, type: string, count: number): Float32Array {
    if (type === 'Float32') return decodeFloat32(raw, count);
    return new Float32Array(decodeToFloat64(raw, type, count));
}

export function decodeConnectivity(raw: Uint8Array, type: string): Int32Array {
    if (type === 'Int64' || type === 'UInt64') {
        return decodeInt64AsInt32(raw, raw.byteLength / 8);
    }
    if (type === 'Int32' || type === 'UInt32') {
        if (raw.byteOffset % 4 === 0) return new Int32Array(raw.buffer, raw.byteOffset, raw.byteLength / 4);
        const aligned = raw.slice();
        return new Int32Array(aligned.buffer, 0, aligned.byteLength / 4);
    }
    throw new Error(`Unsupported connectivity type: "${type}". Expected Int32, UInt32, Int64, or UInt64.`);
}
