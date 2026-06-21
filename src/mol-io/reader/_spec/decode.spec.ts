/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import {
    B64_DECODE,
    decodeB64Bytes, decodeB64Str,
    decodeFloat32, decodeInt64AsInt32, decodeToFloat64,
    decodePositions, decodeConnectivity,
} from '../common/decode';

// Known base64: "Man" → [0x4D, 0x61, 0x6E]
const MAN_B64 = 'TWFu';
const MAN_BYTES = new Uint8Array([0x4D, 0x61, 0x6E]);

describe('B64_DECODE', () => {
    it('has -1 for non-base64 chars', () => {
        expect(B64_DECODE['!'.charCodeAt(0)]).toBe(-1);
        expect(B64_DECODE[' '.charCodeAt(0)]).toBe(-1);
    });
    it('maps A=0, Z=25, a=26, z=51, 0=52, 9=61, +=62, /=63', () => {
        expect(B64_DECODE['A'.charCodeAt(0)]).toBe(0);
        expect(B64_DECODE['Z'.charCodeAt(0)]).toBe(25);
        expect(B64_DECODE['a'.charCodeAt(0)]).toBe(26);
        expect(B64_DECODE['z'.charCodeAt(0)]).toBe(51);
        expect(B64_DECODE['0'.charCodeAt(0)]).toBe(52);
        expect(B64_DECODE['9'.charCodeAt(0)]).toBe(61);
        expect(B64_DECODE['+'.charCodeAt(0)]).toBe(62);
        expect(B64_DECODE['/'.charCodeAt(0)]).toBe(63);
    });
});

describe('decodeB64Str', () => {
    it('decodes "Man"', () => {
        const out = decodeB64Str(MAN_B64, 0, 3);
        expect(Array.from(out)).toEqual(Array.from(MAN_BYTES));
    });
    it('returns empty for nOut=0', () => {
        expect(decodeB64Str(MAN_B64, 0, 0).length).toBe(0);
    });
    it('respects srcStart offset', () => {
        const padded = 'AAAA' + MAN_B64;
        const out = decodeB64Str(padded, 4, 3);
        expect(Array.from(out)).toEqual(Array.from(MAN_BYTES));
    });
});

describe('decodeB64Bytes', () => {
    it('decodes "Man" from Uint8Array', () => {
        const src = new Uint8Array(MAN_B64.split('').map(c => c.charCodeAt(0)));
        const out = decodeB64Bytes(src, 0, 3);
        expect(Array.from(out)).toEqual(Array.from(MAN_BYTES));
    });
    it('returns empty for nOut=0', () => {
        const src = new Uint8Array(4);
        expect(decodeB64Bytes(src, 0, 0).length).toBe(0);
    });
    it('respects srcStart offset', () => {
        const prefix = new Uint8Array(4).fill('A'.charCodeAt(0));
        const rest = new Uint8Array(MAN_B64.split('').map(c => c.charCodeAt(0)));
        const src = new Uint8Array([...prefix, ...rest]);
        const out = decodeB64Bytes(src, 4, 3);
        expect(Array.from(out)).toEqual(Array.from(MAN_BYTES));
    });
});

describe('decodeFloat32', () => {
    it('reuses buffer when byteOffset is aligned', () => {
        const raw = new Uint8Array(new Float32Array([1.5, 2.5, -3.0]).buffer);
        const out = decodeFloat32(raw, 3);
        expect(out[0]).toBeCloseTo(1.5);
        expect(out[1]).toBeCloseTo(2.5);
        expect(out[2]).toBeCloseTo(-3.0);
        expect(out.buffer).toBe(raw.buffer);
    });
    it('copies when byteOffset is unaligned', () => {
        const base = new Uint8Array(new Float32Array([7.0, 8.0]).buffer);
        const unaligned = new Uint8Array(base.buffer, 0, 8); // aligned here actually
        // Create truly unaligned view via offset=1 in a larger buffer
        const big = new Uint8Array(9);
        big.set(base, 1);
        const unalignedView = new Uint8Array(big.buffer, 1, 8);
        const out = decodeFloat32(unalignedView, 2);
        expect(out[0]).toBeCloseTo(7.0);
        expect(out[1]).toBeCloseTo(8.0);
    });
});

describe('decodeInt64AsInt32', () => {
    it('takes low 32 bits of each int64 (little-endian)', () => {
        const buf = new ArrayBuffer(16);
        const dv = new DataView(buf);
        dv.setUint32(0, 42, true); dv.setUint32(4, 0, true);
        dv.setUint32(8, 999, true); dv.setUint32(12, 0, true);
        const out = decodeInt64AsInt32(new Uint8Array(buf), 2);
        expect(out[0]).toBe(42);
        expect(out[1]).toBe(999);
    });
});

describe('decodeToFloat64', () => {
    it('Float32', () => {
        const raw = new Uint8Array(new Float32Array([1.0, -2.0]).buffer);
        const out = decodeToFloat64(raw, 'Float32', 2);
        expect(out[0]).toBeCloseTo(1.0);
        expect(out[1]).toBeCloseTo(-2.0);
    });
    it('Float64', () => {
        const raw = new Uint8Array(new Float64Array([Math.PI, Math.E]).buffer);
        const out = decodeToFloat64(raw, 'Float64', 2);
        expect(out[0]).toBeCloseTo(Math.PI, 10);
        expect(out[1]).toBeCloseTo(Math.E, 10);
    });
    it('Int8', () => {
        const raw = new Uint8Array(new Int8Array([-1, 127]).buffer);
        const out = decodeToFloat64(raw, 'Int8', 2);
        expect(out[0]).toBe(-1);
        expect(out[1]).toBe(127);
    });
    it('UInt8', () => {
        const raw = new Uint8Array([0, 255]);
        const out = decodeToFloat64(raw, 'UInt8', 2);
        expect(out[0]).toBe(0);
        expect(out[1]).toBe(255);
    });
    it('Int16', () => {
        const raw = new Uint8Array(new Int16Array([-32768, 32767]).buffer);
        const out = decodeToFloat64(raw, 'Int16', 2);
        expect(out[0]).toBe(-32768);
        expect(out[1]).toBe(32767);
    });
    it('UInt16', () => {
        const raw = new Uint8Array(new Uint16Array([0, 65535]).buffer);
        const out = decodeToFloat64(raw, 'UInt16', 2);
        expect(out[0]).toBe(0);
        expect(out[1]).toBe(65535);
    });
    it('Int32', () => {
        const raw = new Uint8Array(new Int32Array([-1, 100]).buffer);
        const out = decodeToFloat64(raw, 'Int32', 2);
        expect(out[0]).toBe(-1);
        expect(out[1]).toBe(100);
    });
    it('UInt32', () => {
        const raw = new Uint8Array(new Uint32Array([0, 4294967295]).buffer);
        const out = decodeToFloat64(raw, 'UInt32', 2);
        expect(out[0]).toBe(0);
        expect(out[1]).toBe(4294967295);
    });
    it('Int64 small values', () => {
        const buf = new ArrayBuffer(16);
        const dv = new DataView(buf);
        dv.setUint32(0, 7, true); dv.setInt32(4, 0, true);
        dv.setUint32(8, 0, true); dv.setInt32(12, -1, true); // -4294967296
        const out = decodeToFloat64(new Uint8Array(buf), 'Int64', 2);
        expect(out[0]).toBe(7);
        expect(out[1]).toBe(-4294967296);
    });
    it('throws for unknown type', () => {
        expect(() => decodeToFloat64(new Uint8Array(4), 'Bogus', 1)).toThrow();
    });
});

describe('decodePositions', () => {
    it('returns Float32 directly for Float32 type', () => {
        const raw = new Uint8Array(new Float32Array([0.5, 1.5, 2.5]).buffer);
        const out = decodePositions(raw, 'Float32', 3);
        expect(out).toBeInstanceOf(Float32Array);
        expect(out[1]).toBeCloseTo(1.5);
    });
    it('converts Float64 to Float32', () => {
        const raw = new Uint8Array(new Float64Array([1.0, 2.0]).buffer);
        const out = decodePositions(raw, 'Float64', 2);
        expect(out).toBeInstanceOf(Float32Array);
        expect(out[0]).toBeCloseTo(1.0);
    });
});

describe('decodeConnectivity', () => {
    it('Int32 aligned', () => {
        const raw = new Uint8Array(new Int32Array([0, 1, 2]).buffer);
        const out = decodeConnectivity(raw, 'Int32');
        expect(Array.from(out)).toEqual([0, 1, 2]);
    });
    it('UInt32', () => {
        const raw = new Uint8Array(new Uint32Array([3, 4, 5]).buffer);
        const out = decodeConnectivity(raw, 'UInt32');
        expect(Array.from(out)).toEqual([3, 4, 5]);
    });
    it('Int64 takes low 32 bits', () => {
        const buf = new ArrayBuffer(24);
        const dv = new DataView(buf);
        dv.setUint32(0, 10, true); dv.setUint32(4, 0, true);
        dv.setUint32(8, 20, true); dv.setUint32(12, 0, true);
        dv.setUint32(16, 30, true); dv.setUint32(20, 0, true);
        const out = decodeConnectivity(new Uint8Array(buf), 'Int64');
        expect(Array.from(out)).toEqual([10, 20, 30]);
    });
    it('throws for unsupported type', () => {
        expect(() => decodeConnectivity(new Uint8Array(4), 'Float32')).toThrow();
    });
});
