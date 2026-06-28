/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import {
    B64_DECODE,
    decodeB64Bytes, decodeB64Str,
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
