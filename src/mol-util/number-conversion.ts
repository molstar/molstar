/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * This code has been modified from https://github.com/mrdoob/three.js/,
 * copyright (c) 2010-2024 three.js authors. MIT License
 */

// Fast Half Float Conversions, http://www.fox-toolkit.org/ftp/fasthalffloatconversion.pdf

import { clamp } from '../mol-math/interpolate';

const Tables = generateTables();

function generateTables() {
    // float32 to float16 helpers
    const buffer = new ArrayBuffer(4);
    const floatView = new Float32Array(buffer);
    const uint32View = new Uint32Array(buffer);

    const baseTable = new Uint32Array(512);
    const shiftTable = new Uint32Array(512);

    for (let i = 0; i < 256; ++i) {
        const e = i - 127;
        if (e < -27) { // very small number (0, -0)
            baseTable[i] = 0x0000;
            baseTable[i | 0x100] = 0x8000;
            shiftTable[i] = 24;
            shiftTable[i | 0x100] = 24;
        } else if (e < -14) { // small number (denorm)
            baseTable[i] = 0x0400 >> (-e - 14);
            baseTable[i | 0x100] = (0x0400 >> (-e - 14)) | 0x8000;
            shiftTable[i] = - e - 1;
            shiftTable[i | 0x100] = -e - 1;
        } else if (e <= 15) { // normal number
            baseTable[i] = (e + 15) << 10;
            baseTable[i | 0x100] = ((e + 15) << 10) | 0x8000;
            shiftTable[i] = 13;
            shiftTable[i | 0x100] = 13;
        } else if (e < 128) { // large number (Infinity, -Infinity)
            baseTable[i] = 0x7c00;
            baseTable[i | 0x100] = 0xfc00;
            shiftTable[i] = 24;
            shiftTable[i | 0x100] = 24;
        } else { // stay (NaN, Infinity, -Infinity)
            baseTable[i] = 0x7c00;
            baseTable[i | 0x100] = 0xfc00;
            shiftTable[i] = 13;
            shiftTable[i | 0x100] = 13;
        }
    }

    // float16 to float32 helpers
    const mantissaTable = new Uint32Array(2048);
    const exponentTable = new Uint32Array(64);
    const offsetTable = new Uint32Array(64);

    for (let i = 1; i < 1024; ++i) {
        let m = i << 13; // zero pad mantissa bits
        let e = 0; // zero exponent

        // normalized
        while ((m & 0x00800000) === 0) {
            m <<= 1;
            e -= 0x00800000; // decrement exponent
        }

        m &= ~ 0x00800000; // clear leading 1 bit
        e += 0x38800000; // adjust bias

        mantissaTable[i] = m | e;
    }

    for (let i = 1024; i < 2048; ++i) {
        mantissaTable[i] = 0x38000000 + ((i - 1024) << 13);
    }

    for (let i = 1; i < 31; ++i) {
        exponentTable[i] = i << 23;
    }

    exponentTable[31] = 0x47800000;
    exponentTable[32] = 0x80000000;

    for (let i = 33; i < 63; ++i) {
        exponentTable[i] = 0x80000000 + ((i - 32) << 23);
    }

    exponentTable[63] = 0xc7800000;

    for (let i = 1; i < 64; ++i) {
        if (i !== 32) {
            offsetTable[i] = 1024;
        }
    }

    return {
        floatView,
        uint32View,
        baseTable,
        shiftTable,
        mantissaTable,
        exponentTable,
        offsetTable
    };
}

/** float32 to float16 */
export function toHalfFloat(val: number) {
    val = clamp(val, -65504, 65504);
    Tables.floatView[0] = val;
    const f = Tables.uint32View[0];
    const e = (f >> 23) & 0x1ff;
    return Tables.baseTable[e] + ((f & 0x007fffff) >> Tables.shiftTable[e]);
}

/** float16 to float32 */
export function fromHalfFloat(val: number) {
    const m = val >> 10;
    Tables.uint32View[0] = Tables.mantissaTable[Tables.offsetTable[m] + (val & 0x3ff)] + Tables.exponentTable[m];
    return Tables.floatView[0];
}
