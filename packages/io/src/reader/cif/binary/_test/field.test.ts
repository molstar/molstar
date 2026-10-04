/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Aniruddha Adak <aniruddhaadak80@gmail.com>
 */

import { ArrayEncoder, ArrayEncoding as E, type EncodedColumn, classifyIntArray, classifyFloatArray } from '../../../../common/binary-cif.js';
import { Field } from '../field.js';

/** mirrors `Column.ValueKinds` */
const Present = 0, NotPresent = 1, Unknown = 2;

function intColumn(data: number[], mask?: number[]): EncodedColumn {
    const array = new Int32Array(data);
    return {
        name: 'test',
        data: classifyIntArray(array).encode(array),
        mask: mask ? ArrayEncoder.by(E.byteArray).encode(new Uint8Array(mask)) : void 0
    };
}

function floatColumn(data: number[], mask?: number[]): EncodedColumn {
    const array = new Float32Array(data);
    return {
        name: 'test',
        data: classifyFloatArray(array).encode(array),
        mask: mask ? ArrayEncoder.by(E.byteArray).encode(new Uint8Array(mask)) : void 0
    };
}

const IntValues = [1, -2, 3, -4, 5];
const FloatValues = [1.5, -2.5, 3.5, -4.5, 5.5];
/** rows 1 (`?`) and 3 (`.`) are missing, the rest is present */
const Mask = [Present, Unknown, Present, NotPresent, Present];

describe('binary-cif Field', () => {
    it('masked int values default to 0', () => {
        const f = Field(intColumn(IntValues, Mask));

        expect(f.rowCount).toBe(5);
        expect(f.int(0)).toBe(1);
        expect(f.int(1)).toBe(0);
        expect(f.int(2)).toBe(3);
        expect(f.int(3)).toBe(0);
        expect(f.int(4)).toBe(5);
    });

    it('masked float values default to 0', () => {
        const f = Field(floatColumn(FloatValues, Mask));

        expect(f.float(0)).toBe(1.5);
        expect(f.float(1)).toBe(0);
        expect(f.float(2)).toBe(3.5);
        expect(f.float(3)).toBe(0);
        expect(f.float(4)).toBe(5.5);
    });

    it('toIntArray defaults masked values to 0 and keeps the source data intact', () => {
        const f = Field(intColumn(IntValues, Mask));

        expect(Array.from(f.toIntArray())).toEqual([1, 0, 3, 0, 5]);
        expect(Array.from(f.toIntArray({ start: 1, end: 4 }))).toEqual([0, 3, 0]);

        expect(Array.from(f.__array as ArrayLike<number>)).toEqual(IntValues);
    });

    it('toFloatArray defaults masked values to 0', () => {
        const f = Field(floatColumn(FloatValues, Mask));

        expect(Array.from(f.toFloatArray())).toEqual([1.5, 0, 3.5, 0, 5.5]);
        expect(Array.from(f.toFloatArray({ start: 1, end: 4 }))).toEqual([0, 3.5, 0]);

        expect(Array.from(f.__array as ArrayLike<number>)).toEqual(FloatValues);
    });

    it('leaves str and valueKind unchanged', () => {
        const f = Field(intColumn([1, 2, 3], [Present, Unknown, NotPresent]));

        expect(f.str(0)).toBe('1');
        expect(f.str(1)).toBe('');
        expect(f.str(2)).toBe('');

        expect(f.valueKind(0)).toBe(Present);
        expect(f.valueKind(1)).toBe(Unknown);
        expect(f.valueKind(2)).toBe(NotPresent);
    });

    it('returns the encoded values for a trivial mask', () => {
        const f = Field(intColumn([1, 2, 3], [Present, Present, Present]));

        expect(f.int(1)).toBe(2);
        expect(Array.from(f.toIntArray())).toEqual([1, 2, 3]);
        expect(f.valueKind(0)).toBe(Present);
    });

    it('returns the encoded values when there is no mask', () => {
        const f = Field(intColumn([1, 2, 3]));

        expect(f.int(1)).toBe(2);
        expect(Array.from(f.toIntArray())).toEqual([1, 2, 3]);
        expect(f.valueKind(0)).toBe(Present);
    });
});
