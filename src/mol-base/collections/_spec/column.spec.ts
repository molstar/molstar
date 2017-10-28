/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column, { isTypedArray } from '../column'

describe('column', () => {
    const cc = Column.ofConst(10, 2, Column.Type.int);
    const arr = Column.ofArray({ array: [1, 2, 3, 4], type: Column.Type.int });
    const arrWindow = Column.window(arr, 1, 3);

    const typed = Column.ofArray({ array: new Int32Array([1, 2, 3, 4]), type: Column.Type.int });
    const typedWindow = Column.window(typed, 1, 3);

    const numStr = Column.ofArray({ array: [1, 2] as any, type: Column.Type.str });

    it('constant', () => {
        expect(cc.rowCount).toBe(2);
        expect(cc.value(0)).toBe(10);
    });

    it('arr', () => {
        expect(arr.rowCount).toBe(4);
        expect(arr.value(1)).toBe(2);
        expect(arrWindow.value(0)).toBe(2);
        expect(arrWindow.rowCount).toBe(2);
    });

    it('typed', () => {
        expect(typedWindow.value(0)).toBe(2);
        expect(typedWindow.rowCount).toBe(2);
        expect(isTypedArray(typedWindow.toArray())).toBe(true);
    });

    it('numStr', () => {
        expect(numStr.value(0)).toBe('1');
        expect(numStr.toArray()).toEqual(['1', '2']);
    });
})