/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from './column';

export function getArrayBounds(rowCount: number, params?: Column.ToArrayParams<any>) {
    const start = params && typeof params.start !== 'undefined' ? Math.max(Math.min(params.start, rowCount - 1), 0) : 0;
    const end = params && typeof params.end !== 'undefined' ? Math.min(params.end, rowCount) : rowCount;
    return { start, end };
}

export function createArray(rowCount: number, params?: Column.ToArrayParams<any>) {
    const c = params && typeof params.array !== 'undefined' ? params.array : Array;
    const { start, end } = getArrayBounds(rowCount, params);
    return { array: new c(end - start) as any[], start, end };
}

export function fillArrayValues(value: (row: number) => any, target: any[], start: number) {
    for (let i = 0, _e = target.length; i < _e; i++) target[i] = value(start + i);
    return target;
}

export function createAndFillArray(rowCount: number, value: (row: number) => any, params?: Column.ToArrayParams<any>) {
    const { array, start } = createArray(rowCount, params);
    return fillArrayValues(value, array, start);
}

export function isTypedArray(data: any): boolean {
    return !!data.buffer && typeof data.byteLength === 'number' && typeof data.BYTES_PER_ELEMENT === 'number';
}

export function typedArrayWindow(data: any, params?: Column.ToArrayParams<any>): ReadonlyArray<number> {
    const { constructor, buffer, length, byteOffset, BYTES_PER_ELEMENT } = data;
    const { start, end } = getArrayBounds(length, params);
    if (start === 0 && end === length) return data;
    return new constructor(buffer, byteOffset + BYTES_PER_ELEMENT * start, Math.min(length, end - start));
}