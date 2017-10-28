/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column, { createArray } from '../../../mol-base/collections/column'

export interface File {
    readonly name?: string,
    readonly blocks: ReadonlyArray<Block>
}

export function File(blocks: ArrayLike<Block>, name?: string): File {
    return { name, blocks: blocks as any };
}

export interface Frame {
    readonly header: string,
    readonly categories: Categories
}

export interface Block extends Frame {
    readonly saveFrames: Frame[]
}

export function Block(categories: Categories, header: string, saveFrames: Frame[] = []): Block {
    if (Object.keys(categories).some(k => k[0] !== '_')) {
        throw new Error(`Category names must start with '_'.`);
    }
    return { header, categories, saveFrames };
}

export function SafeFrame(categories: Categories, header: string): Frame {
    return { header, categories };
}

export type Categories = { readonly [name: string]: Category }

export interface Category {
    readonly rowCount: number,
    getField(name: string): Field | undefined
}

export function Category(rowCount: number, fields: { [name: string]: Field }): Category {
    return { rowCount, getField(name) { return fields[name]; } };
}

export namespace Category {
    export const Empty: Category = { rowCount: 0, getField(name: string) { return void 0; } };
}

/**
 * Implementation note:
 * Always implement this as a "plain" object so that the functions are "closures"
 * by default. This is to ensure that the schema access works without definiting
 * additional closures.
 */
export interface Field {
    readonly '@array': ArrayLike<any> | undefined
    readonly isDefined: boolean,
    readonly rowCount: number,

    str(row: number): string,
    int(row: number): number,
    float(row: number): number,

    valueKind(row: number): Column.ValueKind,

    areValuesEqual(rowA: number, rowB: number): boolean,
    stringEquals(row: number, value: string): boolean,

    toStringArray(params?: Column.ToArrayParams): ReadonlyArray<string>,
    toIntArray(params?: Column.ToArrayParams): ReadonlyArray<number>,
    toFloatArray(params?: Column.ToArrayParams): ReadonlyArray<number>
}

export function DefaultUndefinedField(rowCount: number): Field {
    return {
        '@array': void 0,
        isDefined: false,
        rowCount,
        str: row => '',
        int: row => 0,
        float: row => 0,

        valueKind: row => Column.ValueKind.NotPresent,
        areValuesEqual: (rowA, rowB) => true,
        stringEquals: (row, value) => value === null,

        toStringArray: (p) => createArray(rowCount, p).array,
        toIntArray: (p) => createArray(rowCount, p).array,
        toFloatArray: (p) => createArray(rowCount, p).array
    };
}

export function getMatrix(category: Category, field: string, rows: number, cols: number, row: number) {
    const ret: number[][] = [];
    for (let i = 0; i < rows; i++) {
        const r: number[] = [];
        for (let j = 0; j < cols; j++) {
            const f = category.getField(`${field}[${i + 1}][${j + 1}]`);
            r[j] = f ? f.float(row) : 0.0;
        }
        ret[i] = r;
    }
    return ret;
}

export function getVector(category: Category, field: string, rows: number, row: number) {
    const ret: number[] = [];
    for (let i = 0; i < rows; i++) {
        const f = category.getField(`${field}[${i + 1}]`);
        ret[i] = f ? f.float(row) : 0.0;
    }
    return ret;
}