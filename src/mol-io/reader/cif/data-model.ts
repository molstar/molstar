/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-base/collections/database'

export interface File {
    readonly name?: string,
    readonly blocks: ReadonlyArray<Block>
}

export function File(blocks: ArrayLike<Block>, name?: string): File {
    return { name, blocks: blocks as any };
}

export interface Frame {
    readonly header: string,
    readonly categoryNames: ReadonlyArray<string>,
    readonly categories: Categories
}

export interface Block extends Frame {
    readonly saveFrames: Frame[]
}

export function Block(categoryNames: string[], categories: Categories, header: string, saveFrames: Frame[] = []): Block {
    return { categoryNames, header, categories, saveFrames };
}

export function SafeFrame(categoryNames: string[], categories: Categories, header: string): Frame {
    return { categoryNames, header, categories };
}

export type Categories = { readonly [name: string]: Category }

export interface Category {
    readonly rowCount: number,
    readonly name: string,
    readonly fieldNames: ReadonlyArray<string>,
    getField(name: string): Field | undefined
}

export function Category(name: string, rowCount: number, fieldNames: string[], fields: { [name: string]: Field }): Category {
    return { rowCount, name, fieldNames: [...fieldNames], getField(name) { return fields[name]; } };
}

export namespace Category {
    export function empty(name: string): Category {
        return { rowCount: 0, name, fieldNames: [], getField(name: string) { return void 0; } };
    };
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

    toStringArray(params?: Column.ToArrayParams<string>): ReadonlyArray<string>,
    toIntArray(params?: Column.ToArrayParams<number>): ReadonlyArray<number>,
    toFloatArray(params?: Column.ToArrayParams<number>): ReadonlyArray<number>
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