/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from 'mol-data/db'
import { Tensor } from 'mol-math/linear-algebra'

export interface File {
    readonly name?: string,
    readonly blocks: ReadonlyArray<Block>
}

export function File(blocks: ArrayLike<Block>, name?: string): File {
    return { name, blocks: blocks as any };
}

export interface Frame {
    readonly header: string,
    // Category names stored separately so that the ordering can be preserved.
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
 * Always implement without using "this." in any of the interface functions.
 * This is to ensure that the functions can invoked without having to "bind" them.
 */
export interface Field {
    readonly '@array': ArrayLike<any> | undefined
    readonly isDefined: boolean,
    readonly rowCount: number,

    str(row: number): string,
    int(row: number): number,
    float(row: number): number,
    list<T extends number|string>(row: number): T[],

    valueKind(row: number): Column.ValueKind,

    areValuesEqual(rowA: number, rowB: number): boolean,

    toStringArray(params?: Column.ToArrayParams<string>): ReadonlyArray<string>,
    toIntArray(params?: Column.ToArrayParams<number>): ReadonlyArray<number>,
    toFloatArray(params?: Column.ToArrayParams<number>): ReadonlyArray<number>
    toListArray<T extends number|string>(params?: Column.ToArrayParams<T[]>): ReadonlyArray<T[]>
}

export function getTensor(category: Category, field: string, space: Tensor.Space, row: number): Tensor {
    const ret = space.create();
    if (space.rank === 1) {
        const rows = space.dimensions[0];
        for (let i = 0; i < rows; i++) {
            const f = category.getField(`${field}[${i + 1}]`);
            space.set(ret, i, !!f ? f.float(row) : 0.0);
        }
    } else if (space.rank === 2) {
        const rows = space.dimensions[0], cols = space.dimensions[1];
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const f = category.getField(`${field}[${i + 1}][${j + 1}]`);
                space.set(ret, i, j, !!f ? f.float(row) : 0.0);
            }
        }
    } else if (space.rank === 3) {
        const d0 = space.dimensions[0], d1 = space.dimensions[1], d2 = space.dimensions[2];
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) {
                for (let k = 0; k < d2; k++) {
                    const f = category.getField(`${field}[${i + 1}][${j + 1}][${k + 1}]`);
                    space.set(ret, i, j, k, !!f ? f.float(row) : 0.0);
                }
            }
        }
    } else throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.');
    return ret;
}