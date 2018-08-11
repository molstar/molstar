/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from 'mol-data/db'
import { Tensor } from 'mol-math/linear-algebra'
import { getNumberType, NumberType } from '../common/text/number-parser';

export interface CifFile {
    readonly name?: string,
    readonly blocks: ReadonlyArray<CifBlock>
}

export function CifFile(blocks: ArrayLike<CifBlock>, name?: string): CifFile {
    return { name, blocks: blocks as any };
}

export interface CifFrame {
    readonly header: string,
    // Category names stored separately so that the ordering can be preserved.
    readonly categoryNames: ReadonlyArray<string>,
    readonly categories: CifCategories
}

export interface CifBlock extends CifFrame {
    readonly saveFrames: CifFrame[]
}

export function CifBlock(categoryNames: string[], categories: CifCategories, header: string, saveFrames: CifFrame[] = []): CifBlock {
    return { categoryNames, header, categories, saveFrames };
}

export function CifSafeFrame(categoryNames: string[], categories: CifCategories, header: string): CifFrame {
    return { categoryNames, header, categories };
}

export type CifCategories = { readonly [name: string]: CifCategory }

export interface CifCategory {
    readonly rowCount: number,
    readonly name: string,
    readonly fieldNames: ReadonlyArray<string>,
    getField(name: string): CifField | undefined
}

export function CifCategory(name: string, rowCount: number, fieldNames: string[], fields: { [name: string]: CifField }): CifCategory {
    return { rowCount, name, fieldNames: [...fieldNames], getField(name) { return fields[name]; } };
}

export namespace CifCategory {
    export function empty(name: string): CifCategory {
        return { rowCount: 0, name, fieldNames: [], getField(name: string) { return void 0; } };
    };
}

/**
 * Implementation note:
 * Always implement without using "this." in any of the interface functions.
 * This is to ensure that the functions can invoked without having to "bind" them.
 */
export interface CifField {
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

export function getTensor(category: CifCategory, field: string, space: Tensor.Space, row: number, zeroIndexed: boolean): Tensor.Data {
    const ret = space.create();
    const offset = zeroIndexed ? 0 : 1;

    if (space.rank === 1) {
        const rows = space.dimensions[0];
        for (let i = 0; i < rows; i++) {
            const f = category.getField(`${field}[${i + offset}]`);
            space.set(ret, i, !!f ? f.float(row) : 0.0);
        }
    } else if (space.rank === 2) {
        const rows = space.dimensions[0], cols = space.dimensions[1];
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const f = category.getField(`${field}[${i + offset}][${j + offset}]`);
                space.set(ret, i, j, !!f ? f.float(row) : 0.0);
            }
        }
    } else if (space.rank === 3) {
        const d0 = space.dimensions[0], d1 = space.dimensions[1], d2 = space.dimensions[2];
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) {
                for (let k = 0; k < d2; k++) {
                    const f = category.getField(`${field}[${i + offset}][${j + offset}][${k + offset}]`);
                    space.set(ret, i, j, k, !!f ? f.float(row) : 0.0);
                }
            }
        }
    } else throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.');
    return ret;
}

export function getCifFieldType(field: CifField): Column.Schema.Int | Column.Schema.Float | Column.Schema.Str {
    let floatCount = 0, hasString = false;
    for (let i = 0, _i = field.rowCount; i < _i; i++) {
        const k = field.valueKind(i);
        if (k !== Column.ValueKind.Present) continue
        const type = getNumberType(field.str(i));
        if (type === NumberType.Int) continue;
        else if (type === NumberType.Float) floatCount++;
        else { hasString = true; break; }
    }

    if (hasString) return Column.Schema.str;
    if (floatCount > 0) return Column.Schema.float;
    return Column.Schema.int;
}