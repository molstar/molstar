/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, ColumnHelpers, Table } from '../../../mol-data/db';
import { Tensor } from '../../../mol-math/linear-algebra';
import { getNumberType, NumberType, parseInt as fastParseInt, parseFloat as fastParseFloat } from '../common/text/number-parser';
import { Encoding } from '../../common/binary-cif';
import { Tokens } from '../common/text/tokenizer';
import { areValuesEqualProvider } from '../common/text/column/token';

export interface CifFile {
    readonly name?: string,
    readonly blocks: ReadonlyArray<CifBlock>
}

export function CifFile(blocks: ArrayLike<CifBlock>, name?: string): CifFile {
    return { name, blocks: blocks as any };
}

export interface CifFrame {
    readonly header: string,
    /** Category names, stored separately so that the ordering can be preserved. */
    readonly categoryNames: ReadonlyArray<string>,
    readonly categories: CifCategories
}

export interface CifBlock extends CifFrame {
    readonly saveFrames: CifFrame[]
    getField(name: string): CifField | undefined
}

export function CifBlock(categoryNames: string[], categories: CifCategories, header: string, saveFrames: CifFrame[] = []): CifBlock {
    return {
        categoryNames, header, categories, saveFrames,
        getField(name: string) {
            const [ category, field ] = name.split('.');
            return categories[category].getField(field || '');
        }
    };
}

export function CifSaveFrame(categoryNames: string[], categories: CifCategories, header: string): CifFrame {
    return { categoryNames, header, categories };
}

export type CifAliases = { readonly [name: string]: string[] }
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

    export type SomeFields<S> = { [P in keyof S]?: CifField }
    export type Fields<S> = { [P in keyof S]: CifField }

    export function ofFields(name: string, fields: { [name: string]: CifField | undefined }): CifCategory {
        const fieldNames = Object.keys(fields);
        return {
            rowCount: fieldNames.length > 0 ? fields[fieldNames[0]]!.rowCount : 0,
            name,
            fieldNames,
            getField(name) { return fields[name]; }
        };
    }

    export function ofTable(name: string, table: Table<any>) {
        const fields: { [name: string]: CifField | undefined } = {};
        for (const name of table._columns) {
            fields[name] = CifField.ofColumn(table[name]);
        }
        return ofFields(name, fields);
    }
}

/**
 * Implementation note:
 * Always implement without using "this." in any of the interface functions.
 * This is to ensure that the functions can invoked without having to "bind" them.
 */
export interface CifField {
    readonly __array: ArrayLike<any> | undefined,
    readonly binaryEncoding: Encoding[] | undefined,
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

export namespace CifField {
    export function ofString(value: string) {
        return ofStrings([value]);
    }

    export function ofStrings(values: ArrayLike<string>): CifField {
        const rowCount = values.length;
        const str: CifField['str'] = row => { const ret = values[row]; if (!ret || ret === '.' || ret === '?') return ''; return ret; };
        const int: CifField['int'] = row => { const v = values[row]; return fastParseInt(v, 0, v.length) || 0; };
        const float: CifField['float'] = row => { const v = values[row]; return fastParseFloat(v, 0, v.length) || 0; };
        const valueKind: CifField['valueKind'] = row => {
            const v = values[row], l = v.length;
            if (l > 1) return Column.ValueKind.Present;
            if (l === 0) return Column.ValueKind.NotPresent;
            const c = v.charCodeAt(0);
            if (c === 46 /* . */) return Column.ValueKind.NotPresent;
            if (c === 63 /* ? */) return Column.ValueKind.Unknown;
            return Column.ValueKind.Present;
        };

        return {
            __array: void 0,
            binaryEncoding: void 0,
            isDefined: true,
            rowCount,
            str,
            int,
            float,
            valueKind,
            areValuesEqual: (rowA, rowB) => values[rowA] === values[rowB],
            toStringArray: params => params ? ColumnHelpers.createAndFillArray(rowCount, str, params) : values as string[],
            toIntArray: params => ColumnHelpers.createAndFillArray(rowCount, int, params),
            toFloatArray: params => ColumnHelpers.createAndFillArray(rowCount, float, params)
        };
    }

    export function ofNumbers(values: ArrayLike<number>): CifField {
        const rowCount = values.length;
        const str: CifField['str'] = row => { return '' + values[row]; };
        const float: CifField['float'] = row => values[row];
        const valueKind: CifField['valueKind'] = row => Column.ValueKind.Present;

        const toFloatArray = (params: Column.ToArrayParams<number>) => {
            if (!params || params.array && values instanceof params.array) {
                return values as number[];
            } else {
                return ColumnHelpers.createAndFillArray(rowCount, float, params);
            }
        };

        return {
            __array: void 0,
            binaryEncoding: void 0,
            isDefined: true,
            rowCount,
            str,
            int: float,
            float,
            valueKind,
            areValuesEqual: (rowA, rowB) => values[rowA] === values[rowB],
            toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
            toIntArray: toFloatArray,
            toFloatArray
        };
    }

    export function ofTokens(tokens: Tokens): CifField {
        const { data, indices, count: rowCount } = tokens;

        const str: CifField['str'] = row => {
            const ret = data.substring(indices[2 * row], indices[2 * row + 1]);
            if (ret === '.' || ret === '?') return '';
            return ret;
        };

        const int: CifField['int'] = row => {
            return fastParseInt(data, indices[2 * row], indices[2 * row + 1]) || 0;
        };

        const float: CifField['float'] = row => {
            return fastParseFloat(data, indices[2 * row], indices[2 * row + 1]) || 0;
        };

        const valueKind: CifField['valueKind'] = row => {
            const s = indices[2 * row], l = indices[2 * row + 1] - s;
            if (l > 1) return Column.ValueKind.Present;
            if (l === 0) return Column.ValueKind.NotPresent;
            const v = data.charCodeAt(s);
            if (v === 46 /* . */) return Column.ValueKind.NotPresent;
            if (v === 63 /* ? */) return Column.ValueKind.Unknown;
            return Column.ValueKind.Present;
        };

        return {
            __array: void 0,
            binaryEncoding: void 0,
            isDefined: true,
            rowCount,
            str,
            int,
            float,
            valueKind,
            areValuesEqual: areValuesEqualProvider(tokens),
            toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
            toIntArray: params => ColumnHelpers.createAndFillArray(rowCount, int, params),
            toFloatArray: params => ColumnHelpers.createAndFillArray(rowCount, float, params)
        };
    }

    export function ofColumn(column: Column<any>): CifField {
        const { rowCount, valueKind, areValuesEqual, isDefined } = column;

        let str: CifField['str'];
        let int: CifField['int'];
        let float: CifField['float'];

        switch (column.schema.valueType) {
            case 'float':
            case 'int':
                str = row => { return '' + column.value(row); };
                int = column.value;
                float = column.value;
                break;
            case 'str':
                str = column.value;
                int = row => { const v = column.value(row); return fastParseInt(v, 0, v.length) || 0; };
                float = row => { const v = column.value(row); return fastParseFloat(v, 0, v.length) || 0; };
                break;
            case 'list':
                const { separator } = column.schema;
                str = row => column.value(row).join(separator);
                int = row => NaN;
                float = row => NaN;
                break;
            default:
                throw new Error(`unsupported valueType '${column.schema.valueType}'`);
        }

        return {
            __array: void 0,
            binaryEncoding: void 0,
            isDefined,
            rowCount,
            str,
            int,
            float,
            valueKind,
            areValuesEqual,
            toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
            toIntArray: params => ColumnHelpers.createAndFillArray(rowCount, int, params),
            toFloatArray: params => ColumnHelpers.createAndFillArray(rowCount, float, params)
        };
    }
}

export function tensorFieldNameGetter(field: string, rank: number, zeroIndexed: boolean, namingVariant: 'brackets' | 'underscore') {
    const offset = zeroIndexed ? 0 : 1;
    switch (rank) {
        case 1:
            return namingVariant === 'brackets'
                ? (i: number) => `${field}[${i + offset}]`
                : (i: number) => `${field}_${i + offset}`;
        case 2:
            return namingVariant === 'brackets'
                ? (i: number, j: number) => `${field}[${i + offset}][${j + offset}]`
                : (i: number, j: number) => `${field}_${i + offset}${j + offset}`;
        case 3:
            return namingVariant === 'brackets'
                ? (i: number, j: number, k: number) => `${field}[${i + offset}][${j + offset}][${k + offset}]`
                : (i: number, j: number, k: number) => `${field}_${i + offset}${j + offset}${k + offset}`;
        default:
            throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.');
    }
}

export function getTensor(category: CifCategory, space: Tensor.Space, row: number, getName: (...args: number[]) => string): Tensor.Data {
    const ret = space.create();

    if (space.rank === 1) {
        const rows = space.dimensions[0];
        for (let i = 0; i < rows; i++) {
            const f = category.getField(getName(i));
            space.set(ret, i, !!f ? f.float(row) : 0.0);
        }
    } else if (space.rank === 2) {
        const rows = space.dimensions[0], cols = space.dimensions[1];
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const f = category.getField(getName(i, j));
                space.set(ret, i, j, !!f ? f.float(row) : 0.0);
            }
        }
    } else if (space.rank === 3) {
        const d0 = space.dimensions[0], d1 = space.dimensions[1], d2 = space.dimensions[2];
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) {
                for (let k = 0; k < d2; k++) {
                    const f = category.getField(getName(i, j, k));
                    space.set(ret, i, j, k, !!f ? f.float(row) : 0.0);
                }
            }
        }
    } else {
        throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.');
    }
    return ret;
}

export function getCifFieldType(field: CifField): Column.Schema.Int | Column.Schema.Float | Column.Schema.Str {
    let floatCount = 0, hasStringOrScientific = false, undefinedCount = 0;
    for (let i = 0, _i = field.rowCount; i < _i; i++) {
        const k = field.valueKind(i);
        if (k !== Column.ValueKind.Present) {
            undefinedCount++;
            continue;
        }
        const type = getNumberType(field.str(i));
        if (type === NumberType.Int) continue;
        else if (type === NumberType.Float) floatCount++;
        else { hasStringOrScientific = true; break; }
    }

    // numbers in scientific notation and plain text are not distinguishable
    if (hasStringOrScientific || undefinedCount === field.rowCount) return Column.Schema.str;
    if (floatCount > 0) return Column.Schema.float;
    return Column.Schema.int;
}