/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface Column<T> {
    readonly '@type': Column.Type,
    readonly '@array': ArrayLike<any> | undefined,

    readonly isDefined: boolean,
    readonly rowCount: number,
    value(row: number): T,
    valueKind(row: number): Column.ValueKind,
    toArray(params?: Column.ToArrayParams): ReadonlyArray<T>,
    stringEquals(row: number, value: string): boolean,
    areValuesEqual(rowA: number, rowB: number): boolean
}

namespace Column {
    export type Type = typeof Type.str | typeof Type.int | typeof Type.float | typeof Type.vector | typeof Type.matrix

    export namespace Type {
        export const str = { '@type': '' as string, kind: 'str' as 'str' };
        export const int = { '@type': 0 as number, kind: 'int' as 'int' };
        export const float = { '@type': 0 as number, kind: 'float' as 'float' };
        export const vector = { '@type': [] as number[], kind: 'vector' as 'vector' };
        export const matrix = { '@type': [] as number[][], kind: 'matrix' as 'matrix' };
    }

    export interface ToArrayParams {
        array?: { new(size: number): ArrayLike<number> },
        start?: number,
        /** Last row (exclusive) */
        end?: number
    }

    export interface LambdaSpec<T extends Type> {
        value: (row: number) => T['@type'],
        rowCount: number,
        type: T,
        valueKind?: (row: number) => ValueKind,
    }

    export interface ArraySpec<T extends Type> {
        array: ArrayLike<T['@type']>,
        type: T,
        valueKind?: (row: number) => ValueKind
    }

    export const enum ValueKind { Present = 0, NotPresent = 1, Unknown = 2 }

    export function Undefined<T extends Type>(rowCount: number, type: T): Column<T['@type']> {
        return constColumn(type['@type'], rowCount, type, ValueKind.NotPresent);
    }

    export function ofConst<T extends Type>(v: T['@type'], rowCount: number, type: T): Column<T['@type']> {
        return constColumn(v, rowCount, type, ValueKind.Present);
    }

    export function ofLambda<T extends Type>(spec: LambdaSpec<T>): Column<T['@type']> {
        return lambdaColumn(spec);
    }

    export function ofArray<T extends Column.Type>(spec: Column.ArraySpec<T>): Column<T['@type']> {
        return arrayColumn(spec);
    }

    export function window<T>(column: Column<T>, start: number, end: number) {
        return windowColumn(column, start, end);
    }

    /** Makes the column backned by an array. Useful for columns that accessed often. */
    export function asArrayColumn<T>(c: Column<T>, array?: ToArrayParams['array']): Column<T> {
        if (c['@array']) return c;
        if (!c.isDefined) return Undefined(c.rowCount, c['@type']) as any as Column<T>;
        return arrayColumn({ array: c.toArray({ array }), type: c['@type'] as any, valueKind: c.valueKind });
    }
}

export default Column;

function constColumn<T extends Column.Type>(v: T['@type'], rowCount: number, type: T, valueKind: Column.ValueKind): Column<T['@type']> {
    const value: Column<T['@type']>['value'] = row => v;
    return {
        '@type': type,
        '@array': void 0,
        isDefined: valueKind === Column.ValueKind.Present,
        rowCount,
        value,
        valueKind: row => valueKind,
        toArray: params => {
            const { array } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v;
            return array;
        },
        stringEquals: type.kind === 'str'
            ? (row, value) => value === v
            : type.kind === 'float' || type.kind === 'int'
                ? (row, value) => +value === v
                : (row, value) => false,
        areValuesEqual: (rowA, rowB) => true
    }
}

function lambdaColumn<T extends Column.Type>({ value, valueKind, rowCount, type }: Column.LambdaSpec<T>): Column<T['@type']> {
    return {
        '@type': type,
        '@array': void 0,
        isDefined: true,
        rowCount,
        value,
        valueKind: valueKind ? valueKind : row => Column.ValueKind.Present,
        toArray: params => {
            const { array, start } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = value(i + start);
            return array;
        },
        stringEquals: type.kind === 'str'
            ? (row, v) => value(row) === v
            : type.kind === 'float' || type.kind === 'int'
            ? (row, v) => value(row) === +v
            : (row, value) => false,
        areValuesEqual: (rowA, rowB) => value(rowA) === value(rowB)
    }
}

function arrayColumn<T extends Column.Type>({ array, type, valueKind }: Column.ArraySpec<T>): Column<T['@type']> {
    const rowCount = array.length;
    const value: Column<T['@type']>['value'] = type.kind === 'str'
        ? row => { const v = array[row]; return typeof v === 'string' ? v : '' + v; }
        : row => array[row];

    const isTyped = isTypedArray(array);
    return {
        '@type': type,
        '@array': array,
        isDefined: true,
        rowCount,
        value,
        valueKind: valueKind ? valueKind : row => Column.ValueKind.Present,
        toArray: type.kind === 'str'
            ? params => {
                const { start, end } = getArrayBounds(rowCount, params);
                const ret = new (params && typeof params.array !== 'undefined' ? params.array : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) {
                    const v = array[start + i];
                    ret[i] = typeof v === 'string' ? v : '' + v;
                }
                return ret;
            }
            : isTyped
            ? params => typedArrayWindow(array, params) as any as ReadonlyArray<T>
            : params => {
                const { start, end } = getArrayBounds(rowCount, params);
                if (start === 0 && end === array.length) return array as ReadonlyArray<T['@type']>;
                const ret = new (params && typeof params.array !== 'undefined' ? params.array : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) ret[i] = array[start + i];
                return ret;
            },
        stringEquals: type.kind === 'int' || type.kind === 'float'
            ? (row, value) => (array as any)[row] === +value
            : type.kind === 'str'
            ? (row, value) => { const v = array[row]; return typeof v === 'string' ? v === value : +v === +value; }
            : (row, value) => false,
        areValuesEqual: (rowA, rowB) => array[rowA] === array[rowB]
    }
}

function windowColumn<T>(column: Column<T>, start: number, end: number) {
    if (!column.isDefined) return Column.Undefined(end - start, column['@type']);
    if (column['@array'] && isTypedArray(column['@array'])) return windowTyped(column, start, end);
    return windowFull(column, start, end);
}

function windowTyped<T>(c: Column<T>, start: number, end: number): Column<T> {
    const array = typedArrayWindow(c['@array'], { start, end });
    return arrayColumn({ array, type: c['@type'], valueKind: c.valueKind }) as any;
}

function windowFull<T>(c: Column<T>, start: number, end: number): Column<T> {
    const v = c.value, vk = c.valueKind, se = c.stringEquals, ave = c.areValuesEqual;
    const value: Column<T>['value'] = start === 0 ? v : row => v(row + start);
    const rowCount = end - start;
    return {
        '@type': c['@type'],
        '@array': void 0,
        isDefined: c.isDefined,
        rowCount,
        value,
        valueKind: start === 0 ? vk : row => vk(row + start),
        toArray: params => {
            const { array } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v(i + start);
            return array;
        },
        stringEquals: start === 0 ? se : (row, value) => se(row + start, value),
        areValuesEqual: start === 0 ? ave : (rowA, rowB) => ave(rowA + start, rowB + start)
    };
}

/** A helped function for Column.toArray */
export function getArrayBounds(rowCount: number, params?: Column.ToArrayParams) {
    const start = params && typeof params.start !== 'undefined' ? Math.max(Math.min(params.start, rowCount - 1), 0) : 0;
    const end = params && typeof params.end !== 'undefined' ? Math.min(params.end, rowCount) : rowCount;
    return { start, end };
}

/** A helped function for Column.toArray */
export function createArray(rowCount: number, params?: Column.ToArrayParams) {
    const c = params && typeof params.array !== 'undefined' ? params.array : Array;
    const { start, end } = getArrayBounds(rowCount, params);
    return { array: new c(end - start) as any[], start, end };
}

/** A helped function for Column.toArray */
export function fillArrayValues(value: (row: number) => any, target: any[], start: number) {
    for (let i = 0, _e = target.length; i < _e; i++) target[i] = value(start + i);
    return target;
}

/** A helped function for Column.toArray */
export function createAndFillArray(rowCount: number, value: (row: number) => any, params?: Column.ToArrayParams) {
    const { array, start } = createArray(rowCount, params);
    return fillArrayValues(value, array, start);
}

export function isTypedArray(data: any): boolean {
    return !!data.buffer && typeof data.byteLength === 'number' && typeof data.BYTES_PER_ELEMENT === 'number';
}

export function typedArrayWindow(data: any, params?: Column.ToArrayParams): ReadonlyArray<number> {
    const { constructor, buffer, length, byteOffset, BYTES_PER_ELEMENT } = data;
    const { start, end } = getArrayBounds(length, params);
    if (start === 0 && end === length) return data;
    return new constructor(buffer, byteOffset + BYTES_PER_ELEMENT * start, Math.min(length, end - start));
}