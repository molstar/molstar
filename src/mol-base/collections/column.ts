/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ColumnType = typeof ColumnType.str | typeof ColumnType.pooledStr | typeof ColumnType.int | typeof ColumnType.float | typeof ColumnType.vector | typeof ColumnType.matrix

export namespace ColumnType {
    export const str = { '@type': '' as string, kind: 'str' as 'str', isScalar: false, isString: true };
    export const pooledStr = { '@type': '' as string, kind: 'pooled-str' as 'pooled-str', isScalar: false, isString: true };
    export const int = { '@type': 0 as number, kind: 'int' as 'int', isScalar: true, isString: false };
    export const float = { '@type': 0 as number, kind: 'float' as 'float', isScalar: true, isString: false };
    export const vector = { '@type': [] as number[], kind: 'vector' as 'vector', isScalar: false, isString: false };
    export const matrix = { '@type': [] as number[][], kind: 'matrix' as 'matrix', isScalar: false, isString: false };
}

export interface ToArrayParams {
    array?: { new(size: number): ArrayLike<number> },
    /** First row */
    start?: number,
    /** Last row (exclusive) */
    end?: number
}

export interface Column<T> {
    readonly '@type': ColumnType,
    readonly isDefined: boolean,
    readonly rowCount: number,
    value(row: number): T,
    isValueDefined(row: number): boolean,
    toArray(params?: ToArrayParams): ReadonlyArray<T>,
    stringEquals(row: number, value: string): boolean,
    areValuesEqual(rowA: number, rowB: number): boolean
}

export function UndefinedColumn<T extends ColumnType>(rowCount: number, type: T): Column<T['@type']> {
    const v = type.isString ? '' : 0;
    const value: Column<T['@type']>['value'] = type.isString ? row => '' : row => 0;
    return {
        '@type': type,
        isDefined: false,
        rowCount,
        value,
        isValueDefined: row => false,
        toArray: params => {
            const { array } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v;
            return array;
        },
        stringEquals: (row, value) => !value,
        areValuesEqual: (rowA, rowB) => true
    }
}

export function ConstColumn<T extends ColumnType>(v: T['@type'], rowCount: number, type: T): Column<T['@type']> {
    const value: Column<T['@type']>['value'] = row => v;
    return {
        '@type': type,
        isDefined: true,
        rowCount,
        value,
        isValueDefined: row => false,
        toArray: params => {
            const { array } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v;
            return array;
        },
        stringEquals: type.isString
            ? (row, value) => value === v
            : type.isScalar
            ? (row, value) => +value === v
            : (row, value) => false,
        areValuesEqual: (rowA, rowB) => true
    }
}

export interface ArrayColumnSpec<T extends ColumnType> {
    array: ArrayLike<T['@type']>,
    type: T,
    isValueDefined?: (row: number) => boolean
}

export function ArrayColumn<T extends ColumnType>({ array, type, isValueDefined }: ArrayColumnSpec<T>): Column<T['@type']> {
    const rowCount = array.length;
    const value: Column<T['@type']>['value'] = row => array[row];
    const isTyped = isTypedArray(array);
    return {
        '@type': type,
        isDefined: true,
        rowCount,
        value,
        isValueDefined: isValueDefined ? isValueDefined : row => true,
        toArray: isTyped
            ? params => typedArrayWindow(array, params) as any as ReadonlyArray<T>
            : params => {
                const { start, end } = getArrayBounds(rowCount, params);
                if (start === 0 && end === array.length) return array as ReadonlyArray<T['@type']>;
                const ret = new (params && typeof params.array !== 'undefined' ? params.array : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) ret[i] = array[start + i];
                return ret;
            },
        stringEquals: type.isScalar
            ? (row, value) => (array as any)[row] === +value
            : type.isString
            ? (row, value) => array[row] === value
            : (row, value) => false,
        areValuesEqual: (rowA, rowB) => array[rowA] === array[rowB]
    }
}

/** Makes the column backned by an array. Useful for columns that accessed often. */
export function toArrayColumn<T>(c: Column<T>): Column<T> {
    if (!c.isDefined) return UndefinedColumn(c.rowCount, c['@type']) as any as Column<T>;
    return ArrayColumn({ array: c.toArray(), type: c['@type'] as any, isValueDefined: c.isValueDefined });
}

/** A helped function for Column.toArray */
export function getArrayBounds(rowCount: number, params?: ToArrayParams) {
    const start = params && typeof params.start !== 'undefined' ? Math.max(Math.min(params.start, rowCount - 1), 0) : 0;
    const end = params && typeof params.end !== 'undefined' ? Math.min(params.end, rowCount) : rowCount;
    return { start, end };
}

/** A helped function for Column.toArray */
export function createArray(rowCount: number, params?: ToArrayParams) {
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
export function createAndFillArray(rowCount: number, value: (row: number) => any, params?: ToArrayParams) {
    const { array, start } = createArray(rowCount, params);
    return fillArrayValues(value, array, start);
}

export function isTypedArray(data: any) {
    return data.buffer && typeof data.byteLength === 'number' && data.BYTES_PER_ELEMENT;
}

export function typedArrayWindow(data: any, params?: ToArrayParams): ReadonlyArray<number> {
    const { constructor, buffer, length, byteOffset, BYTES_PER_ELEMENT } = data;
    const { start, end } = getArrayBounds(length, params);
    if (start === 0 && end === length) return data;
    return new constructor(buffer, byteOffset + BYTES_PER_ELEMENT * start, Math.min(length, end - start));
}