/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ColumnType = typeof ColumnType.str | typeof ColumnType.pooledStr | typeof ColumnType.int | typeof ColumnType.float

export namespace ColumnType {
    export const str = { '@type': '' as string, kind: 'str' as 'str' };
    export const pooledStr = { '@type': '' as string, kind: 'pooled-str' as 'pooled-str' };
    export const int = { '@type': 0 as number, kind: 'int' as 'int' };
    export const float = { '@type': 0 as number, kind: 'float' as 'float' };
}

export interface ToArrayParams {
    array?: { new(size: number): ArrayLike<number> },
    /** First row */
    start?: number,
    /** Last row (exclusive) */
    end?: number
}

export interface Column<T> {
    readonly isDefined: boolean,
    readonly rowCount: number,
    value(row: number): T,
    isValueDefined(row: number): boolean,
    toArray(params?: ToArrayParams): ReadonlyArray<T>,
    stringEquals(row: number, value: string): boolean,
    areValuesEqual(rowA: number, rowB: number): boolean
}

export function UndefinedColumn<T extends ColumnType>(rowCount: number, type: T): Column<T['@type']> {
    const value: Column<T['@type']>['value'] = type.kind === 'str' ? row => '' : row => 0;
    return {
        isDefined: false,
        rowCount,
        value,
        isValueDefined(row) { return false; },
        toArray(params) {
            const { array } = createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = value(0)
            return array;
        },
        stringEquals(row, value) { return !value; },
        areValuesEqual(rowA, rowB) { return true; }
    }
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


