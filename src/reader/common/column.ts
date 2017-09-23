/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ArrayType = string[] | number[] | Float32Array | Float64Array | Int8Array | Int16Array | Int32Array | Uint8Array | Uint16Array | Uint32Array
export type ColumnType = typeof ColumnType.str | typeof ColumnType.int | typeof ColumnType.float

export namespace ColumnType {
    export const str = { '@type': '' as string, kind: 'str' as 'str' };
    export const int = { '@type': 0 as number, kind: 'int' as 'int' };
    export const float = { '@type': 0 as number, kind: 'float' as 'float' };
}

export interface Column<T> {
    readonly isColumnDefined: boolean,
    readonly rowCount: number,
    value(row: number): T,
    toArray(ctor?: (size: number) => ArrayType, startRow?: number, endRowExclusive?: number): ReadonlyArray<T>
}

export function UndefinedColumn<T extends ColumnType>(rowCount: number, type: T): Column<T['@type']> {
    const value: Column<T['@type']>['value'] = type.kind === 'str' ? row => '' : row => 0;
    return {
        isColumnDefined: false,
        rowCount,
        value,
        toArray(ctor, s, e) {
            const { array } = createArray(rowCount, ctor, s, e);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = value(0)
            return array;
        }
    }
}

/** A helped function for Column.toArray */
export function createArray(rowCount: number, ctor?: (size: number) => ArrayType, start?: number, end?: number) {
    const c = typeof ctor !== 'undefined' ? ctor : (s: number) => new Array(s);
    const s = typeof start !== 'undefined' ? Math.max(Math.min(start, rowCount - 1), 0) : 0;
    const e = typeof end !== 'undefined' ? Math.min(end, rowCount) : rowCount;
    return { array: c(e - s) as any[], start: s, end: e };
}