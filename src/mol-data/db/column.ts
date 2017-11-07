/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as ColumnHelpers from './column-helpers'
import { Tensor as Tensors } from 'mol-math/linear-algebra'

interface Column<T> {
    readonly schema: Column.Schema,
    readonly '@array': ArrayLike<any> | undefined,

    readonly isDefined: boolean,
    readonly rowCount: number,
    value(row: number): T,
    valueKind(row: number): Column.ValueKind,
    toArray(params?: Column.ToArrayParams<T>): ArrayLike<T>,
    areValuesEqual(rowA: number, rowB: number): boolean
}

namespace Column {
    export type ArrayCtor<T> = { new(size: number): ArrayLike<T> }

    export type Schema<T = any> = Schema.Str | Schema.Int | Schema.Float | Schema.Coordinate | Schema.Aliased<T> | Schema.Tensor

    export namespace Schema {
        type Base<T extends string> = { valueType: T }
        export type Str = { '@type': 'str', T: string } & Base<'str'>
        export type Int = { '@type': 'int', T: number } & Base<'int'>
        export type Float = { '@type': 'float', T: number } & Base<'float'>
        export type Coordinate = { '@type': 'coord', T: number } & Base<'float'>

        export type Tensor = { '@type': 'tensor', T: Tensors, space: Tensors.Space } & Base<'tensor'>
        export type Aliased<T> = { '@type': 'aliased', T: T } & Base<'str' | 'int'>

        export const str: Str = { '@type': 'str', T: '', valueType: 'str' };
        export const int: Int = { '@type': 'int', T: 0, valueType: 'int' };
        export const coord: Coordinate = { '@type': 'coord', T: 0, valueType: 'float' };
        export const float: Float = { '@type': 'float', T: 0, valueType: 'float' };

        export function tensor(space: Tensors.Space): Tensor { return { '@type': 'tensor', T: space.create(), space, valueType: 'tensor' }; }
        export function vector(dim: number): Tensor { return tensor(Tensors.Vector(dim)); }
        export function matrix(rows: number, cols: number): Tensor { return tensor(Tensors.ColumnMajorMatrix(rows, cols)); }

        export function aliased<T>(t: Str | Int): Aliased<T> { return t as any as Aliased<T>; }
    }

    export interface ToArrayParams<T> {
        array?: ArrayCtor<T>,
        start?: number,
        /** Last row (exclusive) */
        end?: number
    }

    export interface LambdaSpec<T extends Schema> {
        value: (row: number) => T['T'],
        rowCount: number,
        schema: T,
        valueKind?: (row: number) => ValueKind,
    }

    export interface ArraySpec<T extends Schema> {
        array: ArrayLike<T['T']>,
        schema: T,
        valueKind?: (row: number) => ValueKind
    }

    export interface MapSpec<S extends Schema, T extends Schema> {
        f: (v: S['T']) => T['T'],
        schema: T,
        valueKind?: (row: number) => ValueKind,
    }

    export const enum ValueKind { Present = 0, NotPresent = 1, Unknown = 2 }

    export function Undefined<T extends Schema>(rowCount: number, schema: T): Column<T['T']> {
        return constColumn(schema['T'], rowCount, schema, ValueKind.NotPresent);
    }

    export function ofConst<T extends Schema>(v: T['T'], rowCount: number, type: T): Column<T['T']> {
        return constColumn(v, rowCount, type, ValueKind.Present);
    }

    export function ofLambda<T extends Schema>(spec: LambdaSpec<T>): Column<T['T']> {
        return lambdaColumn(spec);
    }

    export function ofArray<T extends Column.Schema>(spec: Column.ArraySpec<T>): Column<T['T']> {
        return arrayColumn(spec);
    }

    export function ofIntArray(array: ArrayLike<number>) {
        return arrayColumn({ array, schema: Schema.int });
    }

    export function ofFloatArray(array: ArrayLike<number>) {
        return arrayColumn({ array, schema: Schema.float });
    }

    export function ofStringArray(array: ArrayLike<string>) {
        return arrayColumn({ array, schema: Schema.str });
    }

    export function window<T>(column: Column<T>, start: number, end: number) {
        return windowColumn(column, start, end);
    }

    export function view<T>(column: Column<T>, indices: ArrayLike<number>, checkIndentity = true) {
        return columnView(column, indices, checkIndentity);
    }

    /** A map of the 1st occurence of each value. */
    export function createFirstIndexMap<T>(column: Column<T>) {
        return createFirstIndexMapOfColumn(column);
    }

    export function mapToArray<T, S>(column: Column<T>, f: (v: T) => S, ctor?: ArrayCtor<S>): ArrayLike<S> {
        return mapToArrayImpl(column, f, ctor || Array);
    }

    export function areEqual<T>(a: Column<T>, b: Column<T>) {
        return areColumnsEqual(a, b);
    }

    export function indicesOf<T>(c: Column<T>, test: (e: T) => boolean) {
        return columnIndicesOf(c, test);
    }

    /** Makes the column backned by an array. Useful for columns that accessed often. */
    export function asArrayColumn<T>(c: Column<T>, array?: ArrayCtor<T>): Column<T> {
        if (c['@array']) return c;
        if (!c.isDefined) return Undefined(c.rowCount, c.schema) as any as Column<T>;
        return arrayColumn({ array: c.toArray({ array }), schema: c.schema, valueKind: c.valueKind });
    }
}

export default Column;

function createFirstIndexMapOfColumn<T>(c: Column<T>): Map<T, number> {
    const map = new Map<T, number>();
    for (let i = 0, _i = c.rowCount; i < _i; i++) {
        const v = c.value(i);
        if (!map.has(v)) return map.set(c.value(i), i);
    }
    return map;
}

function constColumn<T extends Column.Schema>(v: T['T'], rowCount: number, schema: T, valueKind: Column.ValueKind): Column<T['T']> {
    const value: Column<T['T']>['value'] = row => v;
    return {
        schema: schema,
        '@array': void 0,
        isDefined: valueKind === Column.ValueKind.Present,
        rowCount,
        value,
        valueKind: row => valueKind,
        toArray: params => {
            const { array } = ColumnHelpers.createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v;
            return array;
        },
        areValuesEqual: (rowA, rowB) => true
    }
}

function lambdaColumn<T extends Column.Schema>({ value, valueKind, rowCount, schema }: Column.LambdaSpec<T>): Column<T['T']> {
    return {
        schema: schema,
        '@array': void 0,
        isDefined: true,
        rowCount,
        value,
        valueKind: valueKind ? valueKind : row => Column.ValueKind.Present,
        toArray: params => {
            const { array, start } = ColumnHelpers.createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = value(i + start);
            return array;
        },
        areValuesEqual: (rowA, rowB) => value(rowA) === value(rowB)
    }
}

function arrayColumn<T extends Column.Schema>({ array, schema, valueKind }: Column.ArraySpec<T>): Column<T['T']> {
    const rowCount = array.length;
    const value: Column<T['T']>['value'] = schema.valueType === 'str'
        ? row => { const v = array[row]; return typeof v === 'string' ? v : '' + v; }
        : row => array[row];

    const isTyped = ColumnHelpers.isTypedArray(array);
    return {
        schema: schema,
        '@array': array,
        isDefined: true,
        rowCount,
        value,
        valueKind: valueKind ? valueKind : row => Column.ValueKind.Present,
        toArray: schema.valueType === 'str'
            ? params => {
                const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
                const ret = new (params && typeof params.array !== 'undefined' ? params.array : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) {
                    const v = array[start + i];
                    ret[i] = typeof v === 'string' ? v : '' + v;
                }
                return ret;
            }
            : isTyped
            ? params => ColumnHelpers.typedArrayWindow(array, params) as any as ReadonlyArray<T>
            : params => {
                const { start, end } = ColumnHelpers.getArrayBounds(rowCount, params);
                if (start === 0 && end === array.length) return array as ReadonlyArray<T['T']>;
                const ret = new (params && typeof params.array !== 'undefined' ? params.array : (array as any).constructor)(end - start) as any;
                for (let i = 0, _i = end - start; i < _i; i++) ret[i] = array[start + i];
                return ret;
            },
        areValuesEqual: (rowA, rowB) => array[rowA] === array[rowB]
    }
}

function windowColumn<T>(column: Column<T>, start: number, end: number) {
    if (!column.isDefined) return Column.Undefined(end - start, column.schema);
    if (!!column['@array'] && ColumnHelpers.isTypedArray(column['@array'])) return windowTyped(column, start, end);
    return windowFull(column, start, end);
}

function windowTyped<T>(c: Column<T>, start: number, end: number): Column<T> {
    const array = ColumnHelpers.typedArrayWindow(c['@array'], { start, end });
    return arrayColumn({ array, schema: c.schema, valueKind: c.valueKind }) as any;
}

function windowFull<T>(c: Column<T>, start: number, end: number): Column<T> {
    const v = c.value, vk = c.valueKind, ave = c.areValuesEqual;
    const value: Column<T>['value'] = start === 0 ? v : row => v(row + start);
    const rowCount = end - start;
    return {
        schema: c.schema,
        '@array': void 0,
        isDefined: c.isDefined,
        rowCount,
        value,
        valueKind: start === 0 ? vk : row => vk(row + start),
        toArray: params => {
            const { array } = ColumnHelpers.createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v(i + start);
            return array;
        },
        areValuesEqual: start === 0 ? ave : (rowA, rowB) => ave(rowA + start, rowB + start)
    };
}

function isIdentity(map: ArrayLike<number>, rowCount: number) {
    if (map.length !== rowCount) return false;
    for (let i = 0, _i = map.length; i < _i; i++) {
        if (map[i] !== i) return false;
    }
    return true;
}

function columnView<T>(c: Column<T>, map: ArrayLike<number>, checkIdentity: boolean): Column<T> {
    if (!c.isDefined) return c;
    if (checkIdentity && isIdentity(map, c.rowCount)) return c;
    if (!!c['@array']) return arrayView(c, map);
    return viewFull(c, map);
}

function arrayView<T>(c: Column<T>, map: ArrayLike<number>): Column<T> {
    const array = c['@array']!;
    const ret = new (array as any).constructor(map.length);
    for (let i = 0, _i = map.length; i < _i; i++) ret[i] = array[map[i]];
    return arrayColumn({ array: ret, schema: c.schema, valueKind: c.valueKind });
}

function viewFull<T>(c: Column<T>, map: ArrayLike<number>): Column<T> {
    const v = c.value, vk = c.valueKind, ave = c.areValuesEqual;
    const value: Column<T>['value'] = row => v(map[row]);
    const rowCount = map.length;
    return {
        schema: c.schema,
        '@array': void 0,
        isDefined: c.isDefined,
        rowCount,
        value,
        valueKind: row => vk(map[row]),
        toArray: params => {
            const { array } = ColumnHelpers.createArray(rowCount, params);
            for (let i = 0, _i = array.length; i < _i; i++) array[i] = v(map[i]);
            return array;
        },
        areValuesEqual: (rowA, rowB) => ave(map[rowA], map[rowB])
    };
}

function mapToArrayImpl<T, S>(c: Column<T>, f: (v: T) => S, ctor: Column.ArrayCtor<S>): ArrayLike<S> {
    const ret = new ctor(c.rowCount) as any;
    for (let i = 0, _i = c.rowCount; i < _i; i++) ret[i] = f(c.value(i));
    return ret;
}

function areColumnsEqual(a: Column<any>, b: Column<any>) {
    if (a.rowCount !== b.rowCount || a.isDefined !== b.isDefined || a.schema.valueType !== b.schema.valueType) return false;
    if (!!a['@array'] && !!b['@array']) return areArraysEqual(a, b);
    return areValuesEqual(a, b);
}

function areArraysEqual(a: Column<any>, b: Column<any>) {
    const xs = a['@array']!, ys = b['@array']!;
    for (let i = 0, _i = a.rowCount; i < _i; i++) {
        if (xs[i] !== ys[i]) return false;
    }
    return true;
}

function areValuesEqual(a: Column<any>, b: Column<any>) {
    const va = a.value, vb = b.value;
    for (let i = 0, _i = a.rowCount; i < _i; i++) {
        if (va(i) !== vb(i)) return false;
    }
    return true;
}

function columnIndicesOf<T>(c: Column<T>, test: (e: T) => boolean) {
    const ret = [], v = c.value;
    for (let i = 0, _i = c.rowCount; i < _i; i++) {
        if (test(v(i))) ret[ret.length] = i;
    }
    return ret;
}