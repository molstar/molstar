/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from './column'
import { sortArray } from '../sort'

/** A collection of columns */
type Table<Schema extends Table.Schema> = {
    readonly _rowCount: number,
    readonly _columns: ReadonlyArray<string>,
    readonly _schema: Schema
} & Table.Columns<Schema>

/** An immutable table */
namespace Table {
    export type Schema = { [field: string]: Column.Schema }
    export type Columns<S extends Schema> = { [C in keyof S]: Column<S[C]['T']> }
    export type Row<S extends Schema> = { [C in keyof S]: S[C]['T'] }
    export type Arrays<S extends Schema> = { [C in keyof S]: ArrayLike<S[C]['T']> }
    export type PartialTable<S extends Table.Schema> = { readonly _rowCount: number, readonly _columns: ReadonlyArray<string> } & { [C in keyof S]?: Column<S[C]['T']> }

    export function is(t: any): t is Table<any> {
        return t && typeof t._rowCount === 'number' && !!t._columns && !!t._schema;
    }

    export function pickColumns<S extends Schema>(schema: S, table: PartialTable<S>, guard: Partial<Columns<S>> = {}): Table<S> {
        const ret = Object.create(null);
        const keys = Object.keys(schema);
        ret._rowCount = table._rowCount;
        ret._columns = keys;
        ret._schema = schema;
        for (const k of keys) {
            if (!!table[k]) ret[k] = table[k];
            else if (!!guard[k]) ret[k] = guard[k];
            else throw Error(`Cannot find column '${k}'.`);
        }
        return ret;
    }

    export function ofColumns<S extends Schema, R extends Table<S> = Table<S>>(schema: S, columns: Columns<S>): R {
        const _columns = Object.keys(columns);
        const _rowCount = columns[_columns[0]].rowCount;
        return { _rowCount, _columns, _schema: schema, ...(columns as any) };
    }

    export function ofRows<S extends Schema, R extends Table<S> = Table<S>>(schema: Schema, rows: ArrayLike<Row<S>>): R {
        const ret = Object.create(null);
        const rowCount = rows.length;
        const columns = Object.keys(schema);
        ret._rowCount = rowCount;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            (ret as any)[k] = Column.ofLambda({
                rowCount,
                type: schema[k],
                value: r => rows[r][k],
                valueKind: r => typeof rows[r][k] === 'undefined' ? Column.ValueKind.NotPresent : Column.ValueKind.Present
            })
        }
        return ret as R;
    }

    export function ofArrays<S extends Schema, R extends Table<S> = Table<S>>(schema: Schema, arrays: Arrays<S>): R {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = arrays[columns[0]].length;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            (ret as any)[k] = Column.ofArray({ array: arrays[k], type: schema[k] })
        }
        return ret as R;
    }

    export function view<S extends R, R extends Schema>(table: Table<S>, schema: R, view: ArrayLike<number>) {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = view.length;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            (ret as any)[k] = Column.view(table[k], view);
        }
        return ret as Table<R>;
    }

    export function columnToArray<S extends Schema>(table: Table<S>, name: keyof S, array?: Column.ArrayCtor<any>) {
        table[name] = Column.asArrayColumn(table[name], array);
    }

    /** Sort and return a new table */
    export function sort<T extends Table<S>, S extends Schema>(table: T, cmp: (i: number, j: number) => number) {
        const indices = new Int32Array(table._rowCount);
        for (let i = 0, _i = indices.length; i < _i; i++) indices[i] = i;
        sortArray(indices, (_, i, j) => cmp(i, j));

        let isIdentity = true;
        for (let i = 0, _i = indices.length; i < _i; i++) {
            if (indices[i] !== i) {
                isIdentity = false;
                break;
            }
        }
        if (isIdentity) return table;

        const ret = Object.create(null);
        ret._rowCount = table._rowCount;
        ret._columns = table._columns;
        ret._schema = table._schema;
        for (const c of table._columns) {
            ret[c] = Column.view((table as any)[c], indices, false);
        }
        return ret;
    }

    export function areEqual<T extends Table<Schema>>(a: T, b: T) {
        if (a._rowCount !== b._rowCount) return false;
        if (a._columns.length !== b._columns.length) return false;
        for (const c of a._columns) {
            if (!b[c]) return false;
        }

        for (const c of a._columns) {
            if (!Column.areEqual(a[c], b[c])) return false;
        }

        return true;
    }
}

export default Table