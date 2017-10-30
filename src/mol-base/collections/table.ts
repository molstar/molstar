/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from './column'
import { sortArray } from './sort'

type Table<Schema extends Table.Schema> = { readonly _rowCount: number, readonly _columns: ReadonlyArray<string> } & Table.Columns<Schema>

/** An immutable table */
namespace Table {
    export type Schema = { [field: string]: Column.Type }
    export type Columns<S extends Schema> = { [C in keyof S]: Column<S[C]['T']> }
    export type Row<S extends Schema> = { [C in keyof S]: S[C]['T'] }
    export type Arrays<S extends Schema> = { [C in keyof S]: ArrayLike<S[C]['T']> }

    export function pickColumns<S extends Schema, T extends S>(schema: S, table: Table<T>): Table<S> {
        const ret = Object.create(null);
        const keys = Object.keys(schema);
        ret._rowCount = table._rowCount;
        ret._columns = keys;
        for (const k of keys) ret[k] = table[k];
        return ret;
    }

    export function ofColumns<S extends Schema, R extends Table<S> = Table<S>>(columns: Columns<S>): R {
        const _columns = Object.keys(columns);
        const _rowCount = columns[_columns[0]].rowCount;
        return { _rowCount, _columns, ...(columns as any) };
    }

    export function ofRows<S extends Schema, R extends Table<S> = Table<S>>(schema: Schema, rows: ArrayLike<Row<S>>): R {
        const ret = Object.create(null);
        const rowCount = rows.length;
        const columns = Object.keys(schema);
        ret._rowCount = rowCount;
        ret._columns = columns;
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
        for (const k of Object.keys(schema)) {
            (ret as any)[k] = Column.ofArray({ array: arrays[k], type: schema[k] })
        }
        return ret as R;
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
        for (const c of table._columns) {
            ret[c] = Column.permutation((table as any)[c], indices, false);
        }
        return ret;
    }
}

export default Table