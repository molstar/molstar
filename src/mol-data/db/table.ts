/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from './column';
import { sortArray } from '../util/sort';
import { StringBuilder } from '../../mol-util';

/** A collection of columns */
type Table<Schema extends Table.Schema = any> = {
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
    export type PartialColumns<S extends Schema> = { [C in keyof S]?: Column<S[C]['T']> }
    export type PartialTable<S extends Table.Schema> = { readonly _rowCount: number, readonly _columns: ReadonlyArray<string> } & PartialColumns<S>

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

    export function ofPartialColumns<S extends Schema, R extends Table<S> = Table<S>>(schema: S, partialColumns: PartialColumns<S>, rowCount: number): R {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = rowCount;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            if (k in partialColumns) ret[k] = partialColumns[k];
            else ret[k] = Column.Undefined(rowCount, schema[k]);
        }
        return ret;
    }

    export function ofUndefinedColumns<S extends Schema, R extends Table<S> = Table<S>>(schema: S, rowCount: number): R {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = rowCount;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            ret[k] = Column.Undefined(rowCount, schema[k]);
        }
        return ret;
    }

    export function ofRows<S extends Schema, R extends Table<S> = Table<S>>(schema: S, rows: ArrayLike<Partial<Row<S>>>): R {
        const ret = Object.create(null);
        const rowCount = rows.length;
        const columns = Object.keys(schema);
        ret._rowCount = rowCount;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            (ret as any)[k] = Column.ofLambda({
                rowCount,
                schema: schema[k],
                value: r => rows[r][k],
                valueKind: r => typeof rows[r][k] === 'undefined' ? Column.ValueKind.NotPresent : Column.ValueKind.Present
            });
        }
        return ret as R;
    }

    export function ofArrays<S extends Schema, R extends Table<S> = Table<S>>(schema: S, arrays: Partial<Arrays<S>>): R {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = 0;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            if (typeof arrays[k] !== 'undefined') {
                (ret as any)[k] = Column.ofArray({ array: arrays[k]!, schema: schema[k] });
                ret._rowCount = arrays[k]?.length;
            } else {
                (ret as any)[k] = Column.Undefined(ret._rowCount, schema[k]);
            }
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

    export function pick<S extends R, R extends Schema>(table: Table<S>, schema: R, test: (i: number) => boolean) {
        const _view: number[] = [];
        for (let i = 0, il = table._rowCount; i < il; ++i) {
            if (test(i)) _view.push(i);
        }
        return view(table, schema, _view);
    }

    export function window<S extends R, R extends Schema>(table: Table<S>, schema: R, start: number, end: number) {
        if (start === 0 && end === table._rowCount) return table;
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = end - start;
        ret._columns = columns;
        ret._schema = schema;
        for (const k of columns) {
            (ret as any)[k] = Column.window(table[k], start, end);
        }
        return ret as Table<R>;
    }

    export function concat<S extends R, R extends Schema>(tables: Table<S>[], schema: R) {
        const ret = Object.create(null);
        const columns = Object.keys(schema);
        ret._rowCount = 0;
        for (const table of tables) {
            ret._rowCount += table._rowCount;
        }
        const arrays: any = {};
        for (const column of columns) {
            arrays[column] = new Array(ret._rowCount);
        }
        ret._columns = columns;
        ret._schema = schema;
        let offset = 0;
        for (const table of tables) {
            for (const k of columns) {
                Column.copyToArray(table[k], arrays[k], offset);
            }
            offset += table._rowCount;
        }
        for (const k of columns) {
            ret[k] = Column.ofArray({ array: arrays[k], schema: schema[k] });
        }
        return ret as Table<R>;
    }

    export function columnToArray<S extends Schema>(table: Table<S>, name: keyof S, array?: Column.ArrayCtor<any>) {
        (table as Columns<S>)[name] = Column.asArrayColumn((table as Columns<S>)[name], array);
    }

    /** Sort and return a new table */
    export function sort<T extends Table>(table: T, cmp: (i: number, j: number) => number) {
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

    export function areEqual<T extends Table<any>>(a: T, b: T) {
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

    /** Allocate a new object with the given row values. */
    export function getRow<S extends Schema>(table: Table<S>, index: number) {
        const row: Row<S> = Object.create(null);
        const { _columns: cols } = table;
        for (let i = 0; i < cols.length; i++) {
            const c = cols[i] as keyof S;
            row[c] = table[c].value(index);
        }
        return row;
    }

    /** Pick the first row for which `test` evaluates to true */
    export function pickRow<S extends Schema>(table: Table<S>, test: (i: number) => boolean) {
        for (let i = 0, il = table._rowCount; i < il; ++i) {
            if (test(i)) return getRow(table, i);
        }
    }

    export function getRows<S extends Schema>(table: Table<S>) {
        const ret: Row<S>[] = [];
        const { _rowCount: c } = table;
        for (let i = 0; i < c; i++) {
            ret[i] = getRow(table, i);
        }
        return ret;
    }

    export function toArrays<S extends Schema>(table: Table<S>) {
        const arrays: { [k: string]: ArrayLike<any> } = {};
        const { _columns } = table;
        for (let i = 0; i < _columns.length; i++) {
            const c = _columns[i];
            arrays[c] = table[c].toArray();
        }
        return arrays as { [k in keyof S]: ArrayLike<S[k]['T']> };
    }

    export function formatToString<S extends Schema>(table: Table<S>) {
        const sb = StringBuilder.create();

        const { _columns: cols, _rowCount } = table;

        let headerLength = 1;
        StringBuilder.write(sb, '|');
        for (let i = 0; i < cols.length; i++) {
            StringBuilder.write(sb, cols[i]);
            StringBuilder.write(sb, '|');
            headerLength += cols[i].length + 1;
        }
        StringBuilder.newline(sb);
        StringBuilder.write(sb, new Array(headerLength + 1).join('-'));
        StringBuilder.newline(sb);

        for (let r = 0; r < _rowCount; r++) {
            StringBuilder.write(sb, '|');
            for (let i = 0; i < cols.length; i++) {
                const c = table[cols[i]];
                if (c.valueKind(r) === Column.ValueKind.Present) {
                    StringBuilder.write(sb, c.value(r));
                    StringBuilder.write(sb, '|');
                } else {
                    StringBuilder.write(sb, '.|');
                }
            }
            StringBuilder.newline(sb);
        }
        return StringBuilder.getString(sb);
    }
}

export default Table;