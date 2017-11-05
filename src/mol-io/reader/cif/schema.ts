/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data-model'
import { Database, Table, Column, ColumnHelpers } from 'mol-base/collections/database'

export function toDatabase<Schema extends Database.Schema, Frame extends Database<Schema> = Database<Schema>>(schema: Schema, frame: Data.Frame): Frame {
    return createDatabase(schema, frame) as Frame;
}

export function toTable<Schema extends Table.Schema, R extends Table<Schema> = Table<Schema>>(schema: Schema, category: Data.Category): R {
    return new CategoryTable(category, schema, true) as any;
}

type ColumnCtor = (field: Data.Field, category: Data.Category, key: string) => Column<any>

function getColumnCtor(t: Column.Schema): ColumnCtor {
    switch (t.kind) {
        case 'str': return (f, c, k) => createColumn(Column.Type.str, f, f.str, f.toStringArray);
        case 'int': return (f, c, k) => createColumn(Column.Type.int, f, f.int, f.toIntArray);
        case 'float': return (f, c, k) => createColumn(Column.Type.float, f, f.float, f.toFloatArray);
        case 'vector': return (f, c, k) => {
            const dim = t.dim;
            const value = (row: number) => Data.getVector(c, k, dim, row);
            return createColumn(t, f, value, params => ColumnHelpers.createAndFillArray(f.rowCount, value, params));
        }
        case 'matrix': return (f, c, k) => {
            const rows = t.rows, cols = t.cols;
            const value = (row: number) => Data.getMatrix(c, k, rows, cols, row);
            return createColumn(t, f, value, params => ColumnHelpers.createAndFillArray(f.rowCount, value, params));
        }
    }
}


function createColumn<T>(type: Column.Type, field: Data.Field, value: (row: number) => T, toArray: Column<T>['toArray']): Column<T> {
    return {
        '@type': type,
        '@array': field['@array'],
        isDefined: field.isDefined,
        rowCount: field.rowCount,
        value,
        valueKind: field.valueKind,
        areValuesEqual: field.areValuesEqual,
        toArray
    };
}

class CategoryTable implements Table<any> { // tslint:disable-line:class-name
    _rowCount: number;
    _columns: ReadonlyArray<string>;
    _schema: any;
    [k: string]: any;

    constructor(category: Data.Category, schema: Table.Schema, public _isDefined: boolean) {
        const fieldKeys = Object.keys(schema);
        this._rowCount = category.rowCount;
        this._columns = fieldKeys;
        this._schema = schema;
        const cache = Object.create(null);
        for (const k of fieldKeys) {
            const cType = schema[k];
            const ctor = getColumnCtor(cType);
            Object.defineProperty(this, k, {
                get: function() {
                    if (cache[k]) return cache[k];
                    const field = category.getField(k);
                    cache[k] = !!field ? ctor(field, category, k) : Column.Undefined(category.rowCount, cType);
                    return cache[k];
                },
                enumerable: true,
                configurable: false
            });
        }
    }
}

function createDatabase(schema: Database.Schema, frame: Data.Frame): Database<any> {
    const tables = Object.create(null);
    for (const k of Object.keys(schema)) {
        tables[k] = createTable(k, (schema as any)[k], frame);
    }
    return Database.ofTables(frame.header, schema, tables);
}

function createTable(key: string, schema: Table.Schema, frame: Data.Frame) {
    const cat = frame.categories[key];
    return new CategoryTable(cat || Data.Category.empty(key), schema, !!cat);
}