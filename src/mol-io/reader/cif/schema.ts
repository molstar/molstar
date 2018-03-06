/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Table, Column, ColumnHelpers } from 'mol-data/db'
import { Tensor } from 'mol-math/linear-algebra'
import * as Data from './data-model'

export function toDatabase<Schema extends Database.Schema, Frame extends Database<Schema> = Database<Schema>>(schema: Schema, frame: Data.Frame): Frame {
    return createDatabase(schema, frame) as Frame;
}

export function toTable<Schema extends Table.Schema, R extends Table<Schema> = Table<Schema>>(schema: Schema, category: Data.Category): R {
    return new CategoryTable(category, schema, true) as any;
}

type ColumnCtor = (field: Data.Field, category: Data.Category, key: string) => Column<any>

function getColumnCtor(t: Column.Schema): ColumnCtor {
    switch (t.valueType) {
        case 'str': return (f, c, k) => createColumn(t, f, f.str, f.toStringArray);
        case 'int': return (f, c, k) => createColumn(t, f, f.int, f.toIntArray);
        case 'float': return (f, c, k) => createColumn(t, f, f.float, f.toFloatArray);
        case 'list': return (f, c, k) => createColumn(t, f, f.list, f.toListArray);
        case 'tensor': throw new Error(`Use createTensorColumn instead.`);
    }
}

function createColumn<T>(schema: Column.Schema, field: Data.Field, value: (row: number) => T, toArray: Column<T>['toArray']): Column<T> {
    return {
        schema,
        '@array': field['@array'],
        isDefined: field.isDefined,
        rowCount: field.rowCount,
        value,
        valueKind: field.valueKind,
        areValuesEqual: field.areValuesEqual,
        toArray
    };
}

function createTensorColumn(schema: Column.Schema.Tensor, category: Data.Category, key: string): Column<Tensor> {
    const space = schema.space;
    const first = category.getField(`${key}[1]`) || Column.Undefined(category.rowCount, schema);
    const value = (row: number) => Data.getTensor(category, key, space, row);
    const toArray: Column<Tensor>['toArray'] = params => ColumnHelpers.createAndFillArray(category.rowCount, value, params)

    return {
        schema,
        '@array': void 0,
        isDefined: first.isDefined,
        rowCount: category.rowCount,
        value,
        valueKind: first.valueKind,
        areValuesEqual: (rowA, rowB) => Tensor.areEqualExact(value(rowA), value(rowB)),
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
            Object.defineProperty(this, k, {
                get: function() {
                    if (cache[k]) return cache[k];
                    const fType = schema[k];
                    if (fType.valueType === 'tensor') {
                        cache[k] = createTensorColumn(fType, category, k);
                    } else {
                        const ctor = getColumnCtor(fType);
                        const field = category.getField(k);
                        cache[k] = !!field ? ctor(field, category, k) : Column.Undefined(category.rowCount, fType);
                    }
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