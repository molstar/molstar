/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DatabaseCollection, Database, Table, Column, ColumnHelpers } from '../../../mol-data/db'
import { Tensor } from '../../../mol-math/linear-algebra'
import { arrayEqual } from '../../../mol-util'
import * as Data from './data-model'

export function toDatabaseCollection<Schema extends Database.Schema>(schema: Schema, file: Data.CifFile): DatabaseCollection<Schema> {
    const dbc: DatabaseCollection<Schema> = {}
    for (const data of file.blocks) {
        dbc[data.header] = toDatabase(schema, data)
    }
    return dbc;
}

export function toDatabase<Schema extends Database.Schema, Frame extends Database<Schema> = Database<Schema>>(schema: Schema, frame: Data.CifFrame): Frame {
    return createDatabase(schema, frame) as Frame;
}

export function toTable<Schema extends Table.Schema, R extends Table<Schema> = Table<Schema>>(schema: Schema, category: Data.CifCategory): R {
    return new CategoryTable(category, schema, true) as any;
}

type ColumnCtor = (field: Data.CifField, category: Data.CifCategory, key: string) => Column<any>

function getColumnCtor(t: Column.Schema): ColumnCtor {
    switch (t.valueType) {
        case 'str': return (f, c, k) => createColumn(t, f, f.str, f.toStringArray);
        case 'int': return (f, c, k) => createColumn(t, f, f.int, f.toIntArray);
        case 'float': return (f, c, k) => createColumn(t, f, f.float, f.toFloatArray);
        case 'list': throw new Error('Use createListColumn instead.');
        case 'tensor': throw new Error('Use createTensorColumn instead.');
    }
}

function createColumn<T>(schema: Column.Schema, field: Data.CifField, value: (row: number) => T, toArray: Column<T>['toArray']): Column<T> {
    return {
        schema,
        __array: field.__array,
        isDefined: field.isDefined,
        rowCount: field.rowCount,
        value,
        valueKind: field.valueKind,
        areValuesEqual: field.areValuesEqual,
        toArray
    };
}

function createListColumn<T extends number|string>(schema: Column.Schema.List<T>, category: Data.CifCategory, key: string): Column<(number|string)[]> {
    const separator = schema.separator;
    const itemParse = schema.itemParse;

    const f = category.getField(key);
    const value = f ? (row: number) => f.str(row).split(separator).map(x => itemParse(x.trim())).filter(x => !!x) : (row: number) => []
    const toArray: Column<T[]>['toArray'] = params => ColumnHelpers.createAndFillArray(category.rowCount, value, params)

    return {
        schema,
        __array: void 0,
        isDefined: !!f,
        rowCount: category.rowCount,
        value,
        valueKind: f ? f.valueKind : () => Column.ValueKind.NotPresent,
        areValuesEqual: (rowA, rowB) => arrayEqual(value(rowA), value(rowB)),
        toArray
    };
}

function createTensorColumn(schema: Column.Schema.Tensor, category: Data.CifCategory, key: string): Column<Tensor.Data> {
    const space = schema.space;
    const zeroOffset = category.fieldNames.indexOf(`${key}[0]`) >= 0;
    const fst = zeroOffset ? 0 : 1;

    let firstFieldName: string;
    switch (space.rank) {
        case 1: firstFieldName = `${key}[${fst}]`; break;
        case 2: firstFieldName = `${key}[${fst}][${fst}]`; break;
        case 3: firstFieldName = `${key}[${fst}][${fst}][${fst}]`; break;
        default: throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.');
    }
    const first = category.getField(firstFieldName) || Column.Undefined(category.rowCount, schema);
    const value = (row: number) => Data.getTensor(category, key, space, row, zeroOffset);
    const toArray: Column<Tensor.Data>['toArray'] = params => ColumnHelpers.createAndFillArray(category.rowCount, value, params)

    return {
        schema,
        __array: void 0,
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

    constructor(category: Data.CifCategory, schema: Table.Schema, public _isDefined: boolean) {
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
                    if (fType.valueType === 'list') {
                        cache[k] = createListColumn(fType, category, k);
                    } else if (fType.valueType === 'tensor') {
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

function createDatabase(schema: Database.Schema, frame: Data.CifFrame): Database<any> {
    const tables = Object.create(null);
    for (const k of Object.keys(schema)) {
        tables[k] = createTable(k, (schema as any)[k], frame);
    }
    return Database.ofTables(frame.header, schema, tables);
}

function createTable(key: string, schema: Table.Schema, frame: Data.CifFrame) {
    const cat = frame.categories[key];
    return new CategoryTable(cat || Data.CifCategory.empty(key), schema, !!cat);
}