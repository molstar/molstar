/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DatabaseCollection, Database, Table, Column, ColumnHelpers } from '../../../mol-data/db';
import { Tensor } from '../../../mol-math/linear-algebra';
import { arrayEqual } from '../../../mol-util';
import * as Data from './data-model';

export namespace FieldPath {
    export function canonical(path: string) {
        return path.replace('.', '_').replace(/\[/, '_').replace(/(\[|\])/g, '');
    }

    export function equal(pathA: string, pathB: string) {
        return canonical(pathA) === canonical(pathB);
    }

    export function create(category: string, field: string, asCanonical = false) {
        const p = `${category}${field ? `.${field}` : ''}`;
        return asCanonical ? canonical(p) : p;
    }
}

export function toDatabaseCollection<Schema extends Database.Schema>(schema: Schema, file: Data.CifFile, aliases?: Data.CifAliases): DatabaseCollection<Schema> {
    const dbc: DatabaseCollection<Schema> = {};
    for (const data of file.blocks) {
        dbc[data.header] = toDatabase(schema, data, aliases);
    }
    return dbc;
}

export function toDatabase<Schema extends Database.Schema, Frame extends Database<Schema> = Database<Schema>>(schema: Schema, frame: Data.CifFrame, aliases?: Data.CifAliases): Frame {
    return createDatabase(schema, frame, aliases) as Frame;
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

function createListColumn<T extends number | string>(schema: Column.Schema.List<T>, category: Data.CifCategory, key: string): Column<(number | string)[]> {
    const separator = schema.separator;
    const itemParse = schema.itemParse;

    const f = category.getField(key);
    const value = f ? (row: number) => f.str(row).split(separator).map(x => itemParse(x.trim())).filter(x => !!x) : (row: number) => [];
    const toArray: Column<T[]>['toArray'] = params => ColumnHelpers.createAndFillArray(category.rowCount, value, params);

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
    const zeroOffset = (
        category.fieldNames.includes(`${key}[0]`) ||
        category.fieldNames.includes(`${key}[0][0]`) ||
        category.fieldNames.includes(`${key}[0][0][0]`)
    );
    const fst = zeroOffset ? 0 : 1;
    const namingVariant = (
        category.fieldNames.includes(`${key}_1`) ||
        category.fieldNames.includes(`${key}_11`) ||
        category.fieldNames.includes(`${key}_111`)
    ) ? 'underscore' : 'brackets';

    const getName = Data.tensorFieldNameGetter(key, space.rank, zeroOffset, namingVariant);
    const first = category.getField(getName(fst, fst, fst)) || Column.Undefined(category.rowCount, schema);
    const value = (row: number) => Data.getTensor(category, space, row, getName);
    const toArray: Column<Tensor.Data>['toArray'] = params => ColumnHelpers.createAndFillArray(category.rowCount, value, params);

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

function createDatabase(schema: Database.Schema, frame: Data.CifFrame, aliases?: Data.CifAliases): Database<any> {
    const tables = Object.create(null);
    for (const k of Object.keys(schema)) {
        tables[k] = createTable(k, schema[k], frame, aliases);
    }
    return Database.ofTables(frame.header, schema, tables);
}

type FlatFrame = { [k: string]: Data.CifField }

function flattenFrame(frame: Data.CifFrame): FlatFrame {
    const flatFrame = Object.create(null);
    for (const c of Object.keys(frame.categories)) {
        for (const f of frame.categories[c].fieldNames) {
            const p =  FieldPath.create(c, f, true);
            flatFrame[p] = frame.categories[c].getField(f);
        }
    }
    return flatFrame;
}

function getField(field: string, category: string, flatFrame: FlatFrame, aliases?: Data.CifAliases) {
    const path = FieldPath.create(category, field);
    const canonicalPath = FieldPath.canonical(path);
    if (canonicalPath in flatFrame) return flatFrame[canonicalPath];
    if (aliases && path in aliases) {
        for (const aliased of aliases[path]) {
            const canonicalAliased = FieldPath.canonical(aliased);
            if (canonicalAliased in flatFrame) return flatFrame[canonicalAliased];
        }
    }
}

function createTable(key: string, schema: Table.Schema, frame: Data.CifFrame, aliases?: Data.CifAliases) {
    let cat = frame.categories[key];
    if (aliases) {
        const flatFrame = flattenFrame(frame);
        const fields: { [k: string]: Data.CifField } = Object.create(null);
        const fieldNames: string[] = [];
        let rowCount = 0;
        for (const k of Object.keys(schema)) {
            const field = getField(k, key, flatFrame, aliases);
            if (field) {
                fields[k] = field;
                fieldNames.push(k);
                rowCount = field.rowCount;
            }
        }
        cat = {
            rowCount,
            name: key,
            fieldNames: [...fieldNames],
            getField(name: string) {
                return fields[name];
            }
        };
    }
    return new CategoryTable(cat || Data.CifCategory.empty(key), schema, !!cat);
}