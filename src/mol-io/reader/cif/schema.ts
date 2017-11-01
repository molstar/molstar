/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data-model'
import Column, { ColumnHelpers } from '../../../mol-base/collections/column'
import Table from '../../../mol-base/collections/table'

export function toTypedFrame<Schema extends FrameSchema, Frame extends TypedFrame<Schema> = TypedFrame<Schema>>(schema: Schema, frame: Data.Frame): Frame {
    return createTypedFrame(schema, frame) as Frame;
}

export function toTable<Schema extends Table.Schema, R extends Table<Schema> = Table<Schema>>(schema: Schema, category: Data.Category): R {
    return new _TypedCategory(category, schema, true) as any;
}

export const Types = Column.Type

export type FrameSchema = { [category: string]: Table.Schema }
export type TypedFrame<Schema extends FrameSchema> = {
    readonly _header?: string,
    readonly _frame: Data.Frame
} & { [C in keyof Schema]: Table<Schema[C]> }

type ColumnCtor = (field: Data.Field, category: Data.Category, key: string) => Column<any>

function getColumnCtor(t: Column.Type): ColumnCtor {
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

class _TypedFrame implements TypedFrame<any> { // tslint:disable-line:class-name
    header = this._frame.header;
    [k: string]: any;
    constructor(public _frame: Data.Frame, schema: FrameSchema) {
        for (const k of Object.keys(schema)) {
            Object.defineProperty(this, k, { value: createTypedCategory(k, schema[k], _frame), enumerable: true, writable: false, configurable: false });
        }
    }
}

class _TypedCategory implements Table<any> { // tslint:disable-line:class-name
    _rowCount = this._category.rowCount;
    _columns: ReadonlyArray<string>;
    [k: string]: any;
    constructor(public _category: Data.Category, schema: Table.Schema, public _isDefined: boolean) {
        const fieldKeys = Object.keys(schema);
        this._columns = fieldKeys;
        const cache = Object.create(null);
        for (const k of fieldKeys) {
            const cType = schema[k];
            const ctor = getColumnCtor(cType);
            Object.defineProperty(this, k, {
                get: function() {
                    if (cache[k]) return cache[k];
                    const field = _category.getField(k);
                    cache[k] = !!field ? ctor(field, _category, k) : Column.Undefined(_category.rowCount, cType);
                    return cache[k];
                },
                enumerable: true,
                configurable: false
            });
        }
    }
}

function createTypedFrame(schema: FrameSchema, frame: Data.Frame): any {
    return new _TypedFrame(frame, schema);
}

function createTypedCategory(key: string, schema: Table.Schema, frame: Data.Frame) {
    const cat = frame.categories[key[0] === '_' ? key : '_' + key];
    return new _TypedCategory(cat || Data.Category.Empty, schema, !!cat);
}