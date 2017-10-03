/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data-model'
import * as Column from '../common/column'
import StringPool from '../../utils/short-string-pool'

/**
 * A schema defines the shape of categories and fields.
 *
 * @example:
 * const atom_site = {
 *   '@alias': '_atom_site',
 *   label_atom_id: Field.str(),
 *   Cartn_x: Field.float(),
 *   Cartn_y: Field.float(),
 *   Cartn_z: Field.float(),
 * }
 *
 * const mmCIF = { atom_site };
 */

//////////////////////////////////////////////

export function apply<Schema extends Block.Schema>(schema: Schema, block: Data.Block): Block.Instance<Schema> {
    return createBlock(schema, block) as Block.Instance<Schema>;
}

export type Block<Categories> = Categories & {
    readonly _header?: string,
    /** For accessing 'non-standard' categories */
    _getCategory(name: string): Data.Category | undefined
}

export namespace Block {
    export type Schema = { [category: string]: Category.Schema }
    export type Instance<T extends Schema> = Block<{ [C in keyof T]: Category.Instance<T[C]> }>
}

export type Category<Fields> = Fields & {
    readonly _rowCount: number,
    readonly _isDefined: boolean,
    /** For accessing 'non-standard' fields */
    _getField(name: string): Data.Field | undefined
}

export namespace Category {
    export type Schema = { [field: string]: Field.Schema<any> }
    export type Instance<T extends Schema> = Category<{ [F in keyof T]: Column.Column<T[F]['type']> }>
}

export namespace Field {
    export interface Schema<T> { type: T, ctor: (field: Data.Field, category: Data.Category, key: string) => Column.Column<T>, undefinedField: (c: number) => Data.Field, alias?: string };
    export interface Spec { undefinedField?: (c: number) => Data.Field, alias?: string }

    export function alias(name: string): Schema<any> { return { alias: name } as any; }
    export function pooledStr(spec?: Spec) { return createSchema(spec, PooledStr); }
    export function str(spec?: Spec) { return createSchema(spec, Str); }
    export function int(spec?: Spec) { return createSchema(spec, Int); }
    export function float(spec?: Spec) { return createSchema(spec, Float); }
    export function vector(rows: number, spec?: Spec) { return createSchema(spec, Vector(rows)); }
    export function matrix(rows: number, cols: number, spec?: Spec) { return createSchema(spec, Matrix(rows, cols)); }

    function create<T>(field: Data.Field, value: (row: number) => T, toArray: Column.Column<T>['toArray']): Column.Column<T> {
        const presence = field.presence;
        return {
            isDefined: field.isDefined,
            rowCount: field.rowCount,
            value,
            isValueDefined: row => presence(row) === Data.ValuePresence.Present,
            areValuesEqual: field.areValuesEqual,
            toArray
        };
    }

    function PooledStr(field: Data.Field) {
        const pool = StringPool.create();
        const value = (row: number) => StringPool.get(pool, field.str(row));
        const array = (params?: Column.ToArrayParams) => Column.createAndFillArray(field.rowCount, value, params);
        return create<string>(field, value, array);
    }
    function Str(field: Data.Field) { return create(field, field.str, field.toStringArray); }
    function Int(field: Data.Field) { return create(field, field.int, field.toIntArray); }
    function Float(field: Data.Field) { return create(field, field.float, field.toFloatArray); }

    function Vector(rows: number) {
        return function(field: Data.Field, category: Data.Category, key: string) {
            const value = (row: number) => Data.getVector(category, key, rows, row);
            return create(field, value, params => Column.createAndFillArray(field.rowCount, value, params));
        }
    }

    function Matrix(rows: number, cols: number) {
        return function(field: Data.Field, category: Data.Category, key: string) {
            const value = (row: number) => Data.getMatrix(category, key, rows, cols, row);
            return create(field, value, params => Column.createAndFillArray(field.rowCount, value, params));
        }
    }

    function createSchema<T>(spec: Spec | undefined, ctor: (field: Data.Field, category: Data.Category, key: string) => Column.Column<T>): Schema<T> {
        return { type: 0 as any, ctor, undefinedField: (spec && spec.undefinedField) || Data.DefaultUndefinedField, alias: spec && spec.alias };
    }
}

class _Block implements Block<any> { // tslint:disable-line:class-name
    header = this._block.header;
    getCategory(name: string) { return this._block.categories[name]; }
    constructor(private _block: Data.Block, schema: Block.Schema) {
        for (const k of Object.keys(schema)) {
            Object.defineProperty(this, k, { value: createCategory(k, schema[k], _block), enumerable: true, writable: false, configurable: false });
        }
    }
}

class _Category implements Category<any> { // tslint:disable-line:class-name
    _rowCount = this._category.rowCount;
    _getField(name: string) { return this._category.getField(name); }
    constructor(private _category: Data.Category, schema: Category.Schema, public _isDefined: boolean) {
        const fieldKeys = Object.keys(schema).filter(k => k !== '@alias');
        const cache = Object.create(null);
        for (const k of fieldKeys) {
            const s = schema[k];
            Object.defineProperty(this, k, {
                get: function() {
                    if (cache[k]) return cache[k];
                    const name = s.alias || k;
                    const field = _category.getField(name) || s.undefinedField(_category.rowCount);
                    cache[k] = s.ctor(field, _category, name);
                    return cache[k];
                },
                enumerable: true,
                configurable: false
            });
        }
    }
}

function createBlock(schema: Block.Schema, block: Data.Block): any {
    return new _Block(block, schema);
}

function createCategory(key: string, schema: Category.Schema, block: Data.Block) {
    const alias = (schema['@alias'] && schema['@alias'].alias) || key;
    const name = alias[0] === '_' ? alias : '_' + alias;
    const cat = block.categories[name];
    return new _Category(cat || Data.Category.Empty, schema, !!cat);
}