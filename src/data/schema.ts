/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data'

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
    /** For accessing 'non-standard' fields */
    _getField(name: string): Data.Field | undefined
}

export namespace Category {
    export type Schema = { '@alias'?: string } & { [field: string]: Field.Schema<any> }
    export type Instance<T extends Schema> = Category<{ [F in keyof T]: Field<T[F]['type']> }>
}

export interface Field<T> {
    readonly isDefined: boolean,
    value(row: number): T,
    presence(row: number): Data.ValuePresence,
    areValuesEqual(rowA: number, rowB: number): boolean,
    stringEquals(row: number, value: string | null): boolean,
    /** Converts the selected row range to an array. ctor might or might not be called depedning on the source data format. */
    toArray(startRow: number, endRowExclusive: number, ctor: (size: number) => Data.FieldArray): ReadonlyArray<T> | undefined
}

export namespace Field {
    export interface Schema<T> { type: T, ctor: (field: Data.Field) => Field<T>, undefinedField: Data.Field, alias?: string };
    export interface Spec { undefinedField?: Data.Field, alias?: string }

    export function str(spec?: Spec) { return createSchema(spec, Str); }
    export function int(spec?: Spec) { return createSchema(spec, Int); }
    export function float(spec?: Spec) { return createSchema(spec, Float); }

    function create<T>(field: Data.Field, value: (row: number) => T, toArray: Field<T>['toArray']): Field<T> {
        return { isDefined: field.isDefined, value, presence: field.presence, areValuesEqual: field.areValuesEqual, stringEquals: field.stringEquals, toArray };
    }

    function Str(field: Data.Field) { return create(field, field.str, field.toStringArray); }
    function Int(field: Data.Field) { return create(field, field.int, field.toNumberArray); }
    function Float(field: Data.Field) { return create(field, field.float, field.toNumberArray); }

    const DefaultUndefined: Data.Field = {
        isDefined: false,
        str: row => null,
        int: row => 0,
        float: row => 0,
        value: row => null,

        presence: row => Data.ValuePresence.NotSpecified,
        areValuesEqual: (rowA, rowB) => true,
        stringEquals: (row, value) => value === null,

        toStringArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = null; }
            return ret;
        },
        toNumberArray: (startRow, endRowExclusive, ctor) => new Uint8Array(endRowExclusive - startRow) as any
    };

    function createSchema<T>(spec: Spec | undefined, ctor: (field: Data.Field) => Field<T>): Schema<T> {
        return { type: 0 as any, ctor, undefinedField: (spec && spec.undefinedField) || DefaultUndefined, alias: spec && spec.alias };
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
    constructor(private _category: Data.Category, schema: Category.Schema) {
        const fieldKeys = Object.keys(schema).filter(k => k !== '@alias');
        const cache = Object.create(null);
        for (const k of fieldKeys) {
            const s = schema[k];
            Object.defineProperty(this, k, {
                get: function() {
                    if (cache[k]) return cache[k];
                    const field = _category.getField(s.alias || k) || s.undefinedField;
                    cache[k] = s.ctor(field);
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
    const cat = block.categories[schema['@alias'] || key] || Data.Category.Empty;
    return new _Category(cat, schema);
}