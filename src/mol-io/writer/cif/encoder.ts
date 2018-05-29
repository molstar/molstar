/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Iterator from 'mol-data/iterator'
import { Column, Table, Database, DatabaseCollection } from 'mol-data/db'
import { Tensor } from 'mol-math/linear-algebra'
import Encoder from '../encoder'
import { ArrayEncoder, ArrayEncoding } from '../../common/binary-cif';

// TODO: support for "coordinate fields", make "coordinate precision" a parameter of the encoder
// TODO: automatically detect "precision" of floating point arrays.
// TODO: automatically detect "best encoding" for integer arrays. This could be used for "fixed-point" as well.
// TODO: add "repeat encoding"? [[1, 2], [1, 2], [1, 2]] --- Repeat ---> [[1, 2], 3]
// TODO: Add "higher level fields"? (i.e. generalization of repeat)
// TODO: align "data blocks" to 8 byte offsets for fast typed array windows? (prolly needs some testing if this is actually the case too)

export interface CIFField<Key = any, Data = any> {
    name: string,
    type: CIFField.Type,
    valueKind?: (key: Key, data: Data) => Column.ValueKind,
    defaultFormat?: CIFField.Format,
    value(key: Key, data: Data): string | number
}

export namespace CIFField {
    export const enum Type { Str, Int, Float }

    export interface Format {
        digitCount?: number,
        encoder?: ArrayEncoder,
        typedArray?: ArrayEncoding.TypedArrayCtor
    }

    export function getDigitCount(field: CIFField) {
        if (field.defaultFormat && typeof field.defaultFormat.digitCount !== 'undefined') return field.defaultFormat.digitCount;
        return 6;
    }

    export function str<K, D = any>(name: string, value: (k: K, d: D) => string, params?: { valueKind?: (k: K, d: D) => Column.ValueKind, encoder?: ArrayEncoder }): CIFField<K, D> {
        return { name, type: Type.Str, value, valueKind: params && params.valueKind, defaultFormat: params && params.encoder ? { encoder: params.encoder } : void 0 };
    }

    export function int<K, D = any>(name: string, value: (k: K, d: D) => number, params?: { valueKind?: (k: K, d: D) => Column.ValueKind, encoder?: ArrayEncoder, typedArray?: ArrayEncoding.TypedArrayCtor }): CIFField<K, D> {
        return {
            name,
            type: Type.Int,
            value,
            valueKind: params && params.valueKind,
            defaultFormat: params ? { encoder: params.encoder, typedArray: params.typedArray } : void 0 };
    }

    export function float<K, D = any>(name: string, value: (k: K, d: D) => number, params?: { valueKind?: (k: K, d: D) => Column.ValueKind, encoder?: ArrayEncoder, typedArray?: ArrayEncoding.TypedArrayCtor, digitCount?: number }): CIFField<K, D> {
        return {
            name,
            type: Type.Float,
            value,
            valueKind: params && params.valueKind,
            defaultFormat: params ? { encoder: params.encoder, typedArray: params.typedArray, digitCount: typeof params.digitCount !== 'undefined' ? params.digitCount : void 0 } : void 0
        };
    }
}

export interface CIFCategory<Key = any, Data = any> {
    name: string,
    fields: CIFField<Key, Data>[],
    data?: Data,
    rowCount: number,
    keys?: () => Iterator<Key>
}

export namespace CIFCategory {
    export interface Provider<Ctx = any> {
        (ctx: Ctx): CIFCategory
    }

    export function ofTable(name: string, table: Table<Table.Schema>): CIFCategory<number, Table<Table.Schema>> {
        return { name, fields: cifFieldsFromTableSchema(table._schema), data: table, rowCount: table._rowCount };
    }
}

export interface CIFEncoder<T = string | Uint8Array> extends Encoder {
    // setFormatter(): void,
    startDataBlock(header: string): void,
    writeCategory<Ctx>(category: CIFCategory.Provider<Ctx>, contexts?: Ctx[]): void,
    getData(): T
}

export namespace CIFEncoder {
    export function writeDatabase(encoder: CIFEncoder, name: string, database: Database<Database.Schema>) {
        encoder.startDataBlock(name);
        for (const table of database._tableNames) {
            encoder.writeCategory(() => CIFCategory.ofTable(table, database[table]));
        }
    }

    export function writeDatabaseCollection(encoder: CIFEncoder, collection: DatabaseCollection<Database.Schema>) {
        for (const name of Object.keys(collection)) {
            writeDatabase(encoder, name, collection[name])
        }
    }
}


function columnValue(k: string) {
    return (i: number, d: any) => d[k].value(i);
}

function columnListValue(k: string) {
    return (i: number, d: any) => d[k].value(i).join(d[k].schema.separator);
}

function columnTensorValue(k: string, ...coords: number[]) {
    return (i: number, d: any) => d[k].schema.space.get(d[k].value(i), ...coords);
}

function columnValueKind(k: string) {
    return (i: number, d: any) => d[k].valueKind(i);
}

function getTensorDefinitions(field: string, space: Tensor.Space) {
    const fieldDefinitions: CIFField[] = []
    const type = CIFField.Type.Float
    const valueKind = columnValueKind(field)
    if (space.rank === 1) {
        const rows = space.dimensions[0]
        for (let i = 0; i < rows; i++) {
            const name = `${field}[${i + 1}]`
            fieldDefinitions.push({ name, type, value: columnTensorValue(field, i), valueKind })
        }
    } else if (space.rank === 2) {
        const rows = space.dimensions[0], cols = space.dimensions[1]
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const name = `${field}[${i + 1}][${j + 1}]`
                fieldDefinitions.push({ name, type, value: columnTensorValue(field, i, j), valueKind })
            }
        }
    } else if (space.rank === 3) {
        const d0 = space.dimensions[0], d1 = space.dimensions[1], d2 = space.dimensions[2]
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) {
                for (let k = 0; k < d2; k++) {
                    const name = `${field}[${i + 1}][${j + 1}][${k + 1}]`
                    fieldDefinitions.push({ name, type, value: columnTensorValue(field, i, j, k), valueKind })
                }
            }
        }
    } else {
        throw new Error('Tensors with rank > 3 or rank 0 are currently not supported.')
    }
    return fieldDefinitions
}

function cifFieldsFromTableSchema(schema: Table.Schema) {
    const fields: CIFField[] = [];
    for (const k of Object.keys(schema)) {
        const t = schema[k];
        if (t.valueType === 'int') {
            fields.push({ name: k, type: CIFField.Type.Int, value: columnValue(k), valueKind: columnValueKind(k) });
        } else if (t.valueType === 'float') {
            fields.push({ name: k, type: CIFField.Type.Float, value: columnValue(k), valueKind: columnValueKind(k) });
        } else if (t.valueType === 'str') {
            fields.push({ name: k, type: CIFField.Type.Str, value: columnValue(k), valueKind: columnValueKind(k) });
        } else if (t.valueType === 'list') {
            fields.push({ name: k, type: CIFField.Type.Str, value: columnListValue(k), valueKind: columnValueKind(k) })
        } else if (t.valueType === 'tensor') {
            fields.push(...getTensorDefinitions(k, t.space))
        } else {
            throw new Error(`Unknown valueType ${t.valueType}`);
        }
    }
    return fields;
}