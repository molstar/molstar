/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Iterator from 'mol-data/iterator'
import { Column, Table } from 'mol-data/db'
import { Tensor } from 'mol-math/linear-algebra'
import Encoder from '../encoder'

// TODO: support for "coordinate fields", make "coordinate precision" a parameter of the encoder
// TODO: automatically detect "precision" of floating point arrays.
// TODO: automatically detect "best encoding" for integer arrays. This could be used for "fixed-point" as well.
// TODO: add "repeat encoding"? [[1, 2], [1, 2], [1, 2]] --- Repeat ---> [[1, 2], 3]
// TODO: Add "higher level fields"? (i.e. generalization of repeat)
// TODO: align "data blocks" to 8 byte offsets for fast typed array windows? (prolly needs some testing if this is actually the case too)
// TODO: "parametric encoders"? Specify encoding as [{ param: 'value1', encoding1 }, { param: 'value2', encoding2 }]
//       then the encoder can specify { param: 'value1' } and the correct encoding will be used.
//       Use case: variable precision encoding for different fields.
//       Perhaps implement this as parameter spaces...

export const enum FieldType {
    Str, Int, Float
}

export interface FieldDefinitionBase<Key, Data> {
    name: string,
    valueKind?: (key: Key, data: Data) => Column.ValueKind,
    // TODO:
    // shouldInclude?: (data: Data) => boolean
}

export type FieldDefinition<Key = any, Data = any> =
    | FieldDefinitionBase<Key, Data> & { type: FieldType.Str, value(key: Key, data: Data): string }
    | FieldDefinitionBase<Key, Data> & { type: FieldType.Int, value(key: Key, data: Data): number }
    | FieldDefinitionBase<Key, Data> & { type: FieldType.Float, value(key: Key, data: Data): number }

export interface FieldFormat {
    // TODO: do we actually need this?
    // textDecimalPlaces: number,
    // stringEncoder: ArrayEncoder,
    // numericEncoder: ArrayEncoder,
    // typedArray?: E.TypedArrayCtor
}

export namespace FieldFormat {
    export const Default: FieldFormat = {
        // textDecimalPlaces: 3,
        // stringEncoder: ArrayEncoder.by(E.stringArray),
        // numericEncoder: ArrayEncoder.by(E.byteArray)
    };
}

export interface CategoryDefinition<Key = any, Data = any> {
    name: string,
    fields: FieldDefinition<Key, Data>[]
}

export interface CategoryInstance<Key = any, Data = any> {
    data: Data,
    definition: CategoryDefinition<Key, Data>,
    formats?: { [name: string]: Partial<FieldFormat> },
    rowCount: number,
    keys(): Iterator<Key>
}

export interface CategoryProvider {
    (ctx: any): CategoryInstance
}

export interface CIFEncoder<T = string | Uint8Array, Context = any> extends Encoder {
    startDataBlock(header: string): void,
    writeCategory(category: CategoryProvider, contexts?: Context[]): void,
    getData(): T
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
    const fieldDefinitions: FieldDefinition[] = []
    const type = FieldType.Float
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

export namespace FieldDefinitions {
    export function ofSchema(schema: Table.Schema) {
        const fields: FieldDefinition[] = [];
        for (const k of Object.keys(schema)) {
            const t = schema[k];
            if (t.valueType === 'int') {
                fields.push({ name: k, type: FieldType.Int, value: columnValue(k), valueKind: columnValueKind(k) });
            } else if (t.valueType === 'float') {
                fields.push({ name: k, type: FieldType.Float, value: columnValue(k), valueKind: columnValueKind(k) });
            } else if (t.valueType === 'str') {
                fields.push({ name: k, type: FieldType.Str, value: columnValue(k), valueKind: columnValueKind(k) });
            } else if (t.valueType === 'list') {
                fields.push({ name: k, type: FieldType.Str, value: columnListValue(k), valueKind: columnValueKind(k) })
            } else if (t.valueType === 'tensor') {
                fields.push(...getTensorDefinitions(k, t.space))
            } else {
                throw new Error(`Unknown valueType ${t.valueType}`);
            }
        }
        return fields;
    }
}

export namespace CategoryDefinition {
    export function ofTable<S extends Table.Schema>(name: string, table: Table<S>): CategoryDefinition<number> {
        return { name, fields: FieldDefinitions.ofSchema(table._schema) }
    }

    export function instanceProviderOfTable(name: string, table: Table<Table.Schema>): CategoryProvider {
        return function (ctx: any) {
            return {
                data: table,
                definition: ofTable(name, table),
                keys: () => Iterator.Range(0, table._rowCount - 1),
                rowCount: table._rowCount
            };
        }
    }
}
