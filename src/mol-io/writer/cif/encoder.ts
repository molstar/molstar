/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from 'mol-base/collections/iterator'
import { Column } from 'mol-base/collections/database'
import Encoder from '../encoder'
//import { ArrayEncoder, ArrayEncoding as E } from '../../common/binary-cif'

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
    // TODO
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

export interface CIFEncoder<T = string | Uint8Array, Context = any> extends Encoder<T> {
    startDataBlock(header: string): void,
    writeCategory(category: CategoryProvider, contexts?: Context[]): void,
    getData(): T
}