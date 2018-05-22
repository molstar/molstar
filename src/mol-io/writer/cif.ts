/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, DatabaseCollection } from 'mol-data/db'

import TextCIFEncoder from './cif/encoder/text'
import BinaryCIFEncoder from './cif/encoder/binary'
import { CategoryDefinition } from './cif/encoder'

export * from './cif/encoder'

export type EncoderInstance = BinaryCIFEncoder<{}> | TextCIFEncoder<{}>

export function create(params?: { binary?: boolean, encoderName?: string }): EncoderInstance {
    const { binary = false, encoderName = 'mol*' } = params || {};
    return binary ? new BinaryCIFEncoder(encoderName) : new TextCIFEncoder();
}

export function writeDatabase(encoder: EncoderInstance, name: string, database: Database<Database.Schema>) {
    encoder.startDataBlock(name);
    for (const table of database._tableNames) {
        encoder.writeCategory(
            CategoryDefinition.instanceProviderOfTable(table, database[table])
        );
    }
}

export function writeDatabaseCollection(encoder: EncoderInstance, collection: DatabaseCollection<Database.Schema>) {
    for (const name of Object.keys(collection)) {
        writeDatabase(encoder, name, collection[name])
    }
}