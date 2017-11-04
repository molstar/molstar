/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from 'mol-base/collections/iterator'
import CIF, { Category } from 'mol-io/reader/cif'
import TextCIFEncoder from 'mol-io/writer/cif/encoder/text'
import BinaryCIFEncoder from 'mol-io/writer/cif/encoder/binary'
import * as Encoder from 'mol-io/writer/cif/encoder'
import * as fs from 'fs'
import classify from './field-classifier'

async function getCIF(path: string) {
    const str = fs.readFileSync(path, 'utf8');
    const parsed = await CIF.parseText(str)();
    if (parsed.isError) {
        throw new Error(parsed.toString());
    }
    return parsed.result;
}

function createDefinition(cat: Category): Encoder.CategoryDefinition {
    return {
        name: cat.name,
        fields: cat.fieldNames.map(n => classify(n, cat.getField(n)!))
    }
}

function getCategoryInstanceProvider(cat: Category): Encoder.CategoryProvider {
    return function (ctx: any) {
        return {
            data: cat,
            definition: createDefinition(cat),
            keys: () => Iterator.Range(0, cat.rowCount - 1),
            rowCount: cat.rowCount
        };
    }
}

export default async function convert(path: string, asText = false) {
    const cif = await getCIF(path);

    const encoder = asText ? new TextCIFEncoder() : new BinaryCIFEncoder('mol* cif2bcif');
    for (const b of cif.blocks) {
        encoder.startDataBlock(b.header);
        for (const _c of Object.keys(b.categories)) {
            encoder.writeCategory(getCategoryInstanceProvider(b.categories[_c]));
        }
    }
    return encoder.getData();
}

