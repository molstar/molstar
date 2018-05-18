/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from 'mol-data/iterator'
import CIF, { CifCategory } from 'mol-io/reader/cif'
import * as Encoder from 'mol-io/writer/cif'
import * as fs from 'fs'
import classify from './field-classifier'

async function getCIF(path: string) {
    const str = fs.readFileSync(path, 'utf8');
    const parsed = await CIF.parseText(str).run();
    if (parsed.isError) {
        throw new Error(parsed.toString());
    }
    return parsed.result;
}

function createDefinition(cat: CifCategory): Encoder.CategoryDefinition {
    return {
        name: cat.name,
        fields: cat.fieldNames.map(n => classify(n, cat.getField(n)!))
    }
}

function getCategoryInstanceProvider(cat: CifCategory): Encoder.CategoryProvider {
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

    const encoder = Encoder.create({ binary: !asText, encoderName: 'mol* cif2bcif' });
    for (const b of cif.blocks) {
        encoder.startDataBlock(b.header);
        for (const c of b.categoryNames) {
            encoder.writeCategory(getCategoryInstanceProvider(b.categories[c]));
        }
    }
    return encoder.getData();
}
