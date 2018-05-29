/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import CIF, { CifCategory } from 'mol-io/reader/cif'
import { CIFCategory, createCIFEncoder } from 'mol-io/writer/cif'
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

function getCategoryInstanceProvider(cat: CifCategory): CIFCategory.Provider {
    return function (ctx: any) {
        return {
            data: cat,
            name: cat.name,
            fields: cat.fieldNames.map(n => classify(n, cat.getField(n)!)),
            rowCount: cat.rowCount
        };
    }
}

export default async function convert(path: string, asText = false) {
    const cif = await getCIF(path);

    const encoder = createCIFEncoder({ binary: !asText, encoderName: 'mol* cif2bcif' });
    for (const b of cif.blocks) {
        encoder.startDataBlock(b.header);
        for (const c of b.categoryNames) {
            encoder.writeCategory(getCategoryInstanceProvider(b.categories[c]));
        }
    }
    return encoder.getData();
}
