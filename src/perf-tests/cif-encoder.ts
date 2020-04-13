/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from '../mol-io/writer/cif';
import * as fs from 'fs';

const category1fields: CifWriter.Field[] = [
    CifWriter.Field.str('f1', i => 'v' + i),
    CifWriter.Field.int('f2', i => i * i),
    CifWriter.Field.float('f3', i => Math.random()),
];

const category2fields: CifWriter.Field[] = [
    CifWriter.Field.str('e1', i => 'v\n' + i),
    CifWriter.Field.int('e2', i => i * i),
    CifWriter.Field.float('e3', i => Math.random()),
];

function getCat(name: string): CifWriter.Category {
    return {
        name,
        instance(ctx: { fields: CifWriter.Field[], rowCount: number }) {
            return { fields: ctx.fields, source: [{ rowCount: ctx.rowCount }] };
        }
    };
}

function testText() {
    const enc = CifWriter.createEncoder();

    const filter: CifWriter.Category.Filter = {
        includeCategory(cat) { return true; },
        includeField(cat, field) { return !(cat === 'cat2' && field === 'e2'); }
    };

    enc.startDataBlock('test');
    enc.setFilter(filter);
    enc.writeCategory(getCat('cat1'), [{ rowCount: 5, fields: category1fields }]);
    enc.writeCategory(getCat('cat2'), [{ rowCount: 1, fields: category2fields }]);
    console.log(enc.getData());
}

testText();


function testBinary() {
    const enc = CifWriter.createEncoder({ binary: true });

    const filter: CifWriter.Category.Filter = {
        includeCategory(cat) { return true; },
        includeField(cat, field) { return !(cat === 'cat2' && field === 'e2'); }
    };

    enc.startDataBlock('test');
    enc.setFilter(filter);
    enc.writeCategory(getCat('cat1'), [{ rowCount: 5, fields: category1fields }]);
    enc.writeCategory(getCat('cat2'), [{ rowCount: 1, fields: category2fields }]);
    enc.encode();
    const data = enc.getData() as Uint8Array;
    fs.writeFileSync('e:/test/mol-star/test.bcif', Buffer.from(data));
    console.log('written binary');
}

testBinary();