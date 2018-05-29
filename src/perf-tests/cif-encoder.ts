/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif'

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

function getInstance(ctx: { name: string, fields: CifWriter.Field[], rowCount: number }): CifWriter.Category {
    return {
        data: void 0,
        name: ctx.name,
        fields: ctx.fields,
        rowCount: ctx.rowCount
    }
}

const enc = CifWriter.createEncoder();

const filter: CifWriter.Category.Filter = {
    includeCategory(cat) { return true; },
    includeField(cat, field) { return !(cat === 'cat2' && field === 'e2') }
}

enc.startDataBlock('test');
enc.setFilter(filter);
enc.writeCategory(getInstance, [{ rowCount: 5, name: 'cat1', fields: category1fields }]);
enc.writeCategory(getInstance, [{ rowCount: 1, name: 'cat2', fields: category2fields  }]);
console.log(enc.getData());
