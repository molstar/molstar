/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

 import { CIFCategory, CIFField, createCIFEncoder } from 'mol-io/writer/cif'

const category1fields: CIFField[] = [
    CIFField.str('f1', i => 'v' + i),
    CIFField.int('f2', i => i * i),
    CIFField.float('f3', i => Math.random()),
];

const category2fields: CIFField[] = [
    CIFField.str('e1', i => 'v\n' + i),
    CIFField.int('e2', i => i * i),
    CIFField.float('e3', i => Math.random()),
];

function getInstance(ctx: { name: string, fields: CIFField[], rowCount: number }): CIFCategory {
    return {
        data: void 0,
        name: ctx.name,
        fields: ctx.fields,
        rowCount: ctx.rowCount
    }
}

const w = createCIFEncoder();

w.startDataBlock('test');
w.writeCategory(getInstance, [{ rowCount: 5, name: 'cat1', fields: category1fields }]);
w.writeCategory(getInstance, [{ rowCount: 1, name: 'cat2', fields: category2fields  }]);
console.log(w.getData());
