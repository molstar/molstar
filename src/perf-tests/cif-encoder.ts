/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from 'mol-base/collections/iterator'
import * as Enc from 'mol-io/writer/cif/encoder'
import CW from 'mol-io/writer/cif/encoder/text'

const category1: Enc.CategoryDefinition<number> = {
    name: 'test',
    fields: [{
        name: 'f1',
        type: Enc.FieldType.Str,
        value: i => 'v' + i
    }, {
        name: 'f2',
        type: Enc.FieldType.Int,
        value: i => i * i
    }, {
        name: 'f3',
        type: Enc.FieldType.Float,
        value: i => Math.random()
    }]
}

const category2: Enc.CategoryDefinition<number> = {
    name: 'test2',
    fields: [{
        name: 'e1',
        type: Enc.FieldType.Str,
        value: i => 'v\n' + i
    }, {
        name: 'e2',
        type: Enc.FieldType.Int,
        value: i => i * i
    }, {
        name: 'e3',
        type: Enc.FieldType.Float,
        value: i => Math.random()
    }]
}

function getInstance(ctx: { cat: Enc.CategoryDefinition<number>, rowCount: number }): Enc.CategoryInstance {
    return {
        data: void 0,
        definition: ctx.cat,
        keys: () => Iterator.Range(0, ctx.rowCount - 1),
        rowCount: ctx.rowCount
    }
}

const w = new CW();

w.startDataBlock('test');
w.writeCategory(getInstance, [{ rowCount: 5, cat: category1 }]);
w.writeCategory(getInstance, [{ rowCount: 1, cat: category2 }]);
console.log(w.getData());
