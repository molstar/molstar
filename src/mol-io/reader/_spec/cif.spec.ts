/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from '../cif/data-model'
import TextField from '../cif/text/field'
import * as Schema from '../cif/schema'
import { Column } from 'mol-base/collections/database'

const columnData = `123abc`;

const intField = TextField({ data: columnData, indices: [0, 1, 1, 2, 2, 3], count: 3 }, 3);
const strField = TextField({ data: columnData, indices: [3, 4, 4, 5, 5, 6], count: 3 }, 3);

const testBlock = Data.Block(['atoms'], {
    atoms: Data.Category('atoms', 3, ['x', 'name'], {
        x: intField,
        name: strField
    })
}, 'test');

namespace TestSchema {
    export const atoms = { x: Column.Schema.int, name: Column.Schema.str }
    export const schema = { atoms }
}

describe('schema', () => {
    const db = Schema.toDatabase(TestSchema.schema, testBlock);
    it('property access', () => {
        const { x, name } = db.atoms;
        expect(x.value(0)).toBe(1);
        expect(name.value(1)).toBe('b');
    });

    it('toArray', () => {
        const ret = db.atoms.x.toArray({ array: Int32Array });
        expect(ret.length).toBe(3);
        expect(ret[0]).toBe(1);
        expect(ret[1]).toBe(2);
        expect(ret[2]).toBe(3);
    })
});