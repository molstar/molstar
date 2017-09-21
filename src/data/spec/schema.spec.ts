/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from '../data'
import * as Schema from '../schema'

function Field(values: any[]): Data.Field {
    return {
        isDefined: true,
        str: row => '' + values[row],
        int: row => +values[row] || 0,
        float: row => +values[row] || 0,
        bin: row => null,

        presence: row => Data.ValuePresence.Present,
        areValuesEqual: (rowA, rowB) => values[rowA] === values[rowB],
        stringEquals: (row, value) => '' + values[row] === value,

        toStringArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = values[startRow + i]; }
            return ret;
        },
        toNumberArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = +values[startRow + i]; }
            return ret;
        }
    }
}

const testBlock = Data.Block({
    'atoms': Data.Category(2, {
        x: Field([1, 2]),
        name: Field(['C', 'O'])
    })
});

namespace TestSchema {
    export const atoms = { x: Schema.Field.float(), name: Schema.Field.str() }
    export const schema = { atoms }
}

describe('schema', () => {
    const data = Schema.apply(TestSchema.schema, testBlock);
    it('property access', () => {
        const { x, name } = data.atoms;
        expect(x.value(0)).toBe(1);
        expect(name.value(1)).toBe('O');
    });

    it('toArray', () => {
        const ret = data.atoms.x.toArray(0, 2, (s) => new Int32Array(s))!;
        expect(ret.length).toBe(2);
        expect(ret[0]).toBe(1);
        expect(ret[1]).toBe(2);
    })
});