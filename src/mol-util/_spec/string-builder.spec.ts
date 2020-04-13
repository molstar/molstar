/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SB from '../string-builder';

describe('string-builder', () => {

    function check(name: string, bb: (sb: SB) => void, expected: string) {
        const sb = SB.create();
        bb(sb);
        it(name, () => expect(SB.getString(sb)).toEqual(expected));
    }

    check('write', sb => SB.write(sb, '123'), '123');
    check('whitespace', sb => SB.whitespace(sb, 3), '   ');
    check('writePadLeft', sb => SB.writePadLeft(sb, '1', 3), '  1');
    check('writePadRight', sb => SB.writePadRight(sb, '1', 3), '1  ');
    check('writeIntegerPadLeft', sb => SB.writeIntegerPadLeft(sb, -125, 5), ' -125');
    check('writeIntegerPadRight', sb => SB.writeIntegerPadRight(sb, -125, 5), '-125 ');
    check('writeFloat', sb => SB.writeFloat(sb, 1.123, 100), '1.12');
    check('writeFloatPadLeft', sb => SB.writeFloatPadLeft(sb, 1.123, 100, 6), '  1.12');
    check('writeFloatPadRight', sb => SB.writeFloatPadRight(sb, -1.123, 100, 6), '-1.12 ');

    it('chunks', () => {
        const sb = SB.create(2);
        SB.write(sb, '1');
        SB.write(sb, '2');
        SB.write(sb, '3');

        expect(SB.getChunks(sb)).toEqual(['12', '3']);
        expect(SB.getString(sb)).toEqual('123');
    });
});