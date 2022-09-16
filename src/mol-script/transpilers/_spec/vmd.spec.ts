/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 */

import * as u from './utils';
import { transpiler } from '../vmd/parser';
import { keywords } from '../vmd/keywords';
import { properties } from '../vmd/properties';
import { operators } from '../vmd/operators';

const general = {
    supported: [
        // trimming
        '    name CA   ',
        'name CA   ',
        '    name CA',
    ],
    unsupported: [
        // variables
        'name $atomname',
        'protein and @myselection',

        // values outside of comparisons
        'foobar',
        '34',
        'name',
        'abs(-42)',
        'abs(21+21)',
        'sqr(3)',
        'sqr(x)',
        'sqr(x+33)',
        'protein or foobar',
        '34 and protein',
        'name or protein',
    ]
};

describe('vmd general', () => {
    general.supported.forEach(str => {
        it(str, () => {
            transpiler(str);
            // compile(expr);
        });
    });
    general.unsupported.forEach(str => {
        it(str, () => {
            const transpileStr = () => transpiler(str);
            expect(transpileStr).toThrow();
            expect(transpileStr).not.toThrowError(RangeError);
        });
    });
});

describe('vmd keywords', () => u.testKeywords(keywords, transpiler));
describe('vmd operators', () => u.testOperators(operators, transpiler));
describe('vmd properties', () => u.testProperties(properties, transpiler));
