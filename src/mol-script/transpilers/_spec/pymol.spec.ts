/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 */

import * as u from './utils';
import { transpiler } from '../pymol/parser';
import { keywords } from '../pymol/keywords';
import { properties } from '../pymol/properties';
import { operators } from '../pymol/operators';

const general = {
    supported: [
        // macros
        '10/cb',
        'a/10-12/ca',
        'lig/b/6+8/c+o',

        // trimming
        '    name CA   ',
        'name CA   ',
        '    name CA',
    ],
    unsupported: [
        // macros
        'pept/enz/c/3/n',
        'pept/enz///n',

        '/pept/lig/',
        '/pept/lig/a',
        '/pept/lig/a/10',
        '/pept/lig/a/10/ca',
        '/pept//a/10',

        // object
        'foobar',
        'protein and bazbar',
    ]
};

describe('pymol general', () => {
    general.supported.forEach(str => {
        it(str, () => {
            transpiler(str);
            //          compile(expr);
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

// check against builder output
// 'not (resi 42 or chain A)'
// '!resi 42 or chain A'
// 'b >= 0.3',
// 'b != 0.3',
// 'b>0.3',
// 'b <0.3',
// 'b <= 0.3',
// 'b = 1',
// 'fc.=.1',

describe('pymol keywords', () => u.testKeywords(keywords, transpiler));
describe('pymol operators', () => u.testOperators(operators, transpiler));
describe('pymol properties', () => u.testProperties(properties, transpiler));
