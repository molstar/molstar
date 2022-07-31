/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 */

import * as u from './utils';
import { parse } from '../../transpile';
//import { transpiler } from '../pymol/parser';
import { keywords } from '../pymol/keywords';
import { properties } from '../pymol/properties';
import { operators } from '../pymol/operators';

/* FAULTY IMPORTS */
//import compile from '../../compiler';

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
            const expr = parse("pymol",str);
	    expect(expr).toThrow();
  //          compile(expr);
        });
    });
    general.unsupported.forEach(str => {
        it(str, () => {
            const transpileStr = () => parse("pymol",str);
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

describe('pymol keywords', () => u.testKeywords(keywords, "pymol"));
describe('pymol operators', () => u.testOperators(operators, "pymol"));
describe('pymol properties', () => u.testProperties(properties, "pymol"));
