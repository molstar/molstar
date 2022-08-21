/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 */

import * as u from './utils';
import { transpiler } from '../rasmol/parser';
import { keywords } from '../rasmol/keywords';
import { properties } from '../rasmol/properties';
import { operators } from '../rasmol/operators';

const general = {
    supported: [
        // value comparison
        'resno > 10',
        // atom expression
        '[LEU]100:A.CA',
        '[LEU]100:A',
        '[LEU]100:.CA',
        '[LEU]:A.CA',
        '[LEU].CA',
        // residue numbering
        '(1-10,15,21-30)',
        // within with parentheses
        '( within(5,[HEM]) ) and backbone',
        // trimming
        '[ALA] and [VAL]  ',
        ' [ALA] and [VAL] ',
        '  [ALA] and [VAL]',
        // within with whitespaces
        'within (   5 ,  [HEM] ) ',
        // un-braketed residue name
        'LEU and ILE',
        // un-parenthesized residue index range
        '100-120,220',
        // un-parenthesized residue index
        '20',
    ],
    unsupported: [
        // within in the head or the middle of sentence
        'within (   5 ,  [HEM] ) and backbone',
    ]
};

describe('rasmol general', () => {
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

describe('rasmol keywords', () => u.testKeywords(keywords, transpiler));
describe('rasmol operators', () => u.testOperators(operators, transpiler));
describe('rasmol properties', () => u.testProperties(properties, transpiler));
