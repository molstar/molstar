/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from MolQL project
 */

import * as u from './utils';
import { transpiler } from '../jmol/parser';
import { keywords } from '../jmol/keywords';
import { properties } from '../jmol/properties';
import { operators } from '../jmol/operators';

const general = {
    supported: [
        // atom expressions
        '123',
        '-42',
        '_C',
        '.CA',
        'ALA',
        '%A',
        '^B',
        ':C',
        '/2',
        '10^A:F.CA%C/0',
        '10^A:F.CA%C',
        '10^A:F.CA',
        '10^A:F',
        '10^A',
        '10:F.CA',
        '10/0',
        '32 or 42',
        '.CA/0 OR 42:A',
        '!23',
        'not ASP',
        '(ASP or .CA)',
        'ASP and .CA',
        '123.CA',
        '(1 or 2) and .CA',
        '(1 or 2) and (.CA or .N)',
        '.CA and (2 or 3)',
        '.CA and (2 or 3) and ^A',
        '!32 or :A and .CA',

        // trimming
        '    atomName = CA   ',
        'atomName = CA   ',
        '    atomName = CA',

        // value comparison
        'resno > 10',
        // atom expression
        '[LEU]100:A.CA',
        '[LEU]100:A',
        '[LEU]100.CA',
        '[LEU]:A.CA',
        '[LEU].CA',
        // comma as OR
        '100, 42, ALA',
        // residue numbering
        '(1-10,15,21-30)',
        // within
        'within(5,[HEM])',
        // within with parentheses
        '(within(5,[HEM])) and backbone',
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
        // within in the head or the middle of sentence
        'within (   5 ,  [HEM] ) and backbone',

        // atom expressions with ranges
        '19-32:A',
        '-2-32:B',
        '-10--2:C',
        '[1FO]19-32:A',
    ],
    unsupported: [
        // values outside of comparisons
        'foobar',
        'protein or foobar',
    ]
};

describe('jmol general', () => {
    general.supported.forEach(str => {
        it(str, () => {
            transpiler(str);
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

describe('jmol keywords', () => u.testKeywords(keywords, transpiler));
describe('jmol properties', () => u.testProperties(properties, transpiler));
describe('jmol operators', () => u.testOperators(operators, transpiler));
