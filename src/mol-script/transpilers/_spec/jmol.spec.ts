/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
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
