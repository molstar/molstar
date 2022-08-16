/*
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 */


import * as P from '../../../mol-util/monadic-parser';
import * as h from '../helper';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { OperatorList } from '../types';
import { Expression } from '../../language/expression';

export const operators: OperatorList = [
    {
        '@desc': 'Selects atoms that are not included in s1.',
        '@examples': ['not [ARG]'],
        name: 'not',
        type: h.prefix,
        rule: P.MonadicParser.alt(P.MonadicParser.regex(/NOT/i).skip(P.MonadicParser.whitespace), P.MonadicParser.string('!').skip(P.MonadicParser.optWhitespace)),
        map: (op, selection) => h.invertExpr(selection),
    },
    {
        '@desc': 'Selects atoms included in both s1 and s2.',
        '@examples': ['[ASP] and :.CA'],
        name: 'and',
        type: h.binaryLeft,
        rule: P.MonadicParser.alt(h.infixOp(/AND|&/i)),
        map: (op, selection, by) => B.struct.modifier.intersectBy({ 0: selection, by })
    },
    {
        '@desc': 'Selects atoms included in either s1 or s2.',
        '@examples': ['[ASP] or [GLU]'],
        name: 'or',
        abbr: [','],
        type: h.binaryLeft,
        rule: P.MonadicParser.alt(h.infixOp(/OR|\||\|\||,/i)),
        map: (op, s1, s2) => B.struct.combinator.merge([s1, s2])
    },
    {
        '@desc': 'Selects atoms in s1 that are within X Angstroms of any atom in s2.',
        '@examples': ['within(5.0, [HEM])'],
        name: 'within',
        abbr: ['w2.'],
        type: h.prefixRemoveKet,
        rule: h.prefixOpNoWhiteSpace(/within\s*\(\s*([-+]?[0-9]*\.?[0-9]+)\s*,/i, 1).map((x: any) => {
	    return parseFloat(x);
        }),
        map: (radius: number, target: Expression) => {
	    return B.struct.filter.within({
                0: B.struct.generator.all(),
                target,
                'max-radius': radius,
	    });
        },
    },
];

