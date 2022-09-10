/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from MolQL project
 */

import * as P from '../../../mol-util/monadic-parser';
import * as h from '../helper';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { OperatorList } from '../types';

export const operators: OperatorList = [
    {
        '@desc': 'Selects atoms that are not included in s1.',
        '@examples': ['not ARG'],
        name: 'not',
        type: h.prefix,
        rule: P.MonadicParser.alt(P.MonadicParser.regex(/NOT/i).skip(P.MonadicParser.whitespace), P.MonadicParser.string('!').skip(P.MonadicParser.optWhitespace)),
        map: (op, selection) => h.invertExpr(selection),
    },
    {
        '@desc': 'Selects atoms included in both s1 and s2.',
        '@examples': ['ASP and .CA'],
        name: 'and',
        type: h.binaryLeft,
        rule: h.infixOp(/AND|&/i),
        map: (op, selection, by) => B.struct.modifier.intersectBy({ 0: selection, by })
    },
    {
        '@desc': 'Selects atoms included in either s1 or s2.',
        '@examples': ['ASP or GLU'],
        name: 'or',
        type: h.binaryLeft,
        rule: h.infixOp(/OR|\||,/i),
        map: (op, s1, s2) => B.struct.combinator.merge([s1, s2])
    }
];

