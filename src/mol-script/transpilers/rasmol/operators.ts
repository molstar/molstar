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
import { macroproperties } from './macroproperties';

const propNames = Object.keys(macroproperties).sort(h.strLenSortFn)
    .filter(name => !macroproperties[name].isUnsupported).join('|');


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
        rule: P.MonadicParser.alt(h.infixOp(/AND|&/i), P.MonadicParser.whitespace),
        // rule: h.infixOp(/AND|&/i),
        map: (op, selection, by) => B.struct.modifier.intersectBy({ 0: selection, by })
    },
    {
        '@desc': 'Selects atoms included in either s1 or s2.',
        '@examples': ['ASP or GLU'],
        name: 'or',
        type: h.binaryLeft,
        rule: h.infixOp(/OR|\||\|\|/i),
        map: (op, s1, s2) => B.struct.combinator.merge([s1, s2])
    },
    /*
    {
        '@desc': 'Selects atoms within a specified distance of a selection',
        '@examples': ['within (5.0, LYS:A.CA)'],
        name: 'aaa',
        type: h.prefixRemoveKet,
        rule: h.prefixOpNoWhiteSpace(/within\s*\(\s*([-+]?[0-9]*\.?[0-9]+)\s*,/, 1).map((x: any) => {
            console.log(x)
	    return parseFloat(x)}),
	map: (radius: number, selection: Expression) => {

            return B.struct.modifier.includeSurroundings({ 0: selection, radius });
        }
    },
  */
    {
        '@desc': 'Selects atoms in s1 that are within X Angstroms of any atom in s2.',
        '@examples': ['chain A WITHIN 3 OF chain B'],
        name: 'within',
        abbr: ['w2.'],
        //   type: h.binaryLeft,
        type: h.prefixRemoveKet,
        rule: h.prefixOpNoWhiteSpace(/within\s*\(\s*([-+]?[0-9]*\.?[0-9]+)\s*,/i, 1).map((x: any) => {
	    console.log(x);
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
    {
        '@desc': 'Selects atoms which have the same keyword as the atoms in a given selection',
        '@examples': ['same resid as name FE'],
        name: 'same',
        type: h.prefix,
        rule: h.prefixOp(new RegExp(`SAME\\s+(${propNames})\\s+AS`, 'i'), 1).map((x: any) => macroproperties[x].property),
        map: (property: Expression, source: Expression) => {
            return B.struct.filter.withSameAtomProperties({
                '0': B.struct.generator.all(),
                source,
                property
            });
        }
    },
    {
        '@desc':
            'Selects atoms in s1 whose identifiers name and resi match atoms in s2.',
        '@examples': ['chain A LIKE chain B'],
        name: 'like',
        type: h.binaryLeft,
        rule: h.infixOp(/LIKE|l\./i),
        map: (op: string, selection: Expression, source: Expression) => {
            return B.struct.filter.withSameAtomProperties({
                0: selection,
                source,
                property: B.core.type.compositeKey([
                    B.ammp('label_atom_id'),
                    B.ammp('label_seq_id'),
                ]),
            });
        },
    },
];

