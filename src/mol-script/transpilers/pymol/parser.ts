/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

// https://pymol.org/dokuwiki/doku.php?id=selection:alpha
// https://pymolwiki.org/index.php/Selection_Algebra
// https://github.com/evonove/pymol/blob/master/pymol/layer3/Selector.cpp

import * as P from '../../../mol-util/monadic-parser';
import * as h from '../helper';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { properties } from './properties';
import { operators } from './operators';
import { keywords } from './keywords';
import { AtomGroupArgs } from '../types';
import { Transpiler } from '../transpiler';

const propertiesDict = h.getPropertyRules(properties);

const slash = P.MonadicParser.string('/');

function orNull(rule: P.MonadicParser<any>) {
    return rule.or(P.MonadicParser.of(null));
}

function atomSelectionQuery(x: any) {
    const tests: AtomGroupArgs = {};
    const props: { [k: string]: any[] } = {};

    for (const k in x) {
        const ps = properties[k];
        if (!ps) {
            throw new Error(`property '${k}' not supported, value '${x[k]}'`);
        }
        if (x[k] === null) continue;
        if (!props[ps.level]) props[ps.level] = [];
        props[ps.level].push(x[k]);
    }

    for (const p in props) {
        tests[p] = h.andExpr(props[p]);
    }

    return B.struct.generator.atomGroups(tests);
}

const lang = P.MonadicParser.createLanguage({
    Parens: function (r: any) {
        return P.MonadicParser.alt(
            r.Parens,
            r.Operator,
            r.Expression
        ).wrap(P.MonadicParser.string('('), P.MonadicParser.string(')'));
    },

    Expression: function (r: any) {
        return P.MonadicParser.alt(
            r.Keywords,
            r.AtomSelectionMacro.map(atomSelectionQuery),
            r.NamedAtomProperties,
            r.Pepseq,
            r.Rep,
            r.Object
        );
    },

    AtomSelectionMacro: function (r: any) {
        return P.MonadicParser.alt(
            slash.then(P.MonadicParser.alt(
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty).skip(slash),
                    orNull(propertiesDict.segi).skip(slash),
                    orNull(propertiesDict.chain).skip(slash),
                    orNull(propertiesDict.resi).skip(slash),
                    orNull(propertiesDict.name)
                ).map(x => { return { object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4] }; }),
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty).skip(slash),
                    orNull(propertiesDict.segi).skip(slash),
                    orNull(propertiesDict.chain).skip(slash),
                    orNull(propertiesDict.resi)
                ).map(x => { return { object: x[0], segi: x[1], chain: x[2], resi: x[3] }; }),
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty).skip(slash),
                    orNull(propertiesDict.segi).skip(slash),
                    orNull(propertiesDict.chain)
                ).map(x => { return { object: x[0], segi: x[1], chain: x[2] }; }),
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty).skip(slash),
                    orNull(propertiesDict.segi)
                ).map(x => { return { object: x[0], segi: x[1] }; }),
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty)
                ).map(x => { return { object: x[0] }; }),
            )),
            P.MonadicParser.alt(
                P.MonadicParser.seq(
                    orNull(r.ObjectProperty).skip(slash),
                    orNull(propertiesDict.segi).skip(slash),
                    orNull(propertiesDict.chain).skip(slash),
                    orNull(propertiesDict.resi).skip(slash),
                    orNull(propertiesDict.name)
                ).map(x => { return { object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.segi).skip(slash),
                    orNull(propertiesDict.chain).skip(slash),
                    orNull(propertiesDict.resi).skip(slash),
                    orNull(propertiesDict.name)
                ).map(x => { return { segi: x[0], chain: x[1], resi: x[2], name: x[3] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.chain).skip(slash),
                    orNull(propertiesDict.resi).skip(slash),
                    orNull(propertiesDict.name)
                ).map(x => { return { chain: x[0], resi: x[1], name: x[2] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.resi).skip(slash),
                    orNull(propertiesDict.name)
                ).map(x => { return { resi: x[0], name: x[1] }; }),
            )
        );
    },

    NamedAtomProperties: function () {
        return P.MonadicParser.alt(...h.getNamedPropertyRules(properties));
    },

    Keywords: () => P.MonadicParser.alt(...h.getKeywordRules(keywords)),

    ObjectProperty: () => {
        const w = h.getReservedWords(properties, keywords, operators)
            .sort(h.strLenSortFn).map(h.escapeRegExp).join('|');
        return P.MonadicParser.regexp(new RegExp(`(?!(${w}))[A-Z0-9_]+`, 'i'));
    },
    Object: (r: any) => {
        return r.ObjectProperty.notFollowedBy(slash)
            .map((x: any) => { throw new Error(`property 'object' not supported, value '${x}'`); });
    },

    // Selects peptide sequence matching upper-case one-letter
    // sequence SEQ (see also FindSeq).
    // PEPSEQ seq
    Pepseq: () => {
        return P.MonadicParser.regexp(/(PEPSEQ|ps\.)\s+([a-z]+)/i, 2)
            .map(h.makeError(`operator 'pepseq' not supported`));
    },

    // Selects atoms which show representation rep.
    // REP rep
    Rep: () => {
        return P.MonadicParser.regexp(/REP\s+(lines|spheres|mesh|ribbon|cartoon|sticks|dots|surface|labels|extent|nonbonded|nb_spheres|slice|extent|slice|dashes|angles|dihedrals|cgo|cell|callback|everything)/i, 1)
            .map(h.makeError(`operator 'rep' not supported`));
    },

    Operator: function (r: any) {
        return h.combineOperators(operators, P.MonadicParser.alt(r.Parens, r.Expression, r.Operator));
    },

    Query: function (r: any) {
        return P.MonadicParser.alt(
            r.Operator,
            r.Parens,
            r.Expression
        ).trim(P.MonadicParser.optWhitespace);
    }
});

export const transpiler: Transpiler = str => lang.Query.tryParse(str);
