/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 **/


import * as P from '../../../mol-util/monadic-parser';
import * as h from '../helper';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { properties } from './properties';
import { macroproperties } from './macroproperties';
import { operators } from './operators';
import { keywords } from './keywords';
import { AtomGroupArgs } from '../types';
import { Transpiler } from '../transpiler';


function listMap(x: string) { return x.split(',').map(x => x.replace(/^["']|["']$/g, '')); }
function rangeMap(x: string) {
    const [min, max] = x.split('-').map(x => parseInt(x));
    return { min, max };
}
function listOrRangeMap(x: string) {
    if (x.includes('-') && x.includes(',')) {
        const pSplit = x.split(',').map(x => x.replace(/^["']|["']$/g, ''));
        const res: number[] = [];
        pSplit.forEach(x => {
	    if (x.includes('-')) {
                const [min, max] = x.split('-').map(x=>parseInt(x));
                for (let i = min; i <= max; i++) {
		    res.push(i);
                }
	    } else {
                res.push(parseInt(x));
	    }
        });
        return res;
    } else if (x.includes('-') && !x.includes(',')) {
	const res: number[] = [];
	const [min, max] = x.split('-').map(x=>parseInt(x));
        for (let i = min; i <= max; i++) {
	    res.push(i);
        }	
        return res;
    } else if (!x.includes('-') && x.includes(',')) {
        return listMap(x).map(x => parseInt(x));
    } else {
        return [parseInt(x)];
    }
}


const propertiesDict = h.getPropertyRules(macroproperties);

const dot = P.MonadicParser.string('.');
const colon = P.MonadicParser.string(':');
const star = P.MonadicParser.string('*');
const bra = P.MonadicParser.string('(');
const ket = P.MonadicParser.string(')');
const commu = P.MonadicParser.string('[');
const tator = P.MonadicParser.string(']');

function orNull(rule: P.MonadicParser<any>) {
    return rule.or(P.MonadicParser.of(null));
}

function atomSelectionQuery2(x: any) {
    const tests: AtomGroupArgs = {};
    const props: { [k: string]: any[] } = {};

    for (const k in x) {
        const ps = macroproperties[k];
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

function atomExpressionQuery(x: any[]) {
    const [resnorange, resno, inscode, chainname, atomname, altloc] = x[1];
    const tests: AtomGroupArgs = {};

    if (chainname) {
    // should be configurable, there is an option in Jmol to use auth or label
        tests['chain-test'] = B.core.rel.eq([B.ammp('auth_asym_id'), chainname]);
    }

    const resnoRangeProps:any = [];
    console.log(resnorange)
    if (resnorange){	
	resnorange.forEach((x:number) =>{
	resnoRangeProps.push(B.core.rel.eq([B.ammp('auth_seq_id'), x]));
	})
	console.log(resnoRangeProps);
    };
    if (resnoRangeProps.length) tests['residue-test'] = h.orExpr(resnoRangeProps);
    
    const resProps:any = [];
    if (resno){
	console.log(resno)
	resProps.push(B.core.rel.eq([B.ammp('auth_seq_id'), resno]));
    }
    if (inscode) resProps.push(B.core.rel.eq([B.ammp('pdbx_PDB_ins_code'), inscode]));
    if (resProps.length) tests['residue-test'] = h.andExpr(resProps);

    const atomProps = [];
    if (atomname) atomProps.push(B.core.rel.eq([B.ammp('auth_atom_id'), atomname]));
    if (altloc) atomProps.push(B.core.rel.eq([B.ammp('label_alt_id'), altloc]));
    if (atomProps.length) tests['atom-test'] = h.andExpr(atomProps);

    return B.struct.generator.atomGroups(tests);
}


const lang = P.MonadicParser.createLanguage({

    Integer: () => P.MonadicParser.regexp(/-?[0-9]+/).map(Number).desc('integer'),

    Parens: function (r: any) {
        return P.MonadicParser.alt(
            r.Parens,
            r.Operator,
            r.Expression
        ).wrap(P.MonadicParser.regexp(/\(\s+/), P.MonadicParser.regexp(/\s+\)/));
    },

    Expression: function (r: any) {
        return P.MonadicParser.alt(
	    // order matters
	    r.Keywords,
	    r.NamedAtomProperties,
	    r.AtomSelectionMacro.map(atomSelectionQuery2),
	    r.AtomExpression.map(atomExpressionQuery),
	    r.Element.map((x: string) => B.struct.generator.atomGroups({
                'atom-test': B.core.rel.eq([B.acp('elementSymbol'), B.struct.type.elementSymbol(x)])
            })),
	    r.Resname.map((x: string) => B.struct.generator.atomGroups({
                'residue-test': B.core.rel.eq([B.ammp('label_comp_id'), x])
	    })),
	    r.Object,
	    r.Object2,
        );
    },

    //    lys:a.ca  -> resn lys and chain A and name ca
    //    lys*a.ca  -> resn lys and chain A and name ca
    //
    //    :a.ca -> chain A and name ca
    //    *a.ca -> chain A and name ca
    //
    //    *.cg -> name ca
    //    :.cg -> name ca
    AtomSelectionMacro: function (r: any) {
        return P.MonadicParser.alt(
	    // :A.CA :.CA :A
            colon.then(P.MonadicParser.alt(
                P.MonadicParser.seq(
                    orNull(propertiesDict.chain).skip(dot),
                    orNull(propertiesDict.name)
                ).map(x => { return { chain: x[0], name: x[1] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.name).skip(dot)
                ).map(x => { return { name: x[0] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.chain)
                ).map(x => { return { chain: x[0] }; }),
            )),
	    // *A.CA *.CA
	    star.then(P.MonadicParser.alt(
                P.MonadicParser.seq(
                    orNull(propertiesDict.chain).skip(dot),
                    orNull(propertiesDict.name)
                ).map(x => { return { chain: x[0], name: x[1] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.name).skip(dot)
                ).map(x => { return { name: x[0] }; }),
                P.MonadicParser.seq(
                    orNull(propertiesDict.chain)
                ).map(x => { return { chain: x[0] }; }),
            )),
	    // 1-100,201
	    bra.then(P.MonadicParser.alt(
                P.MonadicParser.alt(
		    P.MonadicParser.seq(
                        propertiesDict.resi.skip(ket),
		    ).map(x => {
                        return { resi: x[0] }
                        ;
                    })
                )
	    )),
	    //  [lys]10:a.ca [lys]10:a [lys]10 [lys]10.ca
	    //  [lys]:a.ca [lys]:a [lys] [lys].ca
	    commu.then(P.MonadicParser.alt(
                P.MonadicParser.alt(
                    P.MonadicParser.alt(
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator),
			    orNull(propertiesDict.resi).skip(colon),
                            orNull(propertiesDict.chain).skip(dot),
                            orNull(propertiesDict.name)
                        ).map(x => { return { resn: x[0], resi: x[1], chain: x[2], name: x[3] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator),
			    orNull(propertiesDict.resi).skip(colon),
                            orNull(propertiesDict.chain)
                        ).map(x => { return { resn: x[0], resi: x[1], chain: x[2] }; }),
                        P.MonadicParser.seq(
			    orNull(propertiesDict.resn).skip(tator),
			    orNull(propertiesDict.resi).skip(colon).skip(dot),
                            orNull(propertiesDict.name)
                        ).map(x => { return { resn: x[0], resi: x[1], name: x[2] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator),
			    orNull(propertiesDict.resi).skip(dot),
                            orNull(propertiesDict.name)
                        ).map(x => { return { resn: x[0], resi: x[1], name: x[2] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator),
			    orNull(propertiesDict.resi)
                        ).map(x => { return { resn: x[0], resi: x[1] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator).skip(colon),
                            orNull(propertiesDict.chain).skip(dot),
                            orNull(propertiesDict.name)
                        ).map(x => { return { resn: x[0], chain: x[1], name: x[2] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator).skip(star),
                            orNull(propertiesDict.chain).skip(dot),
                            orNull(propertiesDict.name)
                        ).map(x => { return { resn: x[0], chain: x[1], name: x[2] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator).skip(colon),
                            orNull(propertiesDict.chain),
                        ).map(x => { return { resn: x[0], chain: x[1] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator).skip(star),
                            orNull(propertiesDict.chain),
                        ).map(x => { return { resn: x[0], chain: x[1] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator).skip(dot),
                            orNull(propertiesDict.name),
                        ).map(x => { return { resn: x[0], name: x[1] }; }),
                        P.MonadicParser.seq(
                            orNull(propertiesDict.resn).skip(tator),
                        ).map(x => { return { resn: x[0] }; })
		    )
                )
	    )
		      )
        );
    },

    ObjectProperty: () => {
        const w = h.getReservedWords(macroproperties, keywords, operators)
            .sort(h.strLenSortFn).map(h.escapeRegExp).join('|');
        return P.MonadicParser.regexp(new RegExp(`(?!(${w}))[A-Z0-9_]+`, 'i'));
    },

    Object: (r: any) => {
        return r.ObjectProperty
            .map((x: any) => { throw new Error(`property 'object' not supported, value '${x}'`); });
    },


    ObjectProperty2: () => {
        const w = h.getReservedWords(properties, keywords, operators)
            .sort(h.strLenSortFn).map(h.escapeRegExp).join('|');
        return P.MonadicParser.regexp(new RegExp(`(?!(${w}))[A-Z0-9_]+`, 'i'));
    },

    Object2: (r: any) => {
        return r.ObjectProperty2
            .map((x: any) => { throw new Error(`property 'object' not supported, value '${x}'`); });
    },

    NamedAtomProperties: function () {
        return P.MonadicParser.alt(...h.getNamedPropertyRules(properties));
    },

    Operator: function (r: any) {
        return h.combineOperators(operators, P.MonadicParser.alt(r.Parens, r.Expression, r.Operator));
    },

    AtomExpression: function (r: any) {
        return P.MonadicParser.seq(
            P.MonadicParser.lookahead(r.AtomPrefix),
            P.MonadicParser.seq(
		r.ResnoRange.or(P.MonadicParser.of(null)),
                r.Resno.or(P.MonadicParser.of(null)),
                r.Inscode.or(P.MonadicParser.of(null)),
                r.Chainname.or(P.MonadicParser.of(null)),
                r.Atomname.or(P.MonadicParser.of(null)),
                r.Altloc.or(P.MonadicParser.of(null)),
                r.Model.or(P.MonadicParser.of(null))
            )
        );
    },

    AtomPrefix: () => P.MonadicParser.regexp(/[0-9:^%/.]/).desc('atom-prefix'),
    Chainname: () => P.MonadicParser.regexp(/:([A-Za-z]{1,3})/, 1).desc('chainname'),
    Model: () => P.MonadicParser.regexp(/\/([0-9]+)/, 1).map(Number).desc('model'),
    Element: () => P.MonadicParser.regexp(/_([A-Za-z]{1,3})/, 1).desc('element'),
    Atomname: () => P.MonadicParser.regexp(/\.([a-zA-Z0-9]{1,4})/, 1).map(B.atomName).desc('atomname'),
    Resname: () => P.MonadicParser.regexp(/[A-Z0-9]{1,4}/).desc('resname'),
    Resno: (r: any) => r.Integer.desc('resno'),
    Altloc: () => P.MonadicParser.regexp(/%([a-zA-Z0-9])/, 1).desc('altloc'),
    Inscode: () => P.MonadicParser.regexp(/\^([a-zA-Z0-9])/, 1).desc('inscode'),


   ResnoRange: function (r:any) {
       return P.MonadicParser.regex(/[0-9,-]+/).map( listOrRangeMap ).desc('resnorange')
	//   // 123-200
    //   // -12--3
    },

    Keywords: () => P.MonadicParser.alt(...h.getKeywordRules(keywords)),

    Query: function (r: any) {
        return P.MonadicParser.alt(
            r.Operator,
            r.Parens,
            r.Expression
        ).trim(P.MonadicParser.optWhitespace);
    }

});

export const transpiler: Transpiler = str => lang.Query.tryParse(str);
