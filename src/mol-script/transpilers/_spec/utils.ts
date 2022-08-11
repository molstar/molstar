/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panangiot_tourlov@hotmail.com>
 */

import { parse } from '../../transpile';
import { KeywordDict, PropertyDict, OperatorList } from '../types';

/* FAULTY IMPORTS */
// import compile from '../../compiler';

export function testKeywords(keywords: KeywordDict, language: string) {
    for (const name in keywords) {
        it(name, () => {
            const k = keywords[name];
            if (k.map) {
                const expr = parse(language, name);
                //              compile(expr);
                expect(expr).toEqual(k.map());
            } else {
                const transpile = () => parse(language, name);
                expect(transpile).toThrow();
                expect(transpile).not.toThrowError(RangeError);
            }
        });
    }
}

export function testProperties(properties: PropertyDict, language: string) {
    for (const name in properties) {
        const p = properties[name];
        p['@examples'].forEach(example => {
            it(name, () => {
                if (!p.isUnsupported) {
                    const expr = parse(language, example);
		    expect(expr).toBe(p);
                    //                    compile(expr);
                } else {
                    const transpile = () => parse(language, example);
                    expect(transpile).toThrow();
                    expect(transpile).not.toThrowError(RangeError);
                }
            });
        });
        it(name, () => {
            if (!p['@examples'].length) {
                throw Error(`'${name}' property has no example(s)`);
            }
        });
    }
}

export function testOperators(operators: OperatorList, language: string) {
    operators.forEach(o => {
        o['@examples'].forEach(example => {
            it(o.name, () => {
                if (!o.isUnsupported) {
                    const expr = parse(language, example);
		    expect(expr).toBe(o);
                } else {
                    const transpile = () => parse(language, example);
                    expect(transpile).toThrow();
                    expect(transpile).not.toThrowError(RangeError);
                }
            });
        });
        it(o.name, () => {
            if (!o['@examples'].length) {
                throw Error(`'${o.name}' operator has no example(s)`);
            }
        });
    });
}
