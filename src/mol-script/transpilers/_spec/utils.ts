/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panangiot_tourlov@hotmail.com>
 */

import { Transpiler } from '../transpiler';
import { KeywordDict, PropertyDict, OperatorList } from '../types';

/* FAULTY IMPORTS */
import compile from '../../reference-implementation/molql/compiler';

export function testKeywords(keywords: KeywordDict, transpiler: Transpiler) {
    for (const name in keywords) {
        it(name, () => {
            const k = keywords[name];
            if (k.map) {
                const expr = transpiler(name);
                compile(expr);
                expect(expr).toEqual(k.map());
            } else {
                const transpile = () => transpiler(name);
                expect(transpile).toThrow();
                expect(transpile).not.toThrowError(RangeError, 'Maximum call stack size exceeded');
            }
        });
    }
}

export function testProperties(properties: PropertyDict, transpiler: Transpiler) {
    for (const name in properties) {
        const p = properties[name];
        p['@examples'].forEach(example => {
            it(name, () => {
                if (!p.isUnsupported) {
                    const expr = transpiler(example);
                    compile(expr);
                } else {
                    const transpile = () => transpiler(example);
                    expect(transpile).toThrow();
                    expect(transpile).not.toThrowError(RangeError, 'Maximum call stack size exceeded');
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

export function testOperators(operators: OperatorList, transpiler: Transpiler) {
    operators.forEach(o => {
        o['@examples'].forEach(example => {
            it(o.name, () => {
                if (!o.isUnsupported) {
                    const expr = transpiler(example);
                    compile(expr);
                } else {
                    const transpile = () => transpiler(example);
                    expect(transpile).toThrow();
                    expect(transpile).not.toThrowError(RangeError, 'Maximum call stack size exceeded');
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