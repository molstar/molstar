/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panangiot_tourlov@hotmail.com>
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 */

import { Transpiler } from '../transpiler';
import { KeywordDict, PropertyDict, OperatorList } from '../types';

export function testKeywords(keywords: KeywordDict, transpiler: Transpiler) {
    for (const name in keywords) {
        it(name, () => {
            const k = keywords[name];
            if (k.map) {
                const expr = transpiler(name);
                expect(expr).toEqual(k.map());
            } else {
                const transpile = () => transpiler(name);
                expect(transpile).toThrow();
                expect(transpile).not.toThrowError(RangeError);
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
                    transpiler(example);
                } else {
                    const transpile = () => transpiler(example);
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

export function testOperators(operators: OperatorList, transpiler: Transpiler) {
    operators.forEach(o => {
        o['@examples'].forEach(example => {
            it(o.name, () => {
                if (!o.isUnsupported) {
                    transpiler(example);
                } else {
                    const transpile = () => transpiler(example);
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
