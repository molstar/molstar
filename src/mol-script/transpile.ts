/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma
 */
/**
 * Adapted from MolQL src/transpile.ts
 */

import { Transpiler } from './transpilers/transpiler';
import { _transpiler } from './transpilers/all';
const transpiler: {[index: string]: Transpiler} = _transpiler;

export function parse(lang: string, str: string) {
    try {
        const query = transpiler[lang](str);

        console.log(str);
        //	console.log(util.inspect(query, {depth: 20, color: true}))
        console.log('\n');

        return query;
    } catch (e) {
        console.log(str);
        console.log(e.message);
        console.log('\n');
    }
}
