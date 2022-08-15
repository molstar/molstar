/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 *
 * Adapted from MolQL src/transpile.ts
 */

import { Transpiler } from './transpilers/transpiler';
import { _transpiler } from './transpilers/all';
import { Expression } from './language/expression';
import { Script } from './script';
const transpiler: {[index: string]: Transpiler} = _transpiler;

export function parse(lang: Script.Language, str: string): Expression {
    try {
        const query = transpiler[lang](str);

        console.log(str);
        console.log(query);
        //	console.log(util.inspect(query, {depth: 20, color: true}))
        console.log('\n');

        return query;
    } catch (e) {
        console.log(str);
        console.log(e.message);
        console.log('\n');
        throw e;
    }
}
