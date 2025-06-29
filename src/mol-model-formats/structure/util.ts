/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StringLike } from '../../mol-io/common/string-like';
import { TokenBuilder, Tokens } from '../../mol-io/reader/common/text/tokenizer';


export function guessElementSymbolTokens(tokens: Tokens, str: StringLike, start: number, end: number) {
    let s = start, e = end - 1;

    // trim spaces and numbers
    let c = str.charCodeAt(s);
    while ((c === 32 || (c >= 48 && c <= 57)) && s <= e) c = str.charCodeAt(++s);
    c = str.charCodeAt(e);
    while ((c === 32 || (c >= 48 && c <= 57)) && e >= s) c = str.charCodeAt(--e);

    ++e;

    if (s === e) return TokenBuilder.add(tokens, s, e); // empty
    if (s + 1 === e) return TokenBuilder.add(tokens, s, e); // one char

    c = str.charCodeAt(s);

    if (s + 2 === e) { // two chars
        const c2 = str.charCodeAt(s + 1);
        if (
            ((c === 78 || c === 110) && (c2 === 65 || c2 === 97)) || // NA na Na nA
            ((c === 67 || c === 99) && (c2 === 76 || c2 === 108)) || // CL
            ((c === 70 || c === 102) && (c2 === 69 || c2 === 101)) || // FE
            ((c === 83 || c === 115) && (c2 === 73 || c2 === 105)) || // SI
            ((c === 66 || c === 98) && (c2 === 82 || c2 === 114)) || // BR
            ((c === 65 || c === 97) && (c2 === 83 || c2 === 115)) // AS
        ) return TokenBuilder.add(tokens, s, s + 2);
    }

    if (
        c === 67 || c === 99 || // C c
        c === 72 || c === 104 || // H h
        c === 78 || c === 110 || // N n
        c === 79 || c === 111 || // O o
        c === 80 || c === 112 || // P p
        c === 83 || c === 115 // S s
    ) return TokenBuilder.add(tokens, s, s + 1);

    TokenBuilder.add(tokens, s, s); // no reasonable guess, add empty token
}

const TwoCharElementNames = new Set(['NA', 'CL', 'FE', 'SI', 'BR', 'AS', 'LI']);
const OneCharElementNames = new Set(['C', 'H', 'N', 'O', 'P', 'S', 'F', 'B']);

const reTrimSpacesAndNumbers = /^[\s\d]+|[\s\d]+$/g;
export function guessElementSymbolString(atomId: string, compId: string) {
    // trim spaces and numbers, convert to upper case
    atomId = atomId.replace(reTrimSpacesAndNumbers, '').toUpperCase();
    const l = atomId.length;

    if (l === 0) return atomId; // empty
    if (l === 1) return atomId; // one char
    if (TwoCharElementNames.has(atomId)) return atomId; // two chars

    // check for Charmm ion names where component and atom id are the same
    if (l === 3 && compId === atomId) {
        if (atomId === 'SOD') return 'NA';
        if (atomId === 'POT') return 'K';
        if (atomId === 'CES') return 'CS';
        if (atomId === 'CAL') return 'CA';
        if (atomId === 'CLA') return 'CL';
    }

    if (OneCharElementNames.has(atomId[0])) return atomId[0];

    return ''; // no reasonable guess, return empty string
}

