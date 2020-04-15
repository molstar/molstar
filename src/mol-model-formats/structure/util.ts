/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TokenBuilder, Tokens } from '../../mol-io/reader/common/text/tokenizer';

export function guessElementSymbolTokens(tokens: Tokens, str: string, start: number, end: number) {
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
            ((c === 78 || c === 110) && (c2 === 65 || c2 ===  97)) || // NA na Na nA
            ((c === 67 || c ===  99) && (c2 === 76 || c2 === 108)) || // CL
            ((c === 70 || c === 102) && (c2 === 69 || c2 === 101)) || // FE
            ((c === 83 || c === 115) && (c2 === 73 || c2 === 105))    // SI
        ) return TokenBuilder.add(tokens, s, s + 2);
    }

    if (
        c === 67 || c ===  99 || // C c
        c === 72 || c === 104 || // H h
        c === 78 || c === 110 || // N n
        c === 79 || c === 111 || // O o
        c === 80 || c === 112 || // P p
        c === 83 || c === 115    // S s
    ) return TokenBuilder.add(tokens, s, s + 1);

    TokenBuilder.add(tokens, s, s); // no reasonable guess, add empty token
}

const reTrimSpacesAndNumbers = /^[\s\d]+|[\s\d]+$/g;
export function guessElementSymbolString(str: string) {
    // trim spaces and numbers, convert to upper case
    str = str.replace(reTrimSpacesAndNumbers, '').toUpperCase();
    const l = str.length;

    if (l === 0) return str; // empty
    if (l === 1) return str; // one char

    if (l === 2) { // two chars
        if (str === 'NA' || str === 'CL' || str === 'FE' || str === 'SI') return str;
    }

    const c = str[0];
    if (c === 'C' || c === 'H' || c === 'N' || c === 'O' || c === 'P' || c === 'S') return c;

    return ''; // no reasonable guess, return empty string
}