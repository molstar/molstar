/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const elm1 = [ 'H', 'C', 'O', 'N', 'S', 'P' ];
const elm2 = [ 'NA', 'CL', 'FE' ];

function charAtIsNumber(str: string, index: number) {
    const code = str.charCodeAt(index);
    return code >= 48 && code <= 57;
}

export function guessElement(str: string) {
    let at = str.trim().toUpperCase();

    if (charAtIsNumber(at, 0)) at = at.substr(1);
    // parse again to check for a second integer
    if (charAtIsNumber(at, 0)) at = at.substr(1);
    const n = at.length;

    if (n === 0) return '';
    if (n === 1) return at;
    if (n === 2) {
        if (elm2.indexOf(at) !== -1) return at;
        if (elm1.indexOf(at[0]) !== -1) return at[0];
    }
    if (n >= 3) {
        if (elm1.indexOf(at[0]) !== -1) return at[0];
    }
    return '';
}
