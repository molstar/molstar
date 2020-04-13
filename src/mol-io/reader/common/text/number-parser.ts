/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * based in part on https://github.com/dsehnal/CIFTools.js
 */

/**
 * Efficient integer and float parsers.
 *
 * For the purposes of parsing numbers from the mmCIF data representations,
 * up to 4 times faster than JS parseInt/parseFloat.
 */

export function parseIntSkipLeadingWhitespace(str: string, start: number, end: number) {
    while (start < end && str.charCodeAt(start) === 32) start++;
    return parseInt(str, start, end);
}

export function parseInt(str: string, start: number, end: number) {
    let _start = start, ret = 0, neg = 1;
    if (str.charCodeAt(_start) === 45 /* - */) { neg = -1; _start++; }
    for (; _start < end; _start++) {
        const c = str.charCodeAt(_start) - 48;
        if (c > 9 || c < 0) return (neg * ret) | 0;
        else ret = (10 * ret + c) | 0;
    }
    return neg * ret;
}

function parseScientific(main: number, str: string, start: number, end: number) {
    // handle + in '1e+1' separately.
    if (str.charCodeAt(start) === 43 /* + */) start++;
    return main * Math.pow(10.0, parseInt(str, start, end));
}

export function parseFloatSkipLeadingWhitespace(str: string, start: number, end: number) {
    while (start < end && str.charCodeAt(start) === 32) start++;
    return parseFloat(str, start, end);
}

export function parseFloat(str: string, start: number, end: number) {
    let _start = start, neg = 1.0, ret = 0.0, point = 0.0, div = 1.0;

    if (str.charCodeAt(_start) === 45) {
        neg = -1.0;
        ++_start;
    }

    while (_start < end) {
        let c = str.charCodeAt(_start) - 48;
        if (c >= 0 && c < 10) {
            ret = ret * 10 + c;
            ++_start;
        } else if (c === -2) { // .
            ++_start;
            while (_start < end) {
                c = str.charCodeAt(_start) - 48;
                if (c >= 0 && c < 10) {
                    point = 10.0 * point + c;
                    div = 10.0 * div;
                    ++_start;
                } else if (c === 53 || c === 21) { // 'e'/'E'
                    return parseScientific(neg * (ret + point / div), str, _start + 1, end);
                } else {
                    return neg * (ret + point / div);
                }
            }
            return neg * (ret + point / div);
        } else if (c === 53 || c === 21) { // 'e'/'E'
            return parseScientific(neg * ret, str, _start + 1, end);
        } else {
            break;
        }
    }
    return neg * ret;
}

export const enum NumberType {
    Int,
    Float,
    Scientific,
    NaN
}

function isInt(str: string, start: number, end: number) {
    if (str.charCodeAt(start) === 45 /* - */) { start++; }
    for (; start < end; start++) {
        const c = str.charCodeAt(start) - 48;
        if (c > 9 || c < 0) return false;
    }
    return true;
}

// TODO: check for "scientific integers?"
function getNumberTypeScientific(str: string, start: number, end: number) {
    // handle + in '1e+1' separately.
    if (str.charCodeAt(start) === 43 /* + */) start++;
    return isInt(str, start, end) ? NumberType.Scientific : NumberType.NaN;
}

/** The whole range must match, otherwise returns NaN */
export function getNumberType(str: string): NumberType {
    let start = 0, end = str.length;

    if (str.charCodeAt(start) === 45) { // -
        ++start;
    }

    // string is . or -.
    if (str.charCodeAt(start) === 46 && end - start === 1) {
        return NumberType.NaN;
    }

    while (start < end) {
        let c = str.charCodeAt(start) - 48;
        if (c >= 0 && c < 10) {
            ++start;
        } else if (c === -2) { // .
            ++start;
            let hasDigit = false;
            while (start < end) {
                c = str.charCodeAt(start) - 48;
                if (c >= 0 && c < 10) {
                    hasDigit = true;
                    ++start;
                } else if (c === 53 || c === 21) { // 'e'/'E'
                    return getNumberTypeScientific(str, start + 1, end);
                } else {
                    return NumberType.NaN;
                }
            }
            return hasDigit ? NumberType.Float : NumberType.Int;
        } else if (c === 53 || c === 21) { // 'e'/'E'
            if (start === 0 || start === 1 && str.charCodeAt(0) === 45) {
                return NumberType.NaN; // string starts with e/E or -e/-E
            }
            return getNumberTypeScientific(str, start + 1, end);
        } else {
            break;
        }
    }
    return start === end ? NumberType.Int : NumberType.NaN;
}
