/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
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
    let ret = 0, neg = 1;
    if (str.charCodeAt(start) === 45 /* - */) { neg = -1; start++; }
    for (; start < end; start++) {
        let c = str.charCodeAt(start) - 48;
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
    let neg = 1.0, ret = 0.0, point = 0.0, div = 1.0;

    if (str.charCodeAt(start) === 45) {
        neg = -1.0;
        ++start;
    }

    while (start < end) {
        let c = str.charCodeAt(start) - 48;
        if (c >= 0 && c < 10) {
            ret = ret * 10 + c;
            ++start;
        } else if (c === -2) { // .
            ++start;
            while (start < end) {
                c = str.charCodeAt(start) - 48;
                if (c >= 0 && c < 10) {
                    point = 10.0 * point + c;
                    div = 10.0 * div;
                    ++start;
                } else if (c === 53 || c === 21) { // 'e'/'E'
                    return parseScientific(neg * (ret + point / div), str, start + 1, end);
                } else {
                    return neg * (ret + point / div);
                }
            }
            return neg * (ret + point / div);
        } else if (c === 53 || c === 21) { // 'e'/'E'
            return parseScientific(neg * ret, str, start + 1, end);
        }
        else break;
    }
    return neg * ret;
}
