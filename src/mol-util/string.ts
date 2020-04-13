/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

const reLine = /^/mg;
export function indentString(str: string, count: number, indent: string) {
    return count === 0 ? str : str.replace(reLine, indent.repeat(count));
}

/** Add space between camelCase text. */
export function splitCamelCase(str: string, separator = ' ') {
    return str.replace(/([a-z\xE0-\xFF])([A-Z\xC0\xDF])/g, `$1${separator}$2`);
}

/** Split camelCase text and capitalize. */
export function camelCaseToWords(str: string) {
    return capitalize(splitCamelCase(str));
}

export const lowerCase = (str: string) => str.toLowerCase();
export const upperCase = (str: string) => str.toUpperCase();

/** Return upper case if string, otherwise return empty string */
export function upperCaseAny(value: any): string {
    if (!value) return '';
    return typeof value === 'string' ? value.toUpperCase() : `${value}`.toUpperCase();
}

/** Uppercase the first character of each word. */
export function capitalize(str: string) {
    return str.toLowerCase().replace(/^\w|\s\w/g, upperCase);
}

export function splitSnakeCase(str: string) {
    return str.replace(/_/g, ' ');
}

export function snakeCaseToWords(str: string) {
    return capitalize(splitSnakeCase(str));
}

export function splitKebabCase(str: string) {
    return str.replace(/-/g, ' ');
}

export function kebabCaseToWords(str: string) {
    return capitalize(splitKebabCase(str));
}

export function stringToWords(str: string) {
    return capitalize(splitCamelCase(splitSnakeCase(splitKebabCase(str))));
}

export function substringStartsWith(str: string, start: number, end: number, target: string) {
    let len = target.length;
    if (len > end - start) return false;
    for (let i = 0; i < len; i++) {
        if (str.charCodeAt(start + i) !== target.charCodeAt(i)) return false;
    }
    return true;
}

export function interpolate(str: string, params: { [k: string]: any }) {
    const names = Object.keys(params);
    const values = Object.values(params);
    return new Function(...names, `return \`${str}\`;`)(...values);
}

export function trimChar(str: string, char: string) {
    let start = 0;
    let end = str.length;
    while (start < end && str[start] === char) ++start;
    while (end > start && str[end - 1] === char) --end;
    return (start > 0 || end < str.length) ? str.substring(start, end) : str;
}

export function trimCharStart(str: string, char: string) {
    let start = 0;
    const end = str.length;
    while (start < end && str[start] === char) ++start;
    return (start > 0) ? str.substring(start, end) : str;
}

export function trimCharEnd(str: string, char: string) {
    let end = str.length;
    while (end > 0 && str[end - 1] === char) --end;
    return (end < str.length) ? str.substring(0, end) : str;
}

/** Simple function to strip tags from a string */
export function stripTags(str: string) {
    return str.replace(/<\/?[^>]+>/g, '');
}