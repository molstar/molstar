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
    const len = target.length;
    if (len > end - start) return false;
    for (let i = 0; i < len; i++) {
        if (str.charCodeAt(start + i) !== target.charCodeAt(i)) return false;
    }
    return true;
}

export function interpolate(str: string, params: { [k: string]: any }) {
    return str.replace(/\$\{(\w+)\}/g, (_, key) => String(params[key] ?? ''));
}

const idPathExpressionRe = /\$\{([^}]+)\}/g;

/** Replicate `substr(start, length)` using `substring`, clamping to available characters. */
function idSubstrLike(id: string, start: number, length: number): string {
    if (length <= 0) return '';
    let from = start;
    if (from < 0) from = Math.max(0, id.length + from);
    return id.substring(from, from + length);
}

function evaluateIdPathExpression(expr: string, id: string): string {
    const trimmed = expr.trim();
    if (trimmed === 'id') return id;
    if (trimmed === 'id.toLowerCase()') return id.toLowerCase();
    if (trimmed === 'id.toUpperCase()') return id.toUpperCase();
    const substrMatch = /^id\.(?:substr|substring)\((\d+),\s*(\d+)\)$/.exec(trimmed);
    if (substrMatch) return idSubstrLike(id, +substrMatch[1], +substrMatch[2]);
    const substringMatch = /^id\.substring\((\d+)\)$/.exec(trimmed);
    if (substringMatch) return id.substring(+substringMatch[1]);
    const sliceMatch = /^id\.slice\((\d+)(?:,\s*(\d+))?\)$/.exec(trimmed);
    if (sliceMatch) return id.slice(+sliceMatch[1], sliceMatch[2] !== undefined ? +sliceMatch[2] : undefined);
    throw new Error(`Unsupported id path expression: \${${expr}}`);
}

/** Validate that a path template only uses supported `${id…}` expressions. */
export function validateIdPathTemplate(template: string) {
    for (const match of template.matchAll(idPathExpressionRe)) {
        evaluateIdPathExpression(match[1], 'test-id');
    }
}

export function interpolateIdPath(template: string, id: string): string {
    return template.replace(idPathExpressionRe, (_, expr) => evaluateIdPathExpression(expr, id));
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

/**
 * Escape string for use in Javascript regex
 *
 * From https://stackoverflow.com/questions/3446170/escape-string-for-use-in-javascript-regex/6969486#6969486
 */
export function escapeRegExp(str: string) {
    return str.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); // $& means the whole matched string
}