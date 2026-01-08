/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */


/** Formatting template created from a Python-like f-string.
 * Supports simplified f-string functionality:
 * - only variable names (no expressions);
 * - supported types: `deEfF%cs`;
 * - not supported: types `bgGnoxX`, options `z` and `#`. */
export interface FormatTemplate {
    /** Source f-string for this template */
    fstring: string,
    /** Apply format template to values obtained by calling `valueGetter` on variable names. If any of the obtained values is `undefined`, return `undefined`. */
    format: (valueGetter: (name: string) => string | undefined) => string | undefined,
}

export function FormatTemplate(fstring: string): FormatTemplate {
    const { n, varNames, varFormats, literals } = parseFString(fstring);
    let _out: string[];
    return {
        fstring,
        format(valueGetter) {
            const out = _out ??= [];
            out.length = 0;
            for (let i = 0; i < n; i++) {
                out.push(literals[i]);
                const value = valueGetter(varNames[i]);
                if (value === undefined) return undefined;
                out.push(formatValue(value, varFormats[i]));
            }
            out.push(literals[n]);
            return out.join('');
        },
    };
}

interface ParsedFstring {
    /** Number of variables to be inserted in the template */
    n: number,
    /** Names of the n variables */
    varNames: string[],
    /** Formats for the n variables */
    varFormats: FormatSpec[],
    /** n+1 literal strings to be placed before, between, and after the variables */
    literals: string[],
}

/** Parse Python-like f-string */
function parseFString(fstring: string): ParsedFstring {
    const literals: string[] = [];
    const varNames: string[] = [];
    const varFormats: FormatSpec[] = [];
    /** Non-negative = where current literal started; negative = -where current tag started (excluding braces). */
    let start = 0;
    for (let i = 0; i < fstring.length; i++) {
        const char = fstring[i];
        if (start >= 0) {
            // In literal
            if (char === '{') {
                if (fstring[i + 1] === '{') {
                    i++; // skip the other {
                } else {
                    literals.push(fstring.slice(start, i).replace(/{{/g, '{').replace(/}}/g, '}'));
                    start = -(i + 1); // start tag
                }
            } else if (char === '}') {
                if (fstring[i + 1] === '}') {
                    i++; // skip the other }
                } else {
                    throw new Error('ValueError: Invalid format template (unmatched "}")');
                }
            } // else do nothing
        } else {
            // In tag
            if (char === '{') {
                throw new Error('ValueError: Invalid format template ("{" within tag)');
            } else if (char === '}') {
                const [varName, formatSpec] = parseFormatTag(fstring.slice(-start, i));
                varNames.push(varName);
                varFormats.push(formatSpec);
                start = i + 1; // start literal
            } // else do nothing
        }
    }
    if (start < 0) {
        throw new Error('ValueError: Invalid format template (unmatched "{")');
    }
    literals.push(fstring.slice(start, fstring.length).replace(/{{/g, '{').replace(/}}/g, '}'));

    return { n: varNames.length, varNames, varFormats, literals };
}

/** Parse a single f-string format tag, e.g. `age:.2f` */
function parseFormatTag(formatTag: string): [string, FormatSpec] {
    // cannot use .split because format spec may contain ':'
    let colonIndex = formatTag.indexOf(':');
    if (colonIndex < 0) colonIndex = formatTag.length;
    const varName = formatTag.slice(0, colonIndex);
    const formatSpec = formatTag.slice(colonIndex + 1, undefined);
    return [varName, parseFormatSpec(formatSpec)];
}


// Python 3.13: [[fill]align][sign]["z"]["#"]["0"] [width][grouping]["." precision] [type]
const FORMAT_SPEC_RE = /^(?:(?<fill>.?)(?<align>[<>=^]))?(?<sign>[-+ ]?)(?<z>z?)(?<alt>#?)(?<zeros>0?)(?<width>\d*)(?<grouping>[,_]?)\.?(?<precision>\d*)(?<type>[bdeEfFgGnoxX%cs]?)$/;
const FORMAT_TYPES = ['b', 'd', 'e', 'E', 'f', 'F', 'g', 'G', 'n', 'o', 'x', 'X', '%', 'c', 's'] as const;
type FormatType = typeof FORMAT_TYPES[number]

interface FormatSpec {
    /** Formatting type */
    type: FormatType,
    /** Controls adding explicit sign for non-negative numbers */
    sign: '+' | '-' | ' ' | '',
    /** Converts negative zeros to positive (NOT SUPPORTED) */
    z: 'z' | '',
    /** Use alternative form (NOT SUPPORTED) */
    alt: '#' | '',
    /** Sets fill char to '0' and align to '=' */
    zeros: '0' | '',
    /** Min required width of output string */
    width: number,
    /** Controls digit grouping in large numbers */
    grouping: ',' | '_' | '',
    /** Number of decimal digits for types eEfF% */
    precision: number,
    /** Fill character for padding */
    fillChar: string,
    /** Align direction for padding */
    align: '<' | '>' | '=' | '^',
}

/** Parse a single f-string format spec, e.g. `.2f` -> `{ type:'f', precision: 2, ... }` */
function parseFormatSpec(formatSpec: string): FormatSpec {
    const match = formatSpec.match(FORMAT_SPEC_RE);
    if (!match || !match.groups) throw new Error(`ValueError: Invalid formatting "${formatSpec}"`);

    const type = (match.groups.type || 's') as FormatSpec['type'];
    const isNumeric = type !== 's' && type !== 'c';
    const sign = match.groups.sign as '+' | '-' | ' ' | '';
    const z = match.groups.z as 'z' | '';
    const alt = match.groups.alt as '#' | '';
    const zeros = match.groups.zeros as '0' | '';
    const width = Number(match.groups.width);
    const grouping = match.groups.grouping as ',' | '_' | '';
    const precision = match.groups.precision ? Number(match.groups.precision) : 6;
    const fillChar = match.groups.fill || zeros || ' ';
    const align = (match.groups.align as '<' | '>' | '=' | '^' | '') || (isNumeric ? (zeros ? '=' : '<') : '>');

    return { type, sign, z, alt, zeros, width, grouping, precision, fillChar, align };
}

/** Format a value a la f-string, e.g. `formatValue('1.2', '.2f')` -> `'1.20'` */
export function formatValue(value: string, formatSpec: FormatSpec | string): string {
    if (typeof formatSpec === 'string') formatSpec = parseFormatSpec(formatSpec);
    if (formatSpec.z) throw new Error(`NotImplementedError: Formatting option "z" not supported`);
    if (formatSpec.alt) throw new Error(`NotImplementedError: Formatting option "#" not supported`);
    return alignString(
        formatWithoutAligning(value, formatSpec.sign, formatSpec.grouping, formatSpec.precision, formatSpec.type),
        formatSpec.width, formatSpec.align, formatSpec.fillChar,
    );
}

function formatWithoutAligning(value: string, sign: '+' | '-' | ' ' | '', grouping: ',' | '_' | '', precision: number, type: FormatType): string {
    let out = '';
    switch (type) {
        case 's':
            return value;
        case 'c':
            return String.fromCharCode(Number(value));
        case 'd':
            out = `${Math.floor(Number(value))}`;
            break;
        case 'f':
        case 'F':
            out = Number(value).toFixed(precision);
            break;
        case 'e':
        case 'E':
            out = Number(value).toExponential(precision).replace(/(e[+-])(\d)$/, '$10$2'); // pad exponent to 2 chars (like Python)
            if (type === 'E') out = out.replace('e', 'E');
            break;
        case '%':
            out = `${(100 * Number(value)).toFixed(precision)}%`;
            break;
        default:
            if (FORMAT_TYPES.includes(type)) {
                throw new Error(`NotImplementedError: Formatting code "${type}" not supported`);
            } else {
                throw new Error(`ValueError: Invalid format code "${type}"`);
            }
    }
    if (grouping === ',') {
        const match = out.match(/(-?\d+)(.*)/)!;
        if (match) {
            out = match[1].replace(/\B(?=(\d{3})+(?!\d))/g, ',') + match[2];
        }
    }
    if (sign === '+' && out[0] !== '-') out = '+' + out;
    if (sign === ' ' && out[0] !== '-') out = ' ' + out;
    return out;
}

function alignString(str: string, length: number, align: '<' | '>' | '=' | '^', fillChar: string): string {
    switch (align) {
        case '<': return str.padEnd(length, fillChar);
        case '>': return str.padStart(length, fillChar);
        case '=':
            if (str.startsWith('+') || str.startsWith('-') || str.startsWith(' ')) {
                return str[0] + str.slice(1, undefined).padStart(length - 1, fillChar);
            } else {
                return str.padStart(length, fillChar);
            }
        case '^':
            const nStart = Math.max(0, Math.floor((length - str.length) / 2));
            return str.padStart(str.length + nStart, fillChar).padEnd(length, fillChar);
    }
}
