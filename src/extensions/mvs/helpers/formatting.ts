/** Resolve Python-like f-string (simplified functionality: only variable names (no expressions), subset of format specification).
 * If any of the formatted values is `undefined`, return `undefined`. */
export function resolveFString(fstring: string, valueGetter: (name: string) => string | undefined): string | undefined {
    const out: string[] = [];
    /** 0 if currently not in a tag, or positive value i if current tag started at position i (excl. opening brace) */
    let opened = 0;
    for (let i = 0; i < fstring.length; i++) {
        const char = fstring[i];
        if (opened) {
            if (char === '{') {
                throw new Error('ValueError: Invalid format template ("{" within tag)');
            } else if (char === '}') {
                const chunk = resolveFormatTag(fstring.slice(opened, i), valueGetter);
                if (chunk === undefined) return undefined;
                out.push(chunk);
                opened = 0;
            } else {
                // DO NOTHING
            }
        } else {
            if (char === '{') {
                if (fstring[i + 1] === '{') {
                    out.push('{');
                    i++;
                } else {
                    opened = i + 1;
                }
            } else if (char === '}') {
                if (fstring[i + 1] === '}') {
                    out.push('}');
                    i++;
                } else {
                    throw new Error('ValueError: Invalid format template (unmatched "}")');
                }
            } else {
                out.push(char);
            }
        }
    }
    if (opened) {
        throw new Error('ValueError: Invalid format template (unmatched "{")');
    }
    return out.join('');
}

const FORMAT_SPEC_RE = /^(?:(?<fill>.?)(?<align>[<>=^]))?(?<sign>[-+ ]?)(?<z>z?)(?<alt>#?)(?<zeros>0?)(?<width>\d*)(?<grouping>[,_]?)\.?(?<precision>\d*)(?<type>[dfFeE%]?)$/;
// Python 3.13: [[fill]align][sign]["z"]["#"]["0"] [width][grouping]["." precision] [type]

const NUM_FORMAT_TYPES = 'bdeEfFgGnoxX%';
const STR_FORMAT_TYPES = 'cs';

/** Resolve a single f-string format tag, e.g. `age:.2f` -> `1.20` */
function resolveFormatTag(formatTag: string, valueGetter: (name: string) => string | undefined): string | undefined {
    const [name, formatSpec] = formatTag.split(':');
    const value = valueGetter(name);
    if (value === undefined) return undefined;
    if (formatSpec === undefined) return value;
    return formatValue(value, formatSpec);
}

/** Format a value a la f-string, e.g. `formatValue('1.2', '.2f')` -> `'1.20'` */
export function formatValue(value: string, formatSpec: string): string {
    const match = formatSpec.match(FORMAT_SPEC_RE);
    if (!match || !match.groups) throw new Error(`ValueError: Invalid formatting "${formatSpec}"`);

    const type = match.groups.type || 's';
    const isNumeric = NUM_FORMAT_TYPES.includes(type);
    const sign = match.groups.sign as '+' | '-' | ' ' | '';
    const z = match.groups.z as 'z' | '';
    const alt = match.groups.alt as '#' | '';
    const zeros = match.groups.zeros as '0' | '';
    const width = Number(match.groups.width);
    const grouping = match.groups.grouping as ',' | '_' | '';
    const precision = match.groups.precision ? Number(match.groups.precision) : 6;
    const fill = match.groups.fill || zeros || ' ';
    const align = (match.groups.align as '<' | '>' | '=' | '^' | '') || (isNumeric ? (zeros ? '=' : '<') : '>');

    if (z) throw new Error(`NotImplementedError: Formatting option "z" not supported`);
    if (alt) throw new Error(`NotImplementedError: Formatting option "#" not supported`);

    return alignString(
        formatWithoutAligning(value, sign, grouping, precision, type),
        width, align, fill,
    );
}

function formatWithoutAligning(value: string, sign: '+' | '-' | ' ' | '', grouping: ',' | '_' | '', precision: number, type: string): string {
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
            if (NUM_FORMAT_TYPES.includes(type) || STR_FORMAT_TYPES.includes(type)) {
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
