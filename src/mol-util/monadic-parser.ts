/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Adapted from Parsimmon (https://github.com/jneen/parsimmon)
 * Copyright (c) 2011-present J. Adkisson (http://jneen.net).
 */

export class MonadicParser<A> {
    constructor(public _: MonadicParser.Action<A>) { }

    parse(input: string): MonadicParser.ParseResult<A> {
        const result = this.skip(MonadicParser.eof)._(input, 0);
        if (result.status) {
            return { success: true, value: result.value };
        }
        return { success: false, index: makeLineColumnIndex(input, result.furthest), expected: result.expected };
    };

    tryParse(str: string) {
        const result = this.parse(str);
        if (result.success) {
            return result.value;
        } else {
            const msg = formatError(str, result);
            const err = new Error(msg);
            throw err;
        }
    }

    or<B>(alternative: MonadicParser<B>): MonadicParser<A | B> {
        return MonadicParser.alt(this, alternative);
    }

    trim<B>(parser: MonadicParser<B> | string): MonadicParser<A> {
        return this.wrap(parser, parser);
    }

    wrap<L, R>(leftParser: MonadicParser<L> | string, rightParser: MonadicParser<R> | string): MonadicParser<A> {
        return seqPick(1,
            typeof leftParser === 'string' ? MonadicParser.string(leftParser) : leftParser,
            this,
            typeof rightParser === 'string' ? MonadicParser.string(rightParser) : rightParser);
    }

    thru<B>(wrapper: (p: MonadicParser<A>) => MonadicParser<B>) {
        return wrapper(this);
    }

    then<B>(next: MonadicParser<B>): MonadicParser<B> {
        return seqPick(1, this, next);
    }

    many() {
        return new MonadicParser((input, i) => {
            const accum: A[] = [];
            let result: MonadicParser.Result<A> | undefined = void 0;

            while (true) {
                result = mergeReplies(this._(input, i), result);
                if (result.status) {
                    if (i === result.index) {
                        throw new Error('infinite loop detected in .many() parser --- calling .many() on a parser which can accept zero characters is usually the cause');
                    }
                    i = result.index;
                    accum.push(result.value);
                } else {
                    return mergeReplies(makeSuccess(i, accum), result);
                }
            }
        });
    };

    times(min: number, _max?: number): MonadicParser<A[]> {
        const max = typeof _max === 'undefined' ? min : _max;
        return new MonadicParser((input, i) => {
            const accum: A[] = [];
            let result: MonadicParser.Result<A> | undefined = void 0;
            let prevResult: MonadicParser.Result<A> | undefined = void 0;
            let times: number;
            for (times = 0; times < min; times++) {
                result = this._(input, i);
                prevResult = mergeReplies(result, prevResult);
                if (result.status) {
                    i = result.index;
                    accum.push(result.value);
                } else {
                    return prevResult as any;
                }
            }
            for (; times < max; times += 1) {
                result = this._(input, i);
                prevResult = mergeReplies(result, prevResult);
                if (result.status) {
                    i = result.index;
                    accum.push(result.value);
                } else {
                    break;
                }
            }
            return mergeReplies(makeSuccess(i, accum), prevResult);
        });
    };

    result<B>(res: B) {
        return this.map(() => res);
    };

    atMost(n: number) {
        return this.times(0, n);
    };

    atLeast(n: number) {
        return MonadicParser.seq(this.times(n), this.many()).map(r => [...r[0], ...r[1]]);
    };

    map<B>(f: (a: A) => B): MonadicParser<B> {
        return new MonadicParser((input, i) => {
            const result = this._(input, i);
            if (!result.status) {
                return result;
            }
            return mergeReplies(makeSuccess(result.index, f(result.value)), result);
        });
    }

    skip<B>(next: MonadicParser<B>): MonadicParser<A> {
        return seqPick(0, this, next);
    }

    mark(): MonadicParser<MonadicParser.Mark<A>> {
        return MonadicParser.seq(MonadicParser.index, this, MonadicParser.index).map(r => ({ start: r[0], value: r[1], end: r[2] }));
    }

    node(name: string): MonadicParser<MonadicParser.Node<A>> {
        return MonadicParser.seq(MonadicParser.index, this, MonadicParser.index).map(r => ({ name, start: r[0], value: r[1], end: r[2] }));
    };

    sepBy<B>(separator: MonadicParser<B>): MonadicParser<A[]> {
        return MonadicParser.sepBy(this, separator);
    }

    sepBy1<B>(separator: MonadicParser<B>): MonadicParser<A[]> {
        return MonadicParser.sepBy1(this, separator);
    }

    lookahead<B>(x: MonadicParser<B>) {
        return this.skip(MonadicParser.lookahead(x));
    };

    notFollowedBy<B>(x: MonadicParser<B>) {
        return this.skip(MonadicParser.notFollowedBy(x));
    };

    desc(expected: string) {
        return new MonadicParser((input, i) => {
            const reply = this._(input, i);
            if (!reply.status) {
                reply.expected = [expected];
            }
            return reply;
        });
    };

    fallback<B>(result: B) {
        return this.or(MonadicParser.succeed(result));
    };

    ap<B>(other: MonadicParser<(x: A) => B>): MonadicParser<B> {
        return MonadicParser.seq(other, this).map(([f, x]) => f(x));
    };

    chain<B>(f: (a: A) => MonadicParser<B>): MonadicParser<B> {
        return new MonadicParser<B>((input, i) => {
            const result = this._(input, i);
            if (!result.status) {
                return result as any;
            }
            const nextParser = f(result.value);
            return mergeReplies(nextParser._(input, result.index), result);
        });
    };
}

export namespace MonadicParser {
    export type Action<T> = (input: string, i: number) => MonadicParser.Result<T>

    export type ParseResult<T> = ParseSuccess<T> | ParseFailure;

    export interface Index {
        /** zero-based character offset */
        offset: number;
        /** one-based line offset */
        line: number;
        /** one-based column offset */
        column: number;
    }

    export interface ParseSuccess<T> {
        success: true,
        value: T
    }

    export interface ParseFailure {
        success: false,
        index: Index,
        expected: string[],
    }

    export interface Mark<T> {
        start: Index;
        end: Index;
        value: T;
    }

    export interface Node<T> extends Mark<T> {
        name: string
    }

    export interface Success<T> {
        status: true,
        value: T,
        index: number
    }

    export interface Failure {
        status: false,
        furthest: number,
        expected: string[]
    }

    export type Result<T> = Success<T> | Failure

    // export function createLanguage(parsers: any) {
    //     const language: any = {};
    //     for (const key of Object.keys(parsers)) {
    //         (function (key) {
    //             language[key] = lazy(() => parsers[key](language));
    //         })(key);
    //     }
    //     return language;
    // }

    export function seq<A>(a: MonadicParser<A>): MonadicParser<[A]>
    export function seq<A, B>(a: MonadicParser<A>, b: MonadicParser<B>): MonadicParser<[A, B]>
    export function seq<A, B, C>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>): MonadicParser<[A, B, C]>
    export function seq<A, B, C, D>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>, d: MonadicParser<D>): MonadicParser<[A, B, C, D]>
    export function seq<A, B, C, D, E>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>, d: MonadicParser<D>, e: MonadicParser<E>): MonadicParser<[A, B, C, D, E]>
    export function seq<T>(...parsers: MonadicParser<T>[]): MonadicParser<T[]>
    export function seq(...parsers: MonadicParser<any>[]): MonadicParser<any[]> {
        const numParsers = parsers.length;
        return new MonadicParser<any[]>((input, index) => {
            let result: MonadicParser.Result<any> | undefined;
            let accum = new Array(numParsers);
            let i = index;
            for (let j = 0; j < numParsers; j++) {
                result = mergeReplies(parsers[j]._(input, i), result);
                if (!result.status) {
                    return result;
                }
                accum[j] = result.value;
                i = result.index;
            }
            return mergeReplies(makeSuccess(i, accum), result);
        });
    }

    export function alt<A>(a: MonadicParser<A>): MonadicParser<A>
    export function alt<A, B>(a: MonadicParser<A>, b: MonadicParser<B>): MonadicParser<A | B>
    export function alt<A, B, C>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>): MonadicParser<A | B | C>
    export function alt<A, B, C, D>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>, d: MonadicParser<D>): MonadicParser<A | B | C | D>
    export function alt<A, B, C, D, E>(a: MonadicParser<A>, b: MonadicParser<B>, c: MonadicParser<C>, d: MonadicParser<D>, e: MonadicParser<E>): MonadicParser<A | B | C | D | E>
    export function alt<T>(...parsers: MonadicParser<T>[]): MonadicParser<T[]>
    export function alt(...parsers: MonadicParser<any>[]): MonadicParser<any> {
        const numParsers = parsers.length;
        if (numParsers === 0) {
            return fail('zero alternates');
        }
        return new MonadicParser((input, i) => {
            let result: MonadicParser.Result<any> | undefined;
            for (let j = 0; j < parsers.length; j++) {
                result = mergeReplies(parsers[j]._(input, i), result);
                if (result.status) {
                    return result;
                }
            }
            return result!;
        });
    }

    export function sepBy<A, B>(parser: MonadicParser<A>, separator: MonadicParser<B>): MonadicParser<A[]> {
        return sepBy1(parser, separator).or(succeed([]));
    }

    export function sepBy1<A, B>(parser: MonadicParser<A>, separator: MonadicParser<B>) {
        const pairs = separator.then(parser).many();
        return seq(parser, pairs).map(r => [r[0], ...r[1]]);
    }

    export function string(str: string) {
        const expected = `'${str}'`;
        if (str.length === 1) {
            const code = str.charCodeAt(0);
            return new MonadicParser((input, i) => input.charCodeAt(i) === code ? makeSuccess(i + 1, str) : makeFailure(i, expected));
        }

        return new MonadicParser((input, i) => {
            const j = i + str.length;
            if (input.slice(i, j) === str) return makeSuccess(j, str);
            else return makeFailure(i, expected);
        });
    }

    function flags(re: RegExp) {
        const s = '' + re;
        return s.slice(s.lastIndexOf('/') + 1);
    }

    function anchoredRegexp(re: RegExp) {
        return RegExp('^(?:' + re.source + ')', flags(re));
    }

    export function regexp(re: RegExp, group = 0) {
        const anchored = anchoredRegexp(re);
        const expected = '' + re;
        return new MonadicParser(function (input, i) {
            const match = anchored.exec(input.slice(i));
            if (match) {
                if (0 <= group && group <= match.length) {
                    const fullMatch = match[0];
                    const groupMatch = match[group];
                    return makeSuccess(i + fullMatch.length, groupMatch);
                }
                const message = `invalid match group (0 to ${match.length}) in ${expected}`;
                return makeFailure(i, message);
            }
            return makeFailure(i, expected);
        });
    }

    export function succeed<A>(value: A) {
        return new MonadicParser((input, i) => makeSuccess(i, value));
    }

    export function fail(expected: string): MonadicParser<any> {
        return new MonadicParser((input, i) => makeFailure(i, expected));
    }

    export function lookahead<A>(x: MonadicParser<A> | string | RegExp): MonadicParser<null> {
        if (isParser(x)) {
            return new MonadicParser((input, i) => {
                const result = x._(input, i);
                if (result.status) {
                    result.index = i;
                    result.value = null as any;
                }
                return result as any;
            });
        } else if (typeof x === 'string') {
            return lookahead(string(x));
        } else if (x instanceof RegExp) {
            return lookahead(regexp(x));
        }
        throw new Error('not a string, regexp, or parser: ' + x);
    }

    export function notFollowedBy<A>(parser: MonadicParser<A>): MonadicParser<null> {
        return new MonadicParser((input, i) => {
            const result = parser._(input, i);
            return result.status
                ? makeFailure(i, 'not "' + input.slice(i, result.index) + '"')
                : makeSuccess(i, null);
        });
    }

    export function test(predicate: (char: string) => boolean): MonadicParser<string> {
        return new MonadicParser((input, i) => {
            const char = input.charAt(i);
            if (i < input.length && predicate(char)) {
                return makeSuccess(i + 1, char);
            } else {
                return makeFailure(i, 'a character ' + predicate);
            }
        });
    }

    export function oneOf(str: string) {
        return test(ch => str.indexOf(ch) >= 0);
    }

    export function noneOf(str: string) {
        return test(ch => str.indexOf(ch) < 0);
    }

    export function range(begin: string, end: string) {
        return test(ch => begin <= ch && ch <= end).desc(begin + '-' + end);
    }

    export function takeWhile(predicate: (ch: string) => boolean) {
        return new MonadicParser((input, i) => {
            let j = i;
            while (j < input.length && predicate(input.charAt(j))) {
                j++;
            }
            return makeSuccess(j, input.slice(i, j));
        });
    }

    export function lazy<T>(f: () => MonadicParser<T>) {
        const parser = new MonadicParser((input, i) => {
            const a = f()._;
            parser._ = a;
            return a(input, i);
        });
        return parser;
    }

    export function empty() {
        return fail('empty');
    }

    export const index = new MonadicParser(function (input, i) {
        return makeSuccess(i, makeLineColumnIndex(input, i));
    });

    export const anyChar = new MonadicParser<string>((input, i) => {
        if (i >= input.length) {
            return makeFailure(i, 'any character');
        }
        return makeSuccess(i + 1, input.charAt(i));
    });

    export const all = new MonadicParser(function (input, i) {
        return makeSuccess(input.length, input.slice(i));
    });

    export const eof = new MonadicParser(function (input, i) {
        if (i < input.length) {
            return makeFailure(i, 'EOF');
        }
        return makeSuccess(i, null);
    });

    export const digit = regexp(/[0-9]/).desc('a digit');
    export const digits = regexp(/[0-9]*/).desc('optional digits');
    export const letter = regexp(/[a-z]/i).desc('a letter');
    export const letters = regexp(/[a-z]*/i).desc('optional letters');
    export const optWhitespace = regexp(/\s*/).desc('optional whitespace');
    export const whitespace = regexp(/\s+/).desc('whitespace');
    export const cr = string('\r');
    export const lf = string('\n');
    export const crlf = string('\r\n');
    export const newline = alt(crlf, lf, cr).desc('newline');
    export const end = alt(newline, eof);
}

function seqPick(idx: number, ...parsers: MonadicParser<any>[]): MonadicParser<any> {
    const numParsers = parsers.length;
    return new MonadicParser<any[]>((input, index) => {
        let result: MonadicParser.Result<any> | undefined;
        let picked: any;
        let i = index;
        for (let j = 0; j < numParsers; j++) {
            result = mergeReplies(parsers[j]._(input, i), result);
            if (!result.status) {
                return result;
            }
            if (idx === j) picked = result.value;
            i = result.index;
        }
        return mergeReplies(makeSuccess(i, picked), result);
    });
}

function makeSuccess<T>(index: number, value: T): MonadicParser.Success<T> {
    return { status: true, index, value };
}

function makeFailure(index: number, expected: string): MonadicParser.Failure {
    return { status: false, furthest: index, expected: [expected] };
}

function mergeReplies<A, B>(result: MonadicParser.Result<A>, last?: MonadicParser.Result<B>): MonadicParser.Result<A> {
    if (!last || result.status || last.status || result.furthest > last.furthest) {
        return result;
    }
    const expected = result.furthest === last.furthest
        ? unsafeUnion(result.expected, last.expected)
        : last.expected;
    return { status: result.status, furthest: last.furthest, expected };
}

function makeLineColumnIndex(input: string, i: number): MonadicParser.Index {
    const lines = input.slice(0, i).split('\n');
    // Note that unlike the character offset, the line and column offsets are
    // 1-based.
    const lineWeAreUpTo = lines.length;
    const columnWeAreUpTo = lines[lines.length - 1].length + 1;
    return { offset: i, line: lineWeAreUpTo, column: columnWeAreUpTo };
}

function formatExpected(expected: string[]) {
    if (expected.length === 1) {
        return expected[0];
    }
    return 'one of ' + expected.join(', ');
}

function formatGot(input: string, error: MonadicParser.ParseFailure) {
    const index = error.index;
    const i = index.offset;
    if (i === input.length) {
        return ', got the end of the input';
    }
    const prefix = i > 0 ? '\'...' : '\'';
    const suffix = input.length - i > 12 ? '...\'' : '\'';
    return ` at line ${index.line} column ${index.column}, got ${prefix}${input.slice(i, i + 12)}${suffix}`;
}

function formatError(input: string, error: MonadicParser.ParseFailure) {
    return `expected ${formatExpected(error.expected)}${formatGot(input, error)}`;
}

function unsafeUnion(xs: string[], ys: string[]) {
    const xn = xs.length;
    const yn = ys.length;
    if (xn === 0) return ys;
    else if (yn === 0) return xs;

    const set = new Set<string>();
    const ret: string[] = [];
    for (let i = 0; i < xn; i++) {
        if (!set.has(xs[i])) {
            ret[ret.length] = xs[i];
            set.add(xs[i]);
        }
    }
    for (let i = 0; i < yn; i++) {
        if (!set.has(ys[i])) {
            ret[ret.length] = ys[i];
            set.add(ys[i]);
        }
    }
    ret.sort();
    return ret;
}

function isParser(obj: any): obj is MonadicParser<any> {
    return obj instanceof MonadicParser;
}