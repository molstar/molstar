/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Type from '../type';
import { MSymbol, Arguments, Argument } from '../symbol';
import { symbol, normalizeTable, symbolList } from '../helpers';

export namespace Types {
    export type List<T = any> = ArrayLike<T>
    export type Set<T = any> = { has(e: T): boolean }

    export const AnyVar = Type.Variable('a', Type.Any);
    export const AnyValueVar = Type.Variable('a', Type.Any);
    export const ConstrainedVar = Type.Variable('a', Type.Any, true);

    export const Regex = Type.Value<RegExp>('Core', 'Regex');

    export const Set = <T extends Type>(t?: T) => Type.Container<Set<T['@type']>>('Core', 'Set', t || AnyValueVar);
    export const List = <T extends Type>(t?: T) => Type.Container<List<T['@type']>>('Core', 'List', t || AnyVar);
    export const Fn = <T extends Type>(t?: T, alias?: string) => Type.Container<(env: any) => T['@type']>('Core', 'Fn', t || AnyVar, alias);
    export const Flags = <T extends Type>(t: T, alias?: string) => Type.Container<number>('Core', 'Flags', t, alias);

    export const BitFlags = Flags(Type.Num, 'BitFlags');
}

function unaryOp<T extends Type>(type: T, description?: string) {
    return symbol(Arguments.Dictionary({ 0: Argument(type) }), type, description);
}

function binOp<T extends Type>(type: T, description?: string) {
    return symbol(Arguments.List(type, { nonEmpty: true }), type, description);
}

function binRel<A extends Type, T extends Type>(src: A, target: T, description?: string) {
    return symbol(Arguments.Dictionary({
        0: Argument(src),
        1: Argument(src)
    }), target, description);
}

export const TTargs = Arguments.Dictionary({
    0: Argument(Type.Num),
    1: Argument(Type.Num)
});

const XX = { test: Argument(Type.Str) };
const t: Arguments.PropTypes<typeof XX> = 0 as any;
t.test;

const type = {
    '@header': 'Types',
    bool: symbol(Arguments.Dictionary({ 0: Argument(Type.AnyValue) }), Type.Bool, 'Convert a value to boolean.'),
    num: symbol(Arguments.Dictionary({ 0: Argument(Type.AnyValue) }), Type.Num, 'Convert a value to number.'),
    str: symbol(Arguments.Dictionary({ 0: Argument(Type.AnyValue) }), Type.Str, 'Convert a value to string.'),
    regex: symbol(
        Arguments.Dictionary({
            0: Argument(Type.Str, { description: 'Expression' }),
            1: Argument(Type.Str, { isOptional: true, description: `Flags, e.g. 'i' for ignore case` })
        }), Types.Regex, 'Creates a regular expression from a string using the ECMAscript syntax.'),

    list: symbol(Arguments.List(Types.AnyVar), Types.List()),
    set: symbol(Arguments.List(Types.AnyValueVar), Types.Set()),
    bitflags: symbol(Arguments.Dictionary({ 0: Argument(Type.Num) }), Types.BitFlags, 'Interpret a number as bitflags.'),
    compositeKey: symbol(Arguments.List(Type.AnyValue), Type.AnyValue),
};

const logic = {
    '@header': 'Logic',
    not: unaryOp(Type.Bool),
    and: binOp(Type.Bool),
    or: binOp(Type.Bool),
};

const ctrl = {
    '@header': 'Control',
    eval: symbol(Arguments.Dictionary({ 0: Argument(Types.Fn(Types.AnyVar)) }), Types.AnyVar, 'Evaluate a function.'),
    fn: symbol(Arguments.Dictionary({ 0: Argument(Types.AnyVar) }), Types.Fn(Types.AnyVar), 'Wrap an expression to a "lazy" function.'),
    if: symbol(Arguments.Dictionary({
        0: Argument(Type.Bool, { description: 'Condition' }),
        1: Argument(Type.Variable('a', Type.Any), { description: 'If true' }),
        2: Argument(Type.Variable('b', Type.Any), { description: 'If false' })
    }), Type.Union([Type.Variable('a', Type.Any), Type.Variable('b', Type.Any)])),
    assoc: symbol(Arguments.Dictionary({
        0: Argument(Type.Str, { description: 'Name' }),
        1: Argument(Type.Variable('a', Type.Any), {description: 'Value to assign' })
    }), Type.Variable('a', Type.Any))
};

const rel = {
    '@header': 'Relational',
    eq: binRel(Type.Variable('a', Type.AnyValue, true), Type.Bool),
    neq: binRel(Type.Variable('a', Type.AnyValue, true), Type.Bool),
    lt: binRel(Type.Num, Type.Bool),
    lte: binRel(Type.Num, Type.Bool),
    gr: binRel(Type.Num, Type.Bool),
    gre: binRel(Type.Num, Type.Bool),
    inRange: symbol(Arguments.Dictionary({
        0: Argument(Type.Num, { description: 'Value to test' }),
        1: Argument(Type.Num, { description: 'Minimum value' }),
        2: Argument(Type.Num, { description: 'Maximum value' })
    }), Type.Bool, 'Check if the value of the 1st argument is >= 2nd and <= 3rd.'),
};

const math = {
    '@header': 'Math',
    add: binOp(Type.Num),
    sub: binOp(Type.Num),
    mult: binOp(Type.Num),
    div: binRel(Type.Num, Type.Num),
    pow: binRel(Type.Num, Type.Num),
    mod: binRel(Type.Num, Type.Num),

    min: binOp(Type.Num),
    max: binOp(Type.Num),

    floor: unaryOp(Type.Num),
    ceil: unaryOp(Type.Num),
    roundInt: unaryOp(Type.Num),
    abs: unaryOp(Type.Num),
    sqrt: unaryOp(Type.Num),
    cbrt: unaryOp(Type.Num),
    sin: unaryOp(Type.Num),
    cos: unaryOp(Type.Num),
    tan: unaryOp(Type.Num),
    asin: unaryOp(Type.Num),
    acos: unaryOp(Type.Num),
    atan: unaryOp(Type.Num),
    sinh: unaryOp(Type.Num),
    cosh: unaryOp(Type.Num),
    tanh: unaryOp(Type.Num),
    exp: unaryOp(Type.Num),
    log: unaryOp(Type.Num),
    log10: unaryOp(Type.Num),
    atan2: binRel(Type.Num, Type.Num)
};

const str = {
    '@header': 'Strings',
    concat: binOp(Type.Str),
    match: symbol(Arguments.Dictionary({ 0: Argument(Types.Regex), 1: Argument(Type.Str) }), Type.Bool)
};

const list = {
    '@header': 'Lists',
    getAt: symbol(Arguments.Dictionary({ 0: Argument(Types.List()), 1: Argument(Type.Num) }), Types.AnyVar),
    equal: symbol(Arguments.Dictionary({ 0: Argument(Types.List()), 1: Argument(Types.List()) }), Type.Bool)
};

const set = {
    '@header': 'Sets',
    has: symbol(Arguments.Dictionary({ 0: Argument(Types.Set(Types.ConstrainedVar)), 1: Argument(Types.ConstrainedVar) }), Type.Bool, 'Check if the the 1st argument includes the value of the 2nd.'),
    isSubset: symbol(Arguments.Dictionary({ 0: Argument(Types.Set(Types.ConstrainedVar)), 1: Argument(Types.Set(Types.ConstrainedVar)) }), Type.Bool, 'Check if the the 1st argument is a subset of the 2nd.')
};

const flags = {
    '@header': 'Flags',
    hasAny: symbol(Arguments.Dictionary({
        0: Argument(Types.Flags(Types.ConstrainedVar)),
        1: Argument(Types.Flags(Types.ConstrainedVar))
    }), Type.Bool, 'Check if the the 1st argument has at least one of the 2nd one\'s flags.'),
    hasAll: symbol(Arguments.Dictionary({
        0: Argument(Types.Flags(Types.ConstrainedVar)),
        1: Argument(Types.Flags(Types.ConstrainedVar))
    }), Type.Bool, 'Check if the the 1st argument has all 2nd one\'s flags.'),
};

const table = {
    '@header': 'Language Primitives',
    type,
    logic,
    ctrl,
    rel,
    math,
    str,
    list,
    set,
    flags
};

normalizeTable(table);

export const SymbolList = symbolList(table);

export const SymbolMap = (function() {
    const map: { [id: string]: MSymbol | undefined } = Object.create(null);
    for (const s of SymbolList) map[s.id] = s;
    return map;
})();

export default table;
