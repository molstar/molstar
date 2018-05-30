/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexanderose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MonadicParser as P } from 'mol-util/monadic-parser'

import Parser from '../parser'
import Expression from '../../expression'
import { SymbolMap, MolScriptSymbol } from './symbols'
import B from '../../builder'

const ws = P.regexp(/[\n\r\s]*/)

function getSymbolExpression(s: MolScriptSymbol, args?: any) {
    switch (s.kind) {
        case 'alias': return args ? Expression.Apply(s.symbol.id, args) : Expression.Apply(s.symbol.id);
        case 'macro': return s.translate(args);
    }
}

namespace Language {
    const Expr = P.lazy(() => P.seq(Symb, ArgList, NamedArgList));

    const Arg: P<Expression> = P.lazy(() => P.seq(
        P.lookahead(P.regexp(/[^:]/)),
        P.alt(
            // order matters
            AtomName,
            ElementSymbol,
            Bool,
            Num,
            Str,
            QuotedStr,
            ListSymbol,
            SetSymbol,
            List
        )
    ).map((x: any) => x[1]).trim(ws));

    const ArgList = Arg.many();
    const ArgName = P.regexp(/:([a-zA-Z0-9_.-]+)/, 1).trim(ws).desc('arg-name');

    const NamedArg = P.seq(ArgName, Arg).trim(ws);

    const NamedArgList = NamedArg.many().map(xs => {
        const namedArgs: { [key: string]: any } = {}
        xs.forEach((a: any) => { namedArgs[a[0]] = a[1] })
        return namedArgs
    });

    const Symb = P.regexp(/[^\s'`,@()\[\]{}';:]+/)  // /[a-zA-Z_-][a-zA-Z0-9_.-]+/)
        .map(x => {
            const s = SymbolMap[x];
            if (!s) {
                throw new Error(`'${x}': unknown symbol.`);
            }
            return s;
        })
        .desc('symbol');

    const Str = P.regexp(/[a-zA-Z_-]+[a-zA-Z0-9_.-]*/).map(x => {
        const s = SymbolMap[x];
        if (s) return getSymbolExpression(s);
        return x;
    }).desc('string');

    const QuotedStr = P.string('`')
        .then(P.regexp(/[^`]*/))
        .skip(P.string('`'))
        .desc('quoted-string');

    const Num = P.regexp(/-?(0|[1-9][0-9]*)([.][0-9]+)?([eE][+-]?[0-9]+)?/)
        .map(v => +v)
        .desc('number');

    const Bool = P.alt(
        P.regexp(/true/i).result(true),
        P.regexp(/false/i).result(false)
    ).desc('boolean');

    // '[a, b, c]' => core.list([a, b, c])
    const ListSymbol = ArgList
        .wrap(P.string('['), P.string(']'))
        .map(B.core.type.list)
        .desc('list-symbol');

    // '{a, b, c}' => core.set([a, b, c])
    const SetSymbol = ArgList
        .wrap(P.string('{'), P.string('}'))
        .map(B.core.type.set)
        .desc('set-symbol');

    // _XYZ -> type.elementSymbol XYZ
    const ElementSymbol = P.string('_')
        .then(P.regexp(/[0-9a-zA-Z]+/))
        .map(x => B.struct.type.elementSymbol([x]))
        .desc('element-symbol');

    // '.e' => struct.type.atomName(e)
    const AtomName = P.string('.')
        .then(P.alt(Str, QuotedStr, Num))
        .map(v => B.atomName('' + v))
        .desc('identifier');

    const List = Expr
        .wrap(P.string('('), P.string(')'))
        .map(x => {
            const array: any[] = x[1];
            const named: any = x[2];

            if (named && Object.keys(named).length) {
                if (array) {
                    for (let i = 0; i < array.length; i++) named[i] = array[i];
                }
                return getSymbolExpression(x[0], named);
            } else if (array && array.length) {
                return getSymbolExpression(x[0], x[1]);
            } else {
                return getSymbolExpression(x[0])
            }
        })
        .desc('list');

    export const Query = List.trim(ws)
}

const reComment = /;[^\n\r]*[\n\r]/g
const transpiler: Parser = str => Language.Query.tryParse(str.replace(reComment, '\n'))
export default transpiler
