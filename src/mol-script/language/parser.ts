/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MonadicParser as P } from '../../mol-util/monadic-parser';
import Expression from './expression';
import { MolScriptBuilder as B } from './builder';

export function parseMolScript(input: string) {
    return Language.parse(input);
}

namespace Language {
    type AST = ASTNode.Expression[]

    namespace ASTNode {
        export type Expression = Str | Symb | List | Comment

        export interface Str {
            kind: 'string',
            value: string
        }

        export interface Symb {
            kind: 'symbol',
            value: string
        }

        export interface List {
            kind: 'list',
            bracket: '(' | '[' | '{',
            nodes: Expression[]
        }

        export interface Comment {
            kind: 'comment',
            value: string
        }

        export function str(value: string): Str { return { kind: 'string', value }; }
        export function symb(value: string): Symb { return { kind: 'symbol', value }; }
        export function list(bracket: '(' | '[' | '{', nodes: Expression[]): List { return { kind: 'list', bracket, nodes }; }
        export function comment(value: string): Comment { return { kind: 'comment', value }; }
    }

    const ws = P.regexp(/[\n\r\s]*/);
    const Expr: P<ASTNode.Expression> = P.lazy(() => (P.alt(Str, List, Symb, Comment).trim(ws)));
    const Str = P.takeWhile(c => c !== '`').trim('`').map(ASTNode.str);
    const Symb = P.regexp(/[^()\[\]{};`,\n\r\s]+/).map(ASTNode.symb);
    const Comment = P.regexp(/\s*;+([^\n\r]*)\n/, 1).map(ASTNode.comment);
    const Args = Expr.many();
    const List1 = Args.wrap('(', ')').map(args => ASTNode.list('(', args));
    const List2 = Args.wrap('[', ']').map(args => ASTNode.list('[', args));
    const List3 = Args.wrap('{', '}').map(args => ASTNode.list('{', args));
    const List = P.alt(List1, List2, List3);

    const Expressions: P<AST> = Expr.many();

    function getAST(input: string) { return Expressions.tryParse(input); }

    function visitExpr(expr: ASTNode.Expression): Expression {
        switch (expr.kind) {
            case 'string': return expr.value;
            case 'symbol': {
                const value = expr.value;
                if (value.length > 1) {
                    const fst = value.charAt(0);
                    switch (fst) {
                        case '.': return B.atomName(value.substr(1));
                        case '_': return B.struct.type.elementSymbol([value.substr(1)]);
                    }
                }
                if (value === 'true') return true;
                if (value === 'false') return false;
                if (isNumber(value)) return +value;
                return Expression.Symbol(value);
            }
            case 'list': {
                switch (expr.bracket) {
                    case '[': return B.core.type.list(withoutComments(expr.nodes).map(visitExpr));
                    case '{': return B.core.type.set(withoutComments(expr.nodes).map(visitExpr));
                    case '(': {
                        const head = visitExpr(expr.nodes[0]);
                        return Expression.Apply(head, getArgs(expr.nodes));
                    }
                }
                return 0 as any;
            }
            default: {
                throw new Error('should not happen');
            }
        }
    }

    function getArgs(nodes: ASTNode.Expression[]): Expression.Arguments | undefined {
        if (nodes.length <= 1) return void 0;
        if (!hasNamedArgs(nodes)) {
            const args: Expression[] = [];
            for (let i = 1, _i = nodes.length; i < _i; i++) {
                const n = nodes[i];
                if (n.kind === 'comment') continue;
                args[args.length] = visitExpr(n);
            }
            return args;
        }
        const args: { [name: string]: Expression } = {};
        let allNumeric = true;
        let pos = 0;
        for (let i = 1, _i = nodes.length; i < _i; i++) {
            const n = nodes[i];
            if (n.kind === 'comment') continue;
            if (n.kind === 'symbol' && n.value.length > 1 && n.value.charAt(0) === ':') {
                const name = n.value.substr(1);
                ++i;
                while (i < _i && nodes[i].kind === 'comment') { i++; }
                if (i >= _i) throw new Error(`There must be a value foolowed a named arg ':${name}'.`);
                args[name] = visitExpr(nodes[i]);
                if (isNaN(+name)) allNumeric = false;
            } else {
                args[pos++] = visitExpr(n);
            }
        }
        if (allNumeric) {
            const keys = Object.keys(args).map(a => +a).sort((a, b) => a - b);
            let isArray = true;
            for (let i = 0, _i = keys.length; i < _i; i++) {
                if (keys[i] !== i) {
                    isArray = false;
                    break;
                }
            }
            if (isArray) {
                const arrayArgs: Expression[] = [];
                for (let i = 0, _i = keys.length; i < _i; i++) {
                    arrayArgs[i] = args[i];
                }
                return arrayArgs;
            }
        }
        return args;
    }

    function hasNamedArgs(nodes: ASTNode.Expression[]) {
        for (let i = 1, _i = nodes.length; i < _i; i++) {
            const n = nodes[i];
            if (n.kind === 'symbol' && n.value.length > 1 && n.value.charAt(0) === ':') return true;
        }
        return false;
    }

    function withoutComments(nodes: ASTNode.Expression[]) {
        let hasComment = false;
        for (let i = 0, _i = nodes.length; i < _i; i++) {
            if (nodes[i].kind === 'comment') {
                hasComment = true;
                break;
            }
        }
        if (!hasComment) return nodes;
        return nodes.filter(n => n.kind !== 'comment');
    }

    function isNumber(value: string) {
        return /-?(0|[1-9][0-9]*)([.][0-9]+)?([eE][+-]?[0-9]+)?/.test(value) && !isNaN(+value);
    }

    export function parse(input: string): Expression[] {
        const ast = getAST(input);
        const ret: Expression[] = [];
        for (const expr of ast) {
            if (expr.kind === 'comment') continue;
            ret[ret.length] = visitExpr(expr);
        }
        return ret;
    }
}