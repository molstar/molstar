/*
 * Copyright (c) 2018-2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type Expression = Expression.Literal | Expression.Symbol | Expression.Apply;

namespace Expression {
  /** A value of a list of literals: `null` is a missing value (it is what JSON makes of `undefined` in a list) */
  export type ListValue = string | number | boolean | null;
  /** A literal: a value, `null` (a missing value), or a list of values, which the compiler passes through as it is (it does not compile every value of it).
   * A list of values is how a column of a table is given to a symbol which takes tables (`structure.generator.element-schema`) in any size. */
  export type Literal = ListValue | ListValue[];
  export type Symbol = { name: string };
  export type Arguments = Expression[] | { [name: string]: Expression };
  export interface Apply {
    readonly head: Expression;
    readonly args?: Arguments;
  }

  export function Symbol(name: string): Expression.Symbol {
    return { name };
  }
  export function Apply(head: Expression, args?: Arguments): Apply {
    return args ? { head, args } : { head };
  }

  export function isArgumentsArray(e?: Arguments): e is Expression[] {
    return !!e && Array.isArray(e);
  }
  export function isArgumentsMap(e?: Arguments): e is { [name: string]: Expression } {
    return !!e && !Array.isArray(e);
  }
  export function isLiteral(e: Expression): e is Expression.Literal {
    return !isApply(e) && !isSymbol(e);
  }
  export function isApply(e: Expression): e is Expression.Apply {
    return !!e && !!(e as Expression.Apply).head && typeof e === 'object';
  }
  export function isSymbol(e: Expression): e is Expression.Symbol {
    return !!e && typeof (e as any).name === 'string';
  }

  /** Decide if a value has the recursive JSON shape of a MolScript expression. */
  export function is(value: unknown): value is Expression {
    if (typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean') return true;
    if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
    return isInArguments(value);
  }

  /** Like `is`, for a value in the arguments of an application, where `null` and a list of values are literals too */
  function isInArguments(value: unknown): value is Expression {
    if (value === null) return true;
    if (Array.isArray(value)) return value.every(isListValue);
    if (typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean') return true;
    if (!value || typeof value !== 'object') return false;
    if (isSymbol(value as Expression)) return true;
    const apply = value as Apply;
    if (!Object.prototype.hasOwnProperty.call(apply, 'head') || !is(apply.head)) return false;
    if (apply.args === undefined) return true;
    if (Array.isArray(apply.args)) return apply.args.every(isInArguments);
    return !!apply.args && typeof apply.args === 'object' && Object.values(apply.args).every(isInArguments);
  }

  function isListValue(value: unknown): value is ListValue {
    return value === null || typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean';
  }
}

export { Expression };
