/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MSymbol } from './symbol.js';
import type { Type } from './type.js';
import { SymbolMap } from './symbol-table.js';

export interface ExpressionValidationOptions {
  /** Resolve symbol definitions, including any application-specific vocabulary. Defaults to the standard symbol table. */
  getSymbol?: (name: string) => MSymbol | undefined;
  /** Only check expression syntax, without resolving callable symbols or checking their argument names. */
  syntaxOnly?: boolean;
}

/** Validate expression syntax and argument names without compiling or evaluating a molecular query.
 * Does not check required arguments, argument types, or result types. */
export function expressionValidationIssues(
  expression: unknown,
  options: ExpressionValidationOptions = {},
): string[] | undefined {
  const issues: string[] = [];
  const active = new WeakSet<object>();
  const getSymbol = options.getSymbol ?? ((name: string) => SymbolMap[name]);
  const fieldPath = (path: string, key: string) => `${path}[${JSON.stringify(key)}]`;
  const report = (path: string, message: string) => issues.push(`${path}: ${message}`);

  /** `argument`: the value is an argument of an application, where `null` and a list of values are literals if the argument can be a list
   * (`admitsList`, unknown for a symbol which is not looked up) */
  function visit(value: unknown, path: string, argument?: { admitsList: boolean }) {
    if (typeof value === 'string' || typeof value === 'boolean') return;
    if (typeof value === 'number' && Number.isFinite(value)) return;
    if (argument && (value === null || Array.isArray(value)) && !argument.admitsList) {
      report(path, 'This argument is not a list of values, or a missing value (null).');
      return;
    }
    if (argument && value === null) return;
    if (argument && Array.isArray(value)) {
      // (a list of values, e.g. a column of a table)
      const bad = value.findIndex(
        (v) =>
          !(
            v === null ||
            typeof v === 'string' ||
            typeof v === 'boolean' ||
            (typeof v === 'number' && Number.isFinite(v))
          ),
      );
      if (bad >= 0)
        report(`${path}[${bad}]`, 'Expected a string, a finite number, a boolean or null in a list of values.');
      return;
    }
    if (!isObject(value)) {
      report(path, 'Expected a literal, symbol, or application.');
      return;
    }
    if (active.has(value)) {
      report(path, 'Cyclic expression.');
      return;
    }
    active.add(value);
    if (Object.hasOwn(value, 'name')) {
      if (typeof value.name !== 'string') report(`${path}.name`, 'Symbol name must be a string.');
      for (const key of Object.keys(value)) {
        if (key !== 'name') report(fieldPath(path, key), `Unknown symbol field '${key}'.`);
      }
    } else {
      for (const key of Object.keys(value)) {
        if (key !== 'head' && key !== 'args') report(fieldPath(path, key), `Unknown application field '${key}'.`);
      }
      let symbol: MSymbol | undefined;
      if (!isObject(value.head) || typeof value.head.name !== 'string' || Object.hasOwn(value.head, 'head')) {
        report(`${path}.head`, 'Application head must be a symbol.');
      } else {
        visit(value.head, `${path}.head`);
        if (!options.syntaxOnly) {
          symbol = getSymbol(value.head.name);
          if (!symbol) report(`${path}.head.name`, `Unknown callable symbol '${value.head.name}'.`);
        }
      }
      if (value.args !== undefined) {
        if (!Array.isArray(value.args) && !isObject(value.args)) {
          report(`${path}.args`, 'Arguments must be an array or an object.');
        } else {
          const keys = Array.isArray(value.args)
            ? [...new Set([...value.args.keys()].map(String).concat(Object.keys(value.args)))]
            : Object.keys(value.args);
          for (const key of keys) {
            const argumentPath = fieldPath(`${path}.args`, key);
            if (Array.isArray(value.args) && !/^(0|[1-9]\d*)$/.test(key)) {
              report(argumentPath, `Unknown array field '${key}'.`);
            }
            if (symbol) {
              const validKey =
                symbol.args.kind === 'dictionary' ? Object.hasOwn(symbol.args.map, key) : /^(0|[1-9]\d*)$/.test(key);
              if (!validKey) report(argumentPath, `Unknown argument '${key}' for symbol '${symbol.id}'.`);
            }
            visit((value.args as Record<string, unknown>)[key], argumentPath, {
              admitsList: !symbol || argumentAdmitsList(symbol, key),
            });
          }
        }
      }
    }
    active.delete(value);
  }

  visit(expression, 'expression');
  return issues.length ? issues : undefined;
}

/** Can the argument `key` of `symbol` be a list of values (a column of a table) or a missing value? Only if it is declared to be a list (or of any type). */
function argumentAdmitsList(symbol: MSymbol, key: string): boolean {
  const type =
    symbol.args.kind === 'dictionary'
      ? (symbol.args.map as Record<string, { type: Type } | undefined>)[key]?.type
      : symbol.args.type;
  return type === undefined || canBeList(type);
}

function canBeList(type: Type): boolean {
  switch (type.kind) {
    case 'any':
    case 'any-value':
    case 'variable':
      return true;
    case 'container':
      return type.name === 'List' || type.name === 'Set';
    case 'union':
      return type.types.some(canBeList);
    default:
      return false;
  }
}

function isObject(value: unknown): value is Record<string, unknown> {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
}
