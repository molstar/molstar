/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Extensions of the query language (_spec/README.md "Queries", Extensions; extension-openings.md; properties.md): a
// `query-form` registry entry adds a form of a role (a generator, filter, modifier or combinator over results, or an
// operator over values), a `query-property` entry a property every test and filter reads and `elements` takes as a
// field ("basically just adding a few symbols", the author, 2026-11-20). Checked when composed: names `<id>:<name>`
// (aliases too), no built-in shadowed, every conflict across the entries listed in one error. Checked by role when
// evaluated: a filter's result a subset of its prefix, of its type; a combinator's operands one type.

import { KindPattern, type Node, type Value } from '@molstar/core/repr/tree';
import { ExpressionKinds, type ExprValue, type QueryOperator } from './expressions.js';
import type { QueryModel, PropertyRead } from './model.js';
import {
  BuiltInProperties,
  type PropertyDef,
  type PropertyLevel,
  type PropertyTable,
  type PropertyType,
  propertyTable,
} from './properties.js';
import type { QueryResult } from './result.js';

export type QueryFormRole = 'generator' | 'filter' | 'modifier' | 'combinator' | 'operator';

/** What an extension's form is given: the model, its form, what its prefix names, and the evaluation's bounds. */
export interface QueryFormContext {
  readonly model: QueryModel;
  readonly form: Node;
  /** What the prefix names, as a result. */
  readonly prefix: QueryResult;
  /** A nested value's result: a form against the same prefix, an address from the root. */
  evaluate(value: Value): QueryResult;
  /** Counts a visit: the evaluation's signal and deadline are read every so many. */
  tick(): void;
  /** The most matches the result may hold: more is an issue, never a truncation. */
  readonly limit: number;
}

/** A form over results (`af:confident`): its role, its keyword params in their order, its result's type. */
export interface QueryFormEntry {
  readonly kind: 'query-form';
  readonly name: string;
  readonly role: Exclude<QueryFormRole, 'operator'>;
  readonly params: readonly string[];
  readonly result: 'location' | 'keys' | 'nodes';
  resolve(ctx: QueryFormContext): QueryResult;
}

/** An operator over values (`af:near`): its result's type; operands positional, as the built-ins'. */
export interface QueryOperatorEntry {
  readonly kind: 'query-form';
  readonly name: string;
  readonly role: 'operator';
  readonly result: 'value';
  /** Its value's type (`boolean` for a test). */
  readonly type: PropertyType;
  resolve(args: readonly ExprValue[]): ReturnType<QueryOperator['evaluate']>;
}

/** A property (`af:plddt`): its level and type, read by a function over the model or a values annotation's field. */
export interface QueryPropertyEntry {
  readonly kind: 'query-property';
  readonly name: string;
  readonly level: PropertyLevel;
  readonly type: PropertyType;
  readonly aliases?: readonly string[];
  readonly read:
    | ((model: QueryModel, system: string) => PropertyRead)
    | { readonly values: string; readonly field: string };
}

export type QueryExtensionEntry = QueryFormEntry | QueryOperatorEntry | QueryPropertyEntry;

/** The registry's query entries, composed: forms by name, operators by name, the property table. */
export interface QueryExtensions {
  readonly forms: ReadonlyMap<string, QueryFormEntry>;
  readonly operators: ReadonlyMap<string, QueryOperator>;
  readonly properties: PropertyTable;
}

export const NoQueryExtensions: QueryExtensions = Object.freeze({
  forms: new Map(),
  operators: new Map(),
  properties: propertyTable([]),
});

/** A namespaced name: `<id>:<name>`, both words. */
const isNamespaced = (n: string) => n.includes(':') && KindPattern.test(n);

/**
 * The query entries of a registry composed (`seq.context.extension.namespace` and every conflict in one error): a
 * name not namespaced, one a built-in has, or two different entries under one name. The same entry twice is one.
 */
export function queryExtensions(entries: readonly QueryExtensionEntry[], builtIns: readonly string[]): QueryExtensions {
  if (!entries.length) return NoQueryExtensions;
  const problems: string[] = [];
  const taken = new Map<string, QueryExtensionEntry>();
  const reserved = new Set([
    ...builtIns,
    ...ExpressionKinds,
    ...BuiltInProperties.flatMap((d) => [d.name, ...d.aliases]),
  ]);
  for (const e of new Set(entries)) {
    for (const n of [e.name, ...(e.kind === 'query-property' ? (e.aliases ?? []) : [])]) {
      if (!isNamespaced(n)) problems.push(`${n} is not namespaced (<id>:<name>): seq.context.extension.namespace`);
      else if (reserved.has(n)) problems.push(`${n} shadows a built-in`);
      else if (taken.has(n) && taken.get(n) !== e) problems.push(`${n} is registered twice`);
      taken.set(n, e);
    }
  }
  if (problems.length) throw new Error(`The registry's query entries conflict: ${problems.join('; ')}`);
  const forms = new Map<string, QueryFormEntry>();
  const operators = new Map<string, QueryOperator>();
  const props: PropertyDef[] = [];
  for (const e of new Set(entries)) {
    if (e.kind === 'query-property')
      props.push({ name: e.name, aliases: e.aliases ?? [], level: e.level, type: e.type, read: e.read });
    else if (e.role === 'operator')
      operators.set(e.name, { name: e.name, type: e.type, evaluate: (a) => e.resolve(a) });
    else forms.set(e.name, e);
  }
  return { forms, operators, properties: propertyTable(props) };
}
