/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A language's own consistency (_spec/README.md "Schema"): `defineLanguage` throws on its issues, a broken description
// written in code being a programmer error. Each issue carries the address of what it is about in the description
// (`/sequence/alphabet`: a kind, then its params by name).

import {
  ColumnTypes,
  defineIssues,
  exactAddress,
  IdPattern,
  type Issue,
  isFloatType,
  isNode,
  isReservedKey,
  issueOf,
  KindHint,
  KindPattern,
  own,
  printAddress,
} from '../tree/index.js';
import type { Fields, Language, ParamSpec, ParamType, UnitsHook } from './language.js';
import { validateValue } from './validate.js';

/** Whether a value of the type may be a node of another language (its default is checked when a spec is). */
const holdsSlot = (t: ParamType): boolean => t.type === 'node' || (t.type === 'union' && t.of.some(holdsSlot));

export const LanguageIssues = defineIssues('repr.schema.language', {
  'repr.schema.language.param': {
    severity: 'error',
    description: 'A param that is both required and defaulted',
  },
  'repr.schema.language.root': { severity: 'error', description: 'The root names no kind of the language' },
  'repr.schema.language.parent': { severity: 'error', description: 'A parent names no kind of the language' },
  'repr.schema.language.union': {
    severity: 'error',
    description:
      'A union with more than one record-like member (a record, a map or a table): a value would read two ways',
  },
  'repr.schema.language.pattern': { severity: 'error', description: 'A pattern is not a regular expression' },
  'repr.schema.language.name': {
    severity: 'error',
    description:
      'A kind is not a (namespaced) word, a param or the language is not an id, or a param is a reserved key',
  },
  'repr.schema.language.default': { severity: 'error', description: 'A default is not of its param type' },
  'repr.schema.language.range': {
    severity: 'error',
    description: 'A min above its max, a min-length above its max-length, or a repeating data form without its max',
  },
  'repr.schema.language.empty': {
    severity: 'error',
    description: 'A union with no member, an enum with no value, or a data form with no positional value',
  },
  'repr.schema.language.ref': {
    severity: 'error',
    description: "A reference's :in is neither a kind of the language nor an address of ids",
  },
  'repr.schema.language.unit': {
    severity: 'error',
    description: 'A unit that is not an id of the units table (with a units hook), or a unit on an integer',
  },
  'repr.schema.language.column': {
    severity: 'error',
    description: 'A column of no element type, or `missing` on an integer column (its validity list is not built yet)',
  },
  'repr.schema.language.table': {
    severity: 'error',
    description: "A table's column that is neither a column nor a list, or a key naming no column",
  },
});
export type LanguageCode = keyof typeof LanguageIssues.codes;

export interface CheckLanguageOptions {
  /** The units table's tie: with it, a param's `unit` must be one of its ids. */
  readonly units?: UnitsHook;
}

/** Every inconsistency of a language: names, roots and parents, patterns, ranges, defaults, columns, at any depth. */
export function checkLanguage(language: Language, options: CheckLanguageOptions = {}): Issue[] {
  const issues: Issue[] = [];
  const report = (code: LanguageCode, message: string, at: readonly string[], hint?: string) =>
    issues.push(issueOf(LanguageIssues, code, message, { address: addressOf(at), hint }));

  if (!IdPattern.test(language.id)) report('repr.schema.language.name', 'The language id is not an id', []);
  if (!own(language.kinds, language.root)) report('repr.schema.language.root', 'The root names no kind', []);
  for (const [kind, spec] of Object.entries(language.kinds)) {
    if (!KindPattern.test(kind)) report('repr.schema.language.name', 'A kind is not a word', [kind], KindHint);
    for (const p of spec.parents)
      if (!own(language.kinds, p)) report('repr.schema.language.parent', 'A parent names no kind', [kind]);
    fields(spec.params ?? {}, [kind]);
    if (spec.data) {
      const d = spec.data;
      for (const a of d.args) param(a, [kind, a.name]);
      if (!d.args.length) report('repr.schema.language.empty', 'A data form with no positional value', [kind]);
      if (d.rest && d.max === undefined)
        report('repr.schema.language.range', 'A repeating data form needs its max (work is bounded)', [kind]);
      if (d.min !== undefined && d.max !== undefined && d.min > d.max)
        report('repr.schema.language.range', 'min > max', [kind]);
    }
  }
  return issues;

  function fields(specs: Fields, at: readonly string[]) {
    for (const [name, spec] of Object.entries(specs)) {
      const here = [...at, name];
      if (!IdPattern.test(name) || isReservedKey(name))
        report('repr.schema.language.name', 'A param is not an id, or is reserved', here, 'kind, __proto__, $…');
      param(spec, here);
    }
  }

  function param(spec: ParamSpec, at: readonly string[]) {
    type(spec, at);
    if (spec.required && spec.default !== undefined)
      report('repr.schema.language.param', 'A param is required or has a default, not both', at);
    // a node in a slot's default needs its language (resolved late): checked when the spec is validated
    const late = isNode(spec.default) && holdsSlot(spec);
    if (spec.default !== undefined && !late && validateValue(spec, spec.default, options).length)
      report('repr.schema.language.default', 'The default is not of its type', at);
  }

  /** The type's own checks, and its members' (a list's item, a union's members, a record's fields). */
  function type(t: ParamType, at: readonly string[]) {
    switch (t.type) {
      case 'string':
        if (t.enum && !t.enum.length) report('repr.schema.language.empty', 'An enum with no value', at);
        if (t.pattern !== undefined) {
          try {
            new RegExp(t.pattern);
          } catch {
            report('repr.schema.language.pattern', 'Not a regular expression', at);
          }
        }
        return;
      case 'number':
      case 'integer':
        if (t.min !== undefined && t.max !== undefined && t.min > t.max)
          report('repr.schema.language.range', 'min > max', at);
        if (t.unit !== undefined && (t.type === 'integer' || (options.units && !options.units.isUnit(t.unit))))
          report('repr.schema.language.unit', 'Not a unit of the table, or a unit on an integer', at);
        return;
      case 'list':
        if (t['min-length'] !== undefined && t['max-length'] !== undefined && t['min-length'] > t['max-length'])
          report('repr.schema.language.range', 'min-length > max-length', at);
        return type(t.of, at);
      case 'record':
        return fields(t.fields, at);
      case 'map':
        return type(t.of, at);
      case 'union':
        if (!t.of.length) report('repr.schema.language.empty', 'A union with no member', at);
        // a record is read by one member: validation takes the first it fits, so two would read one value two ways
        if (t.of.filter((m) => m.type === 'record' || m.type === 'map' || m.type === 'table').length > 1)
          report('repr.schema.language.union', 'A union with more than one record-like member', at, 'one of them');
        for (const m of t.of) type(m, at);
        return;
      case 'ref':
        if (t.in !== undefined && !own(language.kinds, t.in) && exactAddress(t.in) === undefined)
          report('repr.schema.language.ref', ':in is no kind and no address', at, 'a kind (data) or /views');
        return;
      case 'column':
        if (!ColumnTypes.includes(t.of))
          report('repr.schema.language.column', 'Not an element type', at, ColumnTypes.join(', '));
        else if (t.missing && t.of !== 'category' && !isFloatType(t.of))
          report(
            'repr.schema.language.column',
            'An integer column with missing values is not built yet',
            at,
            'a float column (NaN is missing)',
          );
        if (t['min-length'] !== undefined && t['max-length'] !== undefined && t['min-length'] > t['max-length'])
          report('repr.schema.language.range', 'min-length > max-length', at);
        return;
      case 'table':
        for (const [name, c] of Object.entries(t.columns))
          if (c.type !== 'column' && c.type !== 'list')
            report('repr.schema.language.table', 'A table column is a column or a list', [...at, name]);
        if (t.key !== undefined && !own(t.columns, t.key))
          report('repr.schema.language.table', 'The key names no column', at);
        return fields(t.columns, at);
    }
  }
}

/** The address of a place in the description, or undefined when a name in it is no id (it is reported). */
function addressOf(at: readonly string[]): string | undefined {
  return at.every((s) => IdPattern.test(s)) ? printAddress(at) : undefined;
}

/** A language, checked: a broken description written in code is a programmer error, thrown. */
export function defineLanguage<const L extends Language>(language: L, options: CheckLanguageOptions = {}): L {
  const issues = checkLanguage(language, options);
  if (issues.length) {
    const list = issues.map((i) => `${i.code} at ${i.address ?? '?'}`).join('; ');
    throw new Error(`The language ${language.id} is broken: ${list}`);
  }
  return language;
}
