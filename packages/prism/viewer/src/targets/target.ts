/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  SchemaItemFields,
  hasWrongTypes,
  isSchemaItemField,
} from '@molstar/model/model/structure/structure/element/schema-selection/fields';
import type { Expression } from '@molstar/query-language/language/expression';
import { type AddressResult, parseAddressList } from './addresses.js';
import { LanguageName, checkMolQL, compileQuery, failure, formatMolQL, messageOf } from './compile.js';
import { KEYWORD_EXPRESSIONS, isKeyword, readKeyword } from './keywords.js';
import { defaultQueryLanguage } from './languages.js';
import { isMolQLSelector, itemsExpression } from './selector.js';
import type {
  CompileOptions,
  CompiledTarget,
  ComponentExpression,
  Failure,
  ParseResult,
  ParsedTarget,
  SelectorT,
  Target,
} from './types.js';

// Reading a target (targets.md §1): a string is read in order, a keyword, then an address list, then a query text in the
// given language; a selector (schema items, `{ molql }`) is checked. Everything compiles to MolQL (§2).

/** A read target, with the MolQL of a query text (read once, whoever asks). */
type Read = { ok: true; parsed: ParsedTarget; molql?: Expression } | Failure;

/**
 * Reads a target. Does not resolve `selection` and `focus` (the viewer does; see `CompileOptions.resolvePerson`).
 * A string that is none of keyword, addresses or query text fails with the part that did not read and a hint.
 */
export function parseTarget(target: Target, options: Pick<CompileOptions, 'language'> = {}): ParseResult {
  const read = readTarget(target, options);
  return read.ok ? { ok: true, parsed: read.parsed } : read;
}

function readTarget(target: Target, options: Pick<CompileOptions, 'language'>): Read {
  if (typeof target === 'string') return parseString(target, options.language ?? defaultQueryLanguage());
  return parseSelector(target);
}

function parseString(text: string, language: NonNullable<CompileOptions['language']>): Read {
  if (!text.trim()) {
    return failure(
      'The target is empty.',
      'Name what you mean: a keyword (protein, ligand, selection, ...), an address list ("A:120-130") or PyMOL text.',
    );
  }
  const word = readKeyword(text);
  if (word !== undefined) {
    return { ok: true, parsed: isKeyword(word) ? { kind: 'keyword', keyword: word } : { kind: 'person', which: word } };
  }
  const addresses: AddressResult = parseAddressList(text);
  if (addresses.ok) return { ok: true, parsed: { kind: 'addresses', items: addresses.items } };
  const query = compileQuery(text, language);
  if (query.ok) {
    return { ok: true, parsed: { kind: 'query', language, text }, molql: query.molql };
  }
  // Neither: say what the text looked like.
  return addresses.addressy
    ? failure(addresses.message, `${addresses.hint} (Not read as ${LanguageName[language]} either: ${query.message})`)
    : query;
}

function parseSelector(selector: Exclude<Target, string>): Read {
  if (isMolQLSelector(selector)) {
    const checked = checkMolQL(selector.molql);
    return checked === true
      ? { ok: true, parsed: { kind: 'selector', selector } }
      : failure(checked, 'A { molql } selector holds a MolScript expression as JSON.');
  }
  const rows = Array.isArray(selector) ? (selector as readonly unknown[]) : [selector];
  if (rows.length === 0) {
    return failure(
      'The list of selectors is empty.',
      'Give at least one, e.g. { auth_asym_id: "A", auth_seq_id: 120 }.',
    );
  }
  for (const row of rows) {
    const problem = rowProblem(row);
    if (problem) return failure(problem, `The fields are ${SchemaItemFields.join(', ')}.`);
  }
  return { ok: true, parsed: { kind: 'selector', selector: selector as Exclude<SelectorT, string> } };
}

/** Why a row is not a component expression, if it is not. */
function rowProblem(row: unknown): string | undefined {
  if (!row || typeof row !== 'object' || Array.isArray(row))
    return 'A selector is an object of fields, e.g. { auth_asym_id: "A" }.';
  const keys = Object.keys(row);
  if (keys.length === 0) return 'A selector without fields would select everything; use the keyword "all".';
  const unknown = keys.find((k) => !isSchemaItemField(k));
  if (unknown !== undefined) return `"${unknown}" is not a field of a selector.`;
  if (hasWrongTypes(row as ComponentExpression)) {
    return 'A field has the wrong type: identifiers are strings ("A"), numbers are numbers (120).';
  }
  return undefined;
}

/**
 * A target as MolQL. `selection` and `focus` are the person's: `options.resolvePerson` gives their expressions; without
 * one they fail with a message. Never throws.
 */
export function compileTarget(target: Target, options: CompileOptions = {}): CompiledTarget {
  const read = readTarget(target, options);
  if (!read.ok) return read;
  const { parsed } = read;
  switch (parsed.kind) {
    case 'keyword':
      return { ok: true, parsed, molql: KEYWORD_EXPRESSIONS[parsed.keyword] };
    case 'person': {
      if (!options.resolvePerson) {
        return failure(
          `"${parsed.which}" is the person's current ${parsed.which}, which only the viewer knows.`,
          'Name what you mean instead ("A:120-130", "ligand"), or ask the viewer to resolve it.',
        );
      }
      const molql = options.resolvePerson(parsed.which);
      if (!molql) {
        return failure(
          parsed.which === 'selection' ? 'Nothing is selected.' : 'Nothing is focused.',
          // Read by people and the agent alike: the agent's tools add their own words (agent.md §2).
          `${parsed.which === 'selection' ? 'Select' : 'Focus'} something first, or name what to use.`,
        );
      }
      return { ok: true, parsed, molql };
    }
    case 'query':
      return { ok: true, parsed, molql: read.molql as Expression };
    case 'addresses':
      return expressionOf(parsed, () => itemsExpression(parsed.items));
    case 'selector': {
      const s = parsed.selector;
      return expressionOf(parsed, () =>
        isMolQLSelector(s)
          ? (s.molql as Expression)
          : itemsExpression(Array.isArray(s) ? (s as readonly ComponentExpression[]) : [s as ComponentExpression]),
      );
    }
  }
}

function expressionOf(parsed: ParsedTarget, make: () => Expression): CompiledTarget {
  try {
    return { ok: true, parsed, molql: make() };
  } catch (e) {
    return failure(`Could not make a query of it: ${messageOf(e)}`, 'Check the values of the fields.');
  }
}

/**
 * A target as MolScript text, for showing a stored one ("(structure-query.generator.atom-groups ...)"). A target that
 * does not compile is shown as it was written.
 */
export function formatTarget(target: Target, options: CompileOptions = {}): string {
  const compiled = compileTarget(target, options);
  if (compiled.ok) return formatMolQL(compiled.molql);
  return typeof target === 'string' ? target : JSON.stringify(target);
}

export type { Failure };
