/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { compile } from '@molstar/model/script/runtime/query/compiler';
import { compileScript } from '@molstar/query-language/compile';
import { formatMolScript } from '@molstar/query-language/language/expression-formatter';
import { Expression } from '@molstar/query-language/language/expression';
import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import type { Failure, MolQLSelector, QueryLanguage } from './types.js';

// Ported from Stories Untold's `any-molql` (compile.ts). Failures are values: `{ ok: false, message, hint }`.

export type Compiled = { ok: true; molql: Expression } | Failure;

const QueryHint =
  'Write a keyword (protein, ligand, water, ...), an address list ("A:120-130", "B"), or PyMOL selection text ' +
  '("chain A and resn HIS", "byres polymer within 5 of resn STI").';

const LanguageName: Record<QueryLanguage, string> = {
  pymol: 'PyMOL',
  vmd: 'VMD',
  jmol: 'Jmol',
  molscript: 'MolScript',
};

/**
 * A selection written as text (PyMOL, VMD or Jmol selection syntax, or MolScript) as MolQL, checked by Mol*'s query
 * compiler. Errors are said in words the writer can act on.
 */
export function compileQuery(text: string, language: QueryLanguage = 'pymol'): Compiled {
  const source = text.trim();
  if (!source) return failure('The query is empty.', QueryHint);
  let molql: Expression;
  try {
    molql = compileScript(language === 'molscript' ? 'mol-script' : language, source);
  } catch (e) {
    // Mol*'s PyMOL parser recurses without end on an operator missing its left side ("polymer and within 5 of ...")
    // and on a group opening with two parentheses ("((a) within 5 of b)"); say what to write instead of "Maximum call
    // stack size".
    if (/call stack/i.test(messageOf(e))) {
      return failure(
        `${LanguageName[language]} could not read it: the parser went in circles.`,
        'Give "within" and "around" a left side, and open no group with two parentheses; e.g. "byres (chain A and polymer) within 5 of (resn STI)".',
      );
    }
    return failure(`${LanguageName[language]} could not read it: ${messageOf(e)}`, languageHint(language));
  }
  const checked = checkMolQL(molql);
  return checked === true ? { ok: true, molql } : failure(checked, languageHint(language));
}

function languageHint(language: QueryLanguage): string {
  return language === 'pymol'
    ? 'PyMOL selection syntax, e.g. "chain A and resi 100-120", "resn HIS", "byres polymer within 5 of resn HEM"; or an address list ("A:100-120").'
    : `${LanguageName[language]} selection syntax, or an address list ("A:100-120") or a keyword (protein, ligand, ...).`;
}

/** Whether a MolQL expression compiles (true), or why not. */
export function checkMolQL(molql: unknown): true | string {
  if (!molql || typeof molql !== 'object' || !('head' in molql)) return 'MolQL must be a MolScript application.';
  const issues = expressionValidationIssues(molql);
  if (issues) return `Invalid MolQL: ${messageOf(issues.join('; '))}`;
  try {
    compile(molql as Expression);
    return true;
  } catch (e) {
    return `Mol* cannot run it: ${messageOf(e)}`;
  }
}

/** A MolQL expression as MolScript text, for people ("(structure-query.generator.atom-groups ...)"). */
export function formatMolQL(molql: unknown): string {
  return formatMolScript(molql as Expression);
}

/**
 * A MolQL expression as what is kept (a component, a snapshot, the history, a run's arguments): `{ molql }`, plain
 * JSON, so a session replays without a parser (targets.md §2). A copy: what is stored does not change with the
 * expression.
 */
export function toMolqlJson(molql: Expression): MolQLSelector {
  return { molql: JSON.parse(JSON.stringify(molql)) };
}

export type MolqlFromJson = { ok: true; molql: Expression } | Failure;

/**
 * What was kept, back as MolQL: `{ molql }`, the expression itself, or either as JSON text. Checked as stored values
 * should be (names, arguments, and that Mol* can run it).
 */
export function fromMolqlJson(json: unknown): MolqlFromJson {
  let value = json;
  if (typeof value === 'string') {
    try {
      value = JSON.parse(value);
    } catch (e) {
      return failure(`The stored MolQL is not JSON: ${messageOf(e)}`, 'Store a target as { "molql": <expression> }.');
    }
  }
  if (value && typeof value === 'object' && !Array.isArray(value) && Object.hasOwn(value, 'molql')) {
    value = (value as MolQLSelector).molql;
  }
  if (!Expression.is(value)) {
    return failure('The stored MolQL is not an expression.', 'Store a target as { "molql": <expression> }.');
  }
  const checked = checkMolQL(value);
  return checked === true
    ? { ok: true, molql: value }
    : failure(checked, 'Store a target as { "molql": <expression> }.');
}

export function failure(message: string, hint: string): Failure {
  return { ok: false, message, hint };
}

export function messageOf(e: unknown): string {
  const text = e instanceof Error ? e.message : String(e);
  return text.replace(/\s+/g, ' ').trim().slice(0, 300);
}
