/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Expression } from '@molstar/query-language/language/expression';
import { parseMolScript } from '@molstar/query-language/language/parser';
import { SymbolMap as MolQLSymbols } from '@molstar/query-language/language/symbol-table';
import { expressionValidationIssues } from '@molstar/query-language/language/validation';
import { SymbolMap as MolScriptSymbols, transpileMolScript } from '@molstar/query-language/script/mol-script/symbols';
import type { QueryLanguageEntry } from '../languages.js';

// MolScript text: one expression, its aliases and macros checked before they are expanded (as `compileScript` does).
function transpile(source: string): Expression {
  const parsed = parseMolScript(source);
  if (parsed.length === 0) throw new Error('No query');
  if (parsed.length !== 1) throw new Error('Expected exactly one MolScript expression.');
  const issues = expressionValidationIssues(parsed[0], {
    getSymbol: (name) => MolScriptSymbols[name]?.symbol ?? MolQLSymbols[name],
  });
  if (issues) throw new Error(`Invalid MolScript expression:\n${issues.join('\n')}`);
  return transpileMolScript(parsed[0]);
}

export const MolScriptLanguage: QueryLanguageEntry = {
  id: 'molscript',
  name: 'MolScript',
  hint: 'MolScript, e.g. "(sel.atom.res-name HIS)"; or an address list ("A:100-120").',
  transpile,
};
