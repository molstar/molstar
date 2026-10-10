/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Script } from '@molstar/core/util/script';
import type { Expression } from './language/expression.js';
import { parseMolScript } from './language/parser.js';
import { SymbolMap as MolQLSymbols } from './language/symbol-table.js';
import { expressionValidationIssues } from './language/validation.js';
import { SymbolMap as MolScriptSymbols, transpileMolScript } from './script/mol-script/symbols.js';
import { transpiler as pymol } from './transpilers/pymol/parser.js';
import { transpiler as vmd } from './transpilers/vmd/parser.js';
import { transpiler as jmol } from './transpilers/jmol/parser.js';

/** Translate supported selection text to MolQL and validate its syntax and argument names.
 * MolScript input must contain exactly one expression; aliases and macros are validated before expansion.
 * Does not compile an executable query or evaluate against molecular data. All languages are available without
 * registration imports. Import individual transpiler parsers for single-language consumers. */
export function compileScript(language: Script['language'], source: string): Expression {
  let expression: Expression;
  switch (language) {
    case 'mol-script':
      expression = compileMolScript(source);
      break;
    case 'pymol':
      expression = pymol(source);
      break;
    case 'vmd':
      expression = vmd(source);
      break;
    case 'jmol':
      expression = jmol(source);
      break;
    default:
      throw new Error(`Unsupported script language '${language}'.`);
  }
  const issues = expressionValidationIssues(expression);
  if (issues) throw new Error(`Invalid MolQL expression:\n${issues.join('\n')}`);
  return expression;
}

function compileMolScript(source: string): Expression {
  const parsed = parseMolScript(source);
  if (parsed.length === 0) throw new Error('No query');
  if (parsed.length !== 1) throw new Error('Expected exactly one MolScript expression.');
  const issues = expressionValidationIssues(parsed[0], {
    getSymbol: (name) => MolScriptSymbols[name]?.symbol ?? MolQLSymbols[name],
  });
  if (issues) throw new Error(`Invalid MolScript expression:\n${issues.join('\n')}`);
  return transpileMolScript(parsed[0]);
}
