/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Script as ScriptData } from '@molstar/core/util/script';
import type { Expression } from './language/expression.js';
import { parseMolScript } from './language/parser.js';
import { transpileMolScript } from './script/mol-script/symbols.js';
import { getRegisteredLanguages, parse } from './transpile.js';

export type Script = ScriptData;

export function Script(expression: string, language: Script['language']): Script {
  return { expression, language };
}

/** Text-to-expression operations; molecular compilation and evaluation live in @molstar/model. */
export namespace Script {
  export type Language = ScriptData['language'];
  export const Info: Record<Language, string> = {
    'mol-script': 'Mol-Script',
    pymol: 'PyMOL',
    vmd: 'VMD',
    jmol: 'Jmol',
  };

  export function is(x: any): x is Script {
    return !!x && typeof (x as Script).expression === 'string' && !!(x as Script).language;
  }

  export function areEqual(a: Script, b: Script) {
    return a.language === b.language && a.expression === b.expression;
  }

  /** MolScript plus the languages enabled through explicit transpiler registration imports. */
  export function getAvailableLanguages(): Language[] {
    return ['mol-script', ...getRegisteredLanguages()];
  }

  export function toExpression(script: Script): Expression {
    if (script.language === 'mol-script') {
      const parsed = parseMolScript(script.expression);
      if (parsed.length === 0) throw new Error('No query');
      return transpileMolScript(parsed[0]);
    }
    return parse(script.language, script.expression);
  }
}
