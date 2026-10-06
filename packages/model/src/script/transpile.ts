/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 *
 * Adapted from MolQL src/transpile.ts
 */

import type { Transpiler } from './transpilers/transpiler.js';
import type { Expression } from './language/expression.js';
import type { Script } from './script.js';

const transpilers = new Map<Script.Language, Transpiler>();

/** Enable a script language, called by the `transpilers/<lang>` modules on import */
export function registerTranspiler(lang: Script.Language, transpiler: Transpiler) {
  transpilers.set(lang, transpiler);
}

export function getRegisteredLanguages(): Script.Language[] {
  return Array.from(transpilers.keys());
}

export function parse(lang: Script.Language, str: string): Expression {
  const transpiler = transpilers.get(lang);
  if (!transpiler) throw new Error(`Script language '${lang}' is not available in this build`);
  try {
    const query = transpiler(str);
    return query;
  } catch (e) {
    console.error(e.message);
    throw e;
  }
}
