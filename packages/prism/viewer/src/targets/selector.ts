/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { KEYWORD_EXPRESSIONS } from './keywords.js';
import type { ComponentExpression, Keyword, MolQLSelector, SelectorT } from './types.js';

// Ported from Stories Untold's `any-molql` (selector.ts). MVS's runtime turns rows into an expression through
// `StructureElement.Schema.toExpression` (the generator `element-schema`, which resolves rows through the model's
// indices); that is called here directly, since the MVS runtime module that wraps it imports the plugin.

export function isMolQLSelector(s: unknown): s is MolQLSelector {
  return !!s && typeof s === 'object' && !Array.isArray(s) && Object.hasOwn(s, 'molql');
}

/** Component expressions as one MolScript expression (what any of them selects). */
export function itemsExpression(items: readonly ComponentExpression[]): Expression {
  return StructureElement.Schema.toExpression({ items: items as StructureElement.SchemaItem[] });
}

/** A selector as one MolScript expression (as Mol* does when it loads an MVS selector). Throws on a malformed one. */
export function selectorExpression(selector: SelectorT): Expression {
  if (typeof selector === 'string') {
    const keyword = KEYWORD_EXPRESSIONS[selector as Keyword];
    if (!keyword) throw new Error(`"${selector}" is not a keyword.`);
    return keyword;
  }
  if (isMolQLSelector(selector)) return selector.molql as Expression;
  if (Array.isArray(selector)) return itemsExpression(selector as readonly ComponentExpression[]);
  return itemsExpression([selector as ComponentExpression]);
}

/**
 * Several selectors as one: a single one stays as it is (MVS takes every form); several become one `{ molql }` that
 * selects what any of them selects.
 */
export function joinSelectors(selectors: readonly SelectorT[]): SelectorT | undefined {
  if (selectors.length === 0) return undefined;
  if (selectors.length === 1) return selectors[0];
  // Component expressions alone stay a list (still what MVS calls items).
  if (selectors.every((s) => typeof s === 'object' && !isMolQLSelector(s))) {
    return selectors.flatMap((s) => (Array.isArray(s) ? s : [s as ComponentExpression]));
  }
  return { molql: MS.struct.combinator.merge(selectors.map(selectorExpression)) };
}
