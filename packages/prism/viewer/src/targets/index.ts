/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Targets: how the viewer names what it acts on (spec/targets.md). Permissive on the way in, MolQL underneath.
// DOM-free; imports the model and the query language, never the plugin. Ported from Stories Untold's `any-molql`.

export { compactNumbers, compactResidues, parseAddressList, type AddressResult } from './addresses.js';
export {
  type Compiled,
  type MolqlFromJson,
  checkMolQL,
  compileQuery,
  formatMolQL,
  fromMolqlJson,
  toMolqlJson,
} from './compile.js';
export {
  type Evaluated,
  type EvaluatedTarget,
  type SummarizeOptions,
  evaluateExpression,
  evaluateTarget,
  queryLoci,
  selectLoci,
  summarize,
} from './evaluate.js';
export { KEYWORD_EXPRESSIONS, KEYWORDS, PERSON_KEYWORDS, readKeyword } from './keywords.js';
export { isMolQLSelector, itemsExpression, joinSelectors, selectorExpression } from './selector.js';
export { compileTarget, formatTarget, parseTarget } from './target.js';
export type * from './types.js';
