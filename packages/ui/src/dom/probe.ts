/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewScope } from './scope.js';

/**
 * For tests and benches (dom §9): unset, nothing is wrapped or told (one branch a binding write and a render). Set
 * before rendering:
 * - `on` is told every binding write and `watch` render with its scope (a test counts those in paused scopes);
 * - `run` runs each write and render (`fn`) with where it was made (`site`, a stack) and whether it is a binding's
 *   first write since it started (`first`: made, or resumed after a pause; always false for a render), so a test can
 *   see whether it wrote to the DOM.
 */
export const DomProbe: {
  on?: (kind: 'write' | 'render', scope: ViewScope) => void;
  run?: (kind: 'write' | 'render', site: string, first: boolean, fn: () => void) => void;
} = {};

/** Counts for benches: `watch` renders and binding writes (a plain increment each). */
export const DomCounts = { render: 0, write: 0 };
