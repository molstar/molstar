/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of themes.heap.ts (`test:perf`): the same cases (themes-cases.ts) in one measurement, with its positive
// control. A step reads one element through the next theme of the table; the per-theme detail is the deep guard's.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { CASES, reading } from './themes-cases.js';

/** A step is one element, well under a microsecond: windows of 5,000. */
const ELEMENT: QuickOptions = { warmSteps: 5_000, steps: 5_000, windows: 5, minWarmSteps: 10_000 };

describe('theme values per element, quick', () => {
  it('every built-in theme and rule table, an element a step', ({ task }) => {
    // One measurement: step `i` reads an element through theme `i % K`. The per-theme detail is the deep guard's.
    const readers = CASES.map(({ theme }) => {
      expect(theme).not.toBeNull();
      return reading(theme);
    });
    quickCase(task, 'themes: every theme', (i) => readers[i % readers.length](i), ELEMENT);
  });
});
