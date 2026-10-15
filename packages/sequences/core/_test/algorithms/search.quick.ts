/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of search.heap.ts (`test:perf`): one pass of each matcher over the text (search-cases.ts, the guard's
// cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { SEARCH_CASES, searchStep } from './search-cases.js';

describe('search matchers, quick', () => {
  it.for(SEARCH_CASES)('$name', async (c, { task }) => {
    const { step, hits, calls } = searchStep(c);
    quickCase(task, `search: ${c.name}`, step);
    expect(hits(), 'the case visits hits').toBeGreaterThan(calls() * 100);
  });
});
