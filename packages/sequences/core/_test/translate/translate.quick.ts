/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of translate.heap.ts (`test:perf`): a pointer move through the cached routes and a batch of runs
// translated (translate-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { resetVisits, TRANSLATE_CASES, translatedRuns, visited } from './translate-cases.js';

describe('translate hot paths, quick', () => {
  it.for(TRANSLATE_CASES)('$name', async ({ name, step }, { task }) => {
    resetVisits();
    quickCase(task, `translate: ${name}`, step);
    expect(visited(), 'the steps visit images').toBeGreaterThan(1_000);
  });

  it('translated what it should', () => {
    expect(translatedRuns()).toBe(3);
  });
});
