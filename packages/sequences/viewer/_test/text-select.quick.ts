/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of text-select.heap.ts (`test:perf`): the select drag across four lines, and the scroll past the rows'
// edges, a frame a step (text-select-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { textSelect } from './text-select-cases.js';

describe('a select drag in text, quick', () => {
  it('moves across four lines, past their ends and over their labels bands', async ({ task }) => {
    const d = await textSelect('lines');
    quickCase(task, 'text-select: lines', d.step);
    const { owner, ends, active } = d.state();
    expect([owner, active]).toEqual(['select', true]);
    expect(ends).toBeGreaterThan(1_000);
  });

  it("the scroll past the rows' edges, a frame a step", async ({ task }) => {
    const d = await textSelect('scroll');
    quickCase(task, 'text-select: scroll', d.step);
    const { owner, ends, scrolls } = d.state();
    expect(owner).toBe('select');
    expect(scrolls).toBeGreaterThan(1_000);
    expect(ends).toBeGreaterThan(1_000);
  });
});
