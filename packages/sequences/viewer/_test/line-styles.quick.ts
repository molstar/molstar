/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of line-styles.heap.ts (`test:perf`): the same cases (line-styles-cases.ts), each with its positive
// control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { lineStyles, STYLED } from './line-styles-cases.js';

describe('line styles, quick', () => {
  it('a frame of every style in tracks, panned', async ({ task }) => {
    const d = await lineStyles('tracks');
    quickCase(task, 'line-styles: tracks', d.step, STYLED);
    expect(d.state().paths).toBeGreaterThan(5);
  });

  it('a frame of every style in the text, scrolled', async ({ task }) => {
    const d = await lineStyles('text');
    quickCase(task, 'line-styles: text', d.step, STYLED);
    expect(d.state().paths).toBeGreaterThan(10);
  });
});
