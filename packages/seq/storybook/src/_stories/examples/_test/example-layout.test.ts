/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The example story's layout keeps the canvas inside the page at any width: every grid has an explicit shrinkable
// column and its items `min-width: 0` (no layout engine here: the declarations are what is checked; the look is in the
// Globins rows stories at a wide and a narrow page).

import { describe, expect, it } from 'vitest';
import { Box, FoundBlock, PageGrid, pairGrid } from '../example-layout.js';

describe('the example story grids', () => {
  for (const [name, css] of [
    ['page', PageGrid],
    ['section', Box],
    ['pair', pairGrid(false)],
    ['stacked pair', pairGrid(true)],
  ] as const) {
    it(`${name}: no implicit auto column, shrinkable`, () => {
      expect(css).toMatch(/grid-template-columns:[^;]*minmax\([^;]*1fr\)/);
      expect(css).toMatch(/min-width: 0/);
    });
  }

  it('the pair breaks to one column below 420 px, never wider than its box', () => {
    expect(pairGrid(false)).toContain('minmax(min(420px, 100%), 1fr)');
  });

  it('a find box’s result wraps its long lines (GTF’s first line was cut at the right)', () => {
    expect(FoundBlock).toContain('white-space: pre-wrap');
    expect(FoundBlock).toContain('overflow-wrap: anywhere');
  });
});
