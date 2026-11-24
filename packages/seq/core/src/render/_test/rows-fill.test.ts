/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { rowsFillHeight } from '../rows-fill.js';
import { RowHeight } from '../rows-model.js';

describe('a filling rows band’s height', () => {
  it('takes the room left, at least its default height', () => {
    expect(rowsFillHeight(240, 500.7, undefined, 0)).toBe(500);
    expect(rowsFillHeight(240, 32, undefined, 0)).toBe(240);
  });

  it('is never taller than its rows with their title row: no empty well under them', () => {
    expect(rowsFillHeight(240, 2000, 60, 16)).toBe(60 * RowHeight + 16);
    // four rows: their own height, under the default
    expect(rowsFillHeight(240, 2000, 4, 0)).toBe(4 * RowHeight);
    expect(rowsFillHeight(240, 300, 60, 0)).toBe(300);
  });
});
