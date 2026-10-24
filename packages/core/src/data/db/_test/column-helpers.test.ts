/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '../column.js';
import { getArrayBounds } from '../column-helpers.js';

describe('column helpers', () => {
  it('getArrayBounds clamps start to the row count (perf-audit-2 01/B8)', () => {
    expect(getArrayBounds(4, { start: 3, end: 3 })).toEqual({ start: 3, end: 3 });
    expect(getArrayBounds(4, { start: 4 })).toEqual({ start: 4, end: 4 });
    expect(getArrayBounds(4, { start: 6, end: 8 })).toEqual({ start: 4, end: 4 });
    expect(getArrayBounds(4, { start: 2, end: 1 })).toEqual({ start: 2, end: 2 });
    expect(getArrayBounds(4, { start: -1, end: 2 })).toEqual({ start: 0, end: 2 });
  });

  it('toArray from the end has no rows', () => {
    expect(Column.ofArray({ array: [1, 2, 3, 4], schema: Column.Schema.int }).toArray({ start: 4, end: 4 })).toEqual(
      [],
    );
    expect(Column.ofIntArray([1, 2, 3, 4]).toArray({ start: 4 }).length).toBe(0);
  });
});
