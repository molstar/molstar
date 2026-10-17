/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { insertedRun } from '../window.js';

// A list's reveal animates only one change: a run of rows inserted (a node opened), nothing else.

const id = (x: string) => x;

describe('insertedRun', () => {
  it('finds a run inserted at the start, in the middle and at the end', () => {
    expect(insertedRun(['a', 'b'], ['x', 'y', 'a', 'b'], id)).toEqual({ at: 0, count: 2 });
    expect(insertedRun(['a', 'b', 'c'], ['a', 'x', 'y', 'b', 'c'], id)).toEqual({ at: 1, count: 2 });
    expect(insertedRun(['a', 'b'], ['a', 'b', 'x'], id)).toEqual({ at: 2, count: 1 });
  });

  it('is undefined for any other change, and for the first rows', () => {
    expect(insertedRun([], ['a'], id)).toBeUndefined();
    expect(insertedRun(['a', 'b'], ['a'], id)).toBeUndefined();
    expect(insertedRun(['a', 'b'], ['b', 'a', 'x'], id)).toBeUndefined();
    expect(insertedRun(['a', 'b', 'c'], ['x', 'a', 'b', 'y', 'c'], id)).toBeUndefined();
    expect(insertedRun(['a', 'b'], ['a', 'b'], id)).toBeUndefined();
  });
});
