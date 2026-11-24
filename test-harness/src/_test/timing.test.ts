/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect } from 'vitest';
import { itWithin } from '../timing.js';

// The heavy lane's bounds are the best of `BoundedTries` runs: a first run slowed by a burst runs again, setup and all.

describe('timing: itWithin', () => {
  let runs = 0;
  itWithin('a bound missed once and then held passes', (within) => {
    runs++;
    const first = runs === 1;
    within(25, () => {
      // the first run meets a burst (a busy wait longer than the bound); the next does not
      const end = performance.now() + (first ? 60 : 0);
      while (performance.now() < end);
    });
    expect(runs).toBeGreaterThan(0);
  });
});
