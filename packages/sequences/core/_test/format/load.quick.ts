/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of load.heap.ts (`test:perf`): the load paths' cases (load-cases.ts, the guard's), each with its
// positive control. A tile-addressing step is a handful of `Tiles` calls; a read is measured in bytes a read: the
// columns it returns (garbage at the next read) and nothing more. A delimited file of 160,000 cells is a step of many
// milliseconds; a Stockholm file is read at 3 rows instead of the guard's 2,000 (TABLE).

import { describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { delimited, stockholm, tileStep } from './load-cases.js';

/** A step of many milliseconds: few steps, in small windows. */
const READ = { warmSteps: 20, steps: 5, windows: 3, maxWarmSteps: 200 };

/**
 * A table over 3 rows (900 positions, 20 to 27 KB a read; the guard's file is 2,000 rows): small enough to be quick, with
 * the counts of a fast step. Its control is off: a read that collects every few steps loses or gains more than
 * `leak`'s 88 bytes to the accounting (a Stockholm read of 3, 4 or 10 rows reads the control at +9, -179 or -19 bytes,
 * of 2,000 rows at -25,806, 2026-10-14), so `leaking` cannot be read on it. The tile and delimited cases keep theirs.
 */
const TABLE = {
  warmSteps: 1_000,
  steps: 250,
  windows: 3,
  control: false as const,
  why: "a read that collects every few steps hides leak's 88 bytes in the accounting's noise",
};

describe('the load paths, quick', () => {
  it('tile addressing', async ({ task }) => {
    const { step, sum } = tileStep();
    quickCase(task, 'load: tile addressing', step);
    expect(sum()).toBeGreaterThan(0);
  });

  it.for(['values', 'ranges'] as const)('a Stockholm tag table as %s', async (geometry, { task }) => {
    const read = stockholm(3).table(geometry);
    let kept = 0;
    quickCase(
      task,
      `load: a Stockholm tag table as ${geometry}`,
      () => {
        kept = read() ? 1 : 0;
      },
      TABLE,
    );
    expect(kept).toBe(1);
  });

  it('a delimited file', async ({ task }) => {
    const { read } = delimited();
    let kept = 0;
    quickCase(
      task,
      'load: a delimited file',
      () => {
        kept = read() ? 1 : 0;
      },
      READ,
    );
    expect(kept).toBe(1);
  });
});
