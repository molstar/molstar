/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for validation (VW2): a table cell is checked with nothing allocated, in every column form, on RB0b's
// stepped `heapScaling` (a step is one cell; untimed, `heapScaling` warms by count); and validating the 3M-row coded
// genome table grows the heap by the same small amount as validating a quarter of it (the old check kept about 52 MB of
// objects, a path and a ref per cell). The cases are validate-cases.ts's, shared with the quick twin
// (validate.quick.ts).

import { describe, expect, it } from 'vitest';
import {
  controlCaught,
  describeScaling,
  heapGrowth,
  heapScaling,
  isFlat,
  leak,
  leaking,
  SCALED_LIMIT,
} from '../model/heap.js';
import { CELL_COLUMNS, cellStep, genome } from './validate-cases.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 250_000;

describe('table cells allocate nothing', () => {
  it.each(CELL_COLUMNS.map((c) => [c.name, c] as const))('%s', (_, c) => {
    const { step, scan } = cellStep(c.make());
    const r = heapScaling(N, step);
    expect(scan.bad, 'every cell is valid').toBe(0);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the control: a cell's check with a leak is caught (the template, CF-05)", () => {
    const { step } = cellStep(CELL_COLUMNS[0].make());
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});

describe('validating the coded genome table', () => {
  it('3M rows grow the heap as 750k do: no allocation per row', () => {
    /** `perRow` after each validation, once a row: the control's leak, else nothing. */
    const grow = (n: number, perRow: (k: number) => void = () => {}) => {
      const { step } = genome(n);
      // The least of three: the first runs still compile (a few hundred KB, whatever the size).
      let least = Number.POSITIVE_INFINITY;
      for (let k = 0; k < 3; k++)
        least = Math.min(
          least,
          heapGrowth(() => {
            step();
            for (let r = 0; r < n; r++) perRow(r);
          }),
        );
      return least;
    };
    const small = grow(750_000);
    const large = grow(3_000_000);
    expect(large - small, `${(small / 1024).toFixed(1)} KB, then ${(large / 1024).toFixed(1)} KB`).toBeLessThan(
      SCALED_LIMIT,
    );
    // The template's control (CF-05): an object a row, the old check's fault, is caught by the same comparison.
    expect(
      controlCaught(() => {
        const [short, long] = [grow(750_000, leak), grow(3_000_000, leak)];
        return { short, long, scaled: long - short };
      }),
    ).toBe(true);
  });
});
