/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of validate.heap.ts (`test:perf`): a table cell checked in every column form, and the validation of the
// coded genome table (validate-cases.ts, the guard's cases), each with its positive control. The genome's step is one
// validation of 100,000 rows (the guard compares 750,000 with 3M): bytes a validation, 232 KB at 100 rows and at a
// million alike (2026-10-14), so a byte a row shows. Its control is off: 232 KB a step collects every step, and the
// accounting's noise (about 1 KB a step) hides `leak`'s 88 bytes (it read -759).

import { describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { CELL_COLUMNS, cellStep, genome } from './validate-cases.js';

describe('table validation, quick', () => {
  it.for(CELL_COLUMNS)('a cell of $name', async (c, { task }) => {
    const { step, scan } = cellStep(c.make());
    quickCase(task, `validate: a cell of ${c.name}`, step);
    expect(scan.bad, 'every cell is valid').toBe(0);
  });

  it('the coded genome table, validated', async ({ task }) => {
    const { step } = genome(100_000);
    quickCase(task, 'validate: the coded genome table, 100,000 rows', step, {
      warmSteps: 100,
      steps: 50,
      windows: 3,
      maxWarmSteps: 1_000,
      control: false,
      why: "a validation allocates 232 KB and collects every step: the accounting's noise hides leak's 88 bytes",
    });
  });
});
