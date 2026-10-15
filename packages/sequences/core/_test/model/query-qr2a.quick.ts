/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of query-qr2a.heap.ts (`test:perf`): `onto`, `partners`, `expand`, `group`, `pick`, `first`, `ids` and an
// item test (query-qr2a-cases.ts, the guard's cases), each a run of a compiled query over a million positions, measured
// in bytes a run, with its positive control.
// A run allocates its result's wrappers (the typed arrays' headers, the selection object): a constant a run.

import { describe, expect, it } from 'vitest';
import { prepare, QUERIES } from './query-qr2a-cases.js';
import { quickCase } from './quick.js';

/**
 * A run of several milliseconds, in bytes: 3 to 40 KB a run (the result's wrappers, the selection). Its control is off:
 * a step that collects every few steps loses or gains hundreds of bytes a step to the accounting, `leak`'s 88 bytes
 * unreadable in it (half the cases read the control at -700 to +30, on either run of 2026-10-14), so `leaking` cannot
 * be read on a query run; the deep guard's `(leak Q)` control is the one that catches a query's allocation.
 */
const RUN = {
  warmSteps: 20,
  steps: 5,
  windows: 3,
  maxWarmSteps: 40,
  control: false as const,
  why: "a query run collects every few steps: the accounting's noise hides leak's 88 bytes",
};

describe("QR2a's runtimes, quick", () => {
  it.for(QUERIES.map(([text, groups]) => ({ text, groups })))('$text', async ({ text, groups }, { task }) => {
    const { run, groups: selected } = prepare(text);
    quickCase(task, `query-qr2a: ${text}`, run, RUN);
    expect(selected()).toBeGreaterThanOrEqual(groups);
  });
});
