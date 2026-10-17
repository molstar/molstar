/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of primitives.heap.ts (`test:perf`): the base primitives' bytes a step, each with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../../_test/perf.js';
import {
  bindWrite,
  deriveRecompute,
  deriveSkip,
  frameOffers,
  keyedSelectorMove,
  keyedStateSet,
  stateSet,
  watchRender,
  workFrameRun,
  workWant,
} from './primitives-cases.js';

const cases = {
  stateSet,
  deriveRecompute,
  deriveSkip,
  bindWrite,
  watchRender,
  keyedSelectorMove,
  keyedStateSet,
  frameOffers,
  workWant,
  workFrameRun,
};

for (const [name, make] of Object.entries(cases)) {
  it(name, ({ task }) => {
    const { step, sink } = make();
    quickCase(task, `primitives: ${name}`, step);
    expect(sink()).not.toBe(0.5);
  });
}
