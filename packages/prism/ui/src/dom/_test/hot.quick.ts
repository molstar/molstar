// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hot.heap.ts (`test:perf`): prism-ui's hot paths in jsdom, bytes a step against their baseline.

import { expect, it } from 'vitest';
import { quickCase } from '../../_test/perf.js';
import { commandMenuTyping, forUnchanged } from './hot-cases.js';

it('a For given the same rows in a new array', ({ task }) => {
  const { step, sink } = forUnchanged();
  quickCase(task, 'hot: For, the same rows anew', step);
  expect(sink()).toBe(200);
});

it('typing in the command menu', ({ task }) => {
  const { step, sink } = commandMenuTyping();
  quickCase(task, 'hot: typing in the command menu', step, { steps: 200, warmSteps: 400 });
  expect(sink()).toBeGreaterThan(0);
});
