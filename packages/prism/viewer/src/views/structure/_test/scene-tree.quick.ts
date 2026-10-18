// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of scene-tree.heap.ts (`test:perf`): the scene tree's hot paths, bytes a step against the baseline.

import { expect, it } from 'vitest';
import { quickCase } from '../../../_test/perf.js';
import { jump, pick, scrollStep } from './scene-tree-cases.js';

it('a pick moving between two rows', ({ task }) => {
  const { step, sink } = pick();
  quickCase(task, 'scene tree: a pick', step);
  expect(sink()).toBeGreaterThan(10);
});

it('a scroll step', ({ task }) => {
  const { step, sink } = scrollStep();
  quickCase(task, 'scene tree: a scroll step', step);
  expect(sink()).toBeGreaterThan(10);
});

it('a jump', ({ task }) => {
  const { step, sink } = jump();
  quickCase(task, 'scene tree: a jump', step, { steps: 200, warmSteps: 400 });
  expect(sink()).toBeGreaterThan(10);
});
