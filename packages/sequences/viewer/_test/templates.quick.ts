/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of templates.heap.ts (`test:perf`): a pan's pass and composer frame over the template's 1,000 instance
// tracks (templates-cases.ts, the guard's case), with its positive control. The guard's other case, an instance that
// leaves not being retained (a WeakRef after a collection), is not a step and has no twin.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { INSTANCES, templated } from './templates-cases.js';

it("a pan's pass and frame over the template's instance tracks", async ({ task }) => {
  const g = await templated();
  expect(g.t.ctx.graph.assembled.value?.view?.resolved.tracks.count).toBeGreaterThan(INSTANCES);
  quickCase(task, 'templates: pan pass and frame', g.frame);
});
