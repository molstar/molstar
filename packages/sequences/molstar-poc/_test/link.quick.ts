/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of link.heap.ts (`test:perf`): a Mol* hover to member runs and Seq* runs, the guard's step
// (link-cases.ts), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../../core/_test/model/quick.js';
import { molHover } from './link-cases.js';

it('a Mol* hover to member runs and Seq* runs', async ({ task }) => {
  const { step, sink } = await molHover();
  quickCase(task, 'link: a Mol* hover to member runs and Seq* runs', step);
  expect(sink()).not.toBe(0.5);
});
