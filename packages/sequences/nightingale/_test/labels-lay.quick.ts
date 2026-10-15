/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of labels-lay.heap.ts (`test:perf`): `layLabels` at an unchanged view, the guard's step
// (labels-heap-fixture.ts `laying`), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../../core/_test/model/quick.js';
import { laying } from './labels-heap-fixture.js';

it('layLabels at an unchanged view: laid, compared, nothing written', ({ task }) => {
  const { step, writes } = laying();
  quickCase(task, 'labels-lay: layLabels at an unchanged view', step);
  expect(writes()).toBe(0);
});
