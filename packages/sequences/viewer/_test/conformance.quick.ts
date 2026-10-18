/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of conformance.heap.ts (`test:perf`): the same kit runs (conformance-cases.ts), every frame the kit
// hands its `heap` hook a quick case with its positive control (`KIT` counts, ids by kit-cases.ts `frameIds`): the
// stars', linear and on rings, the rows' and `seqv:map`'s, each at every zoom the kit draws.

import { expect, it } from 'vitest';
import { RN3, rn3Kit, starsKit } from './conformance-cases.js';
import { quickCase } from './heap.js';
import { frameIds, KIT } from './kit-cases.js';

it('a stars frame, at every zoom and on rings, quick', async ({ task }) => {
  const id = frameIds('conformance: stars');
  const heap = (await starsKit((frame, at) => quickCase(task, id(at), frame, KIT))).filter((c) => c.check === 'heap');
  expect(heap.filter((c) => c.geometry === 'circular').length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
});

for (const rn3 of RN3)
  it(`${rn3[0].name} at every zoom the kit draws, quick`, async ({ task }) => {
    const id = frameIds(`conformance: ${rn3[0].name}`);
    const heap = (await rn3Kit(rn3, (frame, at) => quickCase(task, id(at), frame, KIT))).filter(
      (c) => c.check === 'heap',
    );
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  });
