/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tracked scheduler allocates nothing per frame or timer with a stable callback (seqspec:viewer#performance): an
// animation through `ctx.scheduler`, a step a frame requested and fired, over runs of n and 4n (the case in
// scheduler-cases.ts, shared with the quick twin, scheduler.quick.ts). The positive control: the same steps allocating
// an object each (`leaking`), measured alike, are caught.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { looping } from './scheduler-cases.js';

const N = 20_000;

it('a frame loop and a timer with stable callbacks allocate nothing through ctx.scheduler', () => {
  const l = looping();
  const r = heapScaling(N, l.step);
  expect(l.frames).toBe(0);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  l.ctx.dispose();
});

it('the same steps allocating an object each are caught (positive control)', () => {
  const l = looping();
  expect(controlCaught(() => heapScaling(N, leaking(l.step)))).toBe(true);
  l.ctx.dispose();
});
