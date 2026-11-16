/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick allocation tier (`test:perf`) of the runtime's sinks (dom §2): a change marked into a binding or a watch,
// and its write or render, allocate nothing; ten frame offers and their tick neither. Each with its positive control.

import { batch, State } from '@molstar/core/reactive';
import { quickCase } from '@molstar/test-harness/quick';
import { expect, it } from 'vitest';
import { FrameState, flushFrame } from '../frame.js';
import { createRoot } from '../scope.js';
import { bind, flush, watch } from '../watch.js';

/** What the steps write to, so nothing they compute is optimized away; kept a small integer (no number boxed). */
let out = 0;
const add = (v: number) => {
  out = (out + v) & 0xffff;
};

const cases: Record<string, () => (i: number) => void> = {
  /** A state with one binding: a set, the binding's write. */
  bindWrite() {
    const s = new State(0);
    createRoot(() => bind(s, add));
    return (i) => s.set((i & 0xffff) + 1);
  },
  /** A watch of five states: one set marks it, a flush renders it. */
  watchOfFive() {
    const states = Array.from({ length: 5 }, () => new State(0));
    createRoot(() =>
      watch(states, () => {
        for (const s of states) add(s.value);
      }),
    );
    return (i) => {
      states[i % 5].set((i & 0xffff) + 1);
      flush();
    };
  },
  /** A set and a set back in a batch: marked, nothing written or rendered. */
  setBack() {
    const s = new State(0);
    createRoot(() => {
      bind(s, add);
      watch([s], () => add(s.value));
    });
    // One function for every step (a closure a step would be the test's allocation, not the runtime's).
    let next = 0;
    const setAndBack = () => {
      s.set(next);
      s.set(0);
    };
    return (i) => {
      next = (i & 0xffff) + 1;
      batch(setAndBack);
      flush();
    };
  },
  /** Ten offers to a frame state watched, then the frame's tick. */
  frameOffers() {
    const f = new FrameState(0);
    createRoot(() => watch([f], () => add(f.value)));
    return (i) => {
      for (let k = 0; k < 10; k++) f.offer(((i & 0xfff) << 4) + k);
      flushFrame();
    };
  },
};

for (const [name, make] of Object.entries(cases)) {
  it(name, ({ task }) => {
    quickCase(task, `sinks: ${name}`, make());
    expect(out).not.toBe(0.5);
  });
}
