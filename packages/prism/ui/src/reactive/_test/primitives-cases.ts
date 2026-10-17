/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The base primitives' steps (rendering.md §2.2–2.3), shared by their quick cases (primitives.quick.ts) and heap guards
// (primitives.heap.ts): each a step on a set-up state, and a sink that keeps its work from being optimized away. The
// target is no allocation a step: what a binding or a render writes goes to a plain object here, so what is read is the
// runtime's own.

import { bind, createRoot, flush, frameState, watch } from '../../dom/scope.js';
import { flushFrame } from '../frame-state.js';
import { keyedSelector } from '../keyed-selector.js';
import { keyedState } from '../keyed-state.js';
import { derive, State } from '../state.js';
import { createWorkQueue } from '../work-queue.js';

export interface Case {
  step: (i: number) => void;
  sink: () => number;
}

/** A state set, told to 10 listeners. */
export function stateSet(): Case {
  const s = new State(0);
  let sum = 0;
  for (let k = 0; k < 10; k++) s.listen((v) => (sum = (sum + v) & 0xffff));
  return { step: (i) => s.set(i), sink: () => sum };
}

/** A derived value recomputed and told (its source changes, the value with it). */
export function deriveRecompute(): Case {
  const a = new State(0);
  const d = derive(a, (v) => v * 2);
  let sum = 0;
  d.listen!((v) => (sum = (sum + v) & 0xffff));
  return { step: (i) => a.set(i), sink: () => sum };
}

/** A derived value recomputed to the same value: skipped, its listener not told. */
export function deriveSkip(): Case {
  const a = new State(0);
  const d = derive(a, (v) => v >= 0);
  let told = 0;
  d.listen!(() => (told = (told + 1) & 0xffff));
  return { step: (i) => a.set(i + 1), sink: () => told };
}

/** A binding's write (to a plain object: the runtime's own cost). */
export function bindWrite(): Case {
  const a = new State(0);
  const target = { v: 0 };
  createRoot(() => bind(a, (v) => (target.v = v & 0xffff)));
  return { step: (i) => a.set(i), sink: () => target.v };
}

/** A `watch` marked by a set and rendered by `flush`. */
export function watchRender(): Case {
  const a = new State(0);
  const b = new State(0);
  const target = { v: 0 };
  createRoot(() => watch([a, b], () => (target.v = (a.value + b.value) & 0xffff)));
  return {
    step: (i) => {
      a.set(i);
      flush();
    },
    sink: () => target.v,
  };
}

/** The current one of 1,000 keys moving: two keys told. */
export function keyedSelectorMove(): Case {
  const current = new State(0);
  const sel = keyedSelector(current);
  let told = 0;
  for (let k = 0; k < 1000; k++) sel.key(k).listen(() => (told = (told + 1) & 0xffff));
  return { step: (i) => current.set(i % 1000), sink: () => told };
}

/** One of 1,000 keys' values set: that key told. */
export function keyedStateSet(): Case {
  const ks = keyedState<number, number>();
  let told = 0;
  for (let k = 0; k < 1000; k++) ks.key(k).listen(() => (told = (told + 1) & 0xffff));
  return { step: (i) => ks.set(i % 1000, i), sink: () => told };
}

/** Ten offers to a frame state in a frame, and its tick: one publish, one render. */
export function frameOffers(): Case {
  const target = { v: 0, renders: 0 };
  let offer!: (v: number) => void;
  createRoot(() => {
    const pointer = frameState(0);
    offer = (v) => pointer.offer(v);
    watch([pointer], () => {
      target.v = pointer.value & 0xffff;
      target.renders = (target.renders + 1) & 0xffff;
    });
  });
  return {
    step: (i) => {
      for (let k = 0; k < 10; k++) offer(i * 10 + k);
      flushFrame();
    },
    sink: () => target.renders,
  };
}

/** Ten wants of one key in a frame (the last wins), and the frame's run: one run. The work is made once (a caller's). */
export function workWant(): Case {
  const q = createWorkQueue({ lane: 'frame' });
  const target = { v: 0, runs: 0 };
  let next = 0;
  const run = () => {
    target.v = next & 0xffff;
    target.runs = (target.runs + 1) & 0xffff;
  };
  return {
    step: (i) => {
      for (let k = 0; k < 10; k++) {
        next = i * 10 + k;
        q.want('k', run);
      }
      flushFrame();
    },
    sink: () => target.runs,
  };
}

/**
 * A frame's run of 50 keys, each with a priority asked at the run and a `stillWanted` (one in ten no longer wanted):
 * ordered, asked, run.
 */
export function workFrameRun(): Case {
  const q = createWorkQueue({ lane: 'frame' });
  const target = { v: 0 };
  let at = 0;
  const keys = Array.from({ length: 50 }, (_, k) => ({
    run: () => (target.v = (target.v + k) & 0xffff),
    options: { priority: () => 50 - Math.abs(k - at), stillWanted: () => k % 10 !== 0 },
  }));
  return {
    step: (i) => {
      at = i % 50;
      for (let k = 0; k < 50; k++) q.want(k, keys[k].run, keys[k].options);
      flushFrame();
    },
    sink: () => target.v,
  };
}
