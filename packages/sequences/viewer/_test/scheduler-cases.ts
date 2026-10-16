/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tracked scheduler's case, for its heap guard and its quick twin (scheduler.heap.ts, scheduler.quick.ts): a
// context over a fake port that keeps its one pending frame and timer in fields, so only the tracked layer is measured,
// and its step, a frame and a timer with stable callbacks requested through `ctx.scheduler`, then fired. The port's
// clock is fractional, as a browser's is, so a due time boxed as a double shows. The clock lives in a `Float64Array`: a
// double set on a closure's `let` is boxed at each set (16 bytes a step, the fake's own, measured in Node 25). A frame
// gets the clock's whole ms: a fractional double passed to the frame's wrapper, a call V8 does not inline there, was
// boxed at the fake's call (16 bytes a step once the steps settled, so the guard failed every run on the pinned Node,
// 2026-10-14; a browser's rAF hands its own new double each frame).

import { createSeqViewer, type Scheduler } from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/spec.js';

/** One frame and one timer at a time, fired by `tick`. */
function fieldPort() {
  let frame: ((t: number) => void) | null = null;
  let timer: (() => void) | null = null;
  let id = 0;
  const clock = new Float64Array(1);
  const port: Scheduler = {
    now: () => clock[0],
    requestFrame(cb) {
      frame = cb;
      return ++id;
    },
    cancelFrame: () => {
      frame = null;
    },
    setTimer(cb) {
      timer = cb;
      return ++id;
    },
    clearTimer: () => {
      timer = null;
    },
  };
  const tick = () => {
    clock[0] += 16.6667;
    const f = frame;
    const t = timer;
    frame = null;
    timer = null;
    f?.(clock[0] | 0);
    t?.();
  };
  return { port, tick };
}

/** A context over the fake port, and its step: a frame and a timer requested, then fired (`frames` stays 0). */
export function looping() {
  const { port, tick } = fieldPort();
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), { scheduler: port });
  const out = { ctx, frames: 0, step: () => {} };
  const onFrame = () => {
    out.frames++;
  };
  const onTimer = () => {
    out.frames--;
  };
  out.step = () => {
    ctx.scheduler.requestFrame(onFrame);
    ctx.scheduler.setTimer(onTimer, 16);
    tick();
  };
  return out;
}
