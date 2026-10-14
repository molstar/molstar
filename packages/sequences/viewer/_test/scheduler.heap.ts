/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tracked scheduler allocates nothing per frame or timer with a stable callback (seqspec:viewer#performance): an
// animation through `ctx.scheduler`, a step a frame requested and fired, over runs of n and 4n. The port is a fake
// that keeps its one pending frame and timer in fields, so only the tracked layer is measured; its clock is fractional,
// as a browser's is, so a due time boxed as a double shows. The clock lives in a `Float64Array`: a double set on a
// closure's `let` is boxed at each set (16 bytes a step, the fake's own, measured in Node 25). A frame gets the clock's
// whole ms: a fractional double passed to the frame's wrapper, a call V8 does not inline there, was boxed at the fake's
// call (16 bytes a step once the steps settled, so the guard failed every run on the pinned Node, 2026-10-14; a browser's
// rAF hands its own new double each frame). The positive control: the same steps allocating an object each
// (`leaking`), measured alike, are caught.

import { expect, it } from "vitest";
import { createSeqViewer, type Scheduler } from "../src/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";

const N = 20_000;

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

/** A context over the fake port, and its step: a frame and a timer requested, then fired. */
function looping() {
  const { port, tick } = fieldPort();
  const ctx = createSeqViewer(undefined, { scheduler: port });
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

it("a frame loop and a timer with stable callbacks allocate nothing through ctx.scheduler", () => {
  const l = looping();
  const r = heapScaling(N, l.step);
  expect(l.frames).toBe(0);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  l.ctx.dispose();
});

it("the same steps allocating an object each are caught (positive control)", () => {
  const l = looping();
  expect(controlCaught(() => heapScaling(N, leaking(l.step)))).toBe(true);
  l.ctx.dispose();
});
