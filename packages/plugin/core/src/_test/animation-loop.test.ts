/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { PluginAnimationLoop } from '../animation-loop.js';
import { PluginConfig } from '../config.js';
import type { PluginContext } from '../context.js';

const clock = vi.hoisted(() => ({ t: 0 }));
vi.mock('@molstar/core/util/now', () => ({ now: () => clock.t }));

// The loop with `PauseLoopWhenIdle`: it stops asking for frames once the canvas has been idle for a while, and asks
// again on the canvas' `wake` or a state animation starting. Driven by hand: frames are queued, time is faked.

function setup(pauseWhenIdle: boolean) {
  clock.t = 0;
  const frames: FrameRequestCallback[] = [];
  const canvas = {
    isIdle: false,
    wake: new Subject<void>(),
    webgl: { xr: { session: undefined } },
    // As Canvas3D: the last callback asked for is kept, asked for again on `resume`, and cleared by a cancel.
    last: undefined as FrameRequestCallback | undefined,
    resume() {
      if (canvas.last) frames.push(canvas.last);
    },
    pause() {},
    resetTime() {},
    tick: vi.fn(),
    requestAnimationFrame(cb: FrameRequestCallback) {
      canvas.last = cb;
      frames.push(cb);
      return frames.length;
    },
    cancelAnimationFrame() {
      canvas.last = undefined;
    },
  };
  const isAnimating = new BehaviorSubject(false);
  const animation = { isAnimating: false, tick: async () => {} };
  const plugin = {
    canvas3d: canvas,
    config: { get: (item: unknown) => (item === PluginConfig.General.PauseLoopWhenIdle ? pauseWhenIdle : undefined) },
    managers: { animation },
    behaviors: { state: { isAnimating } },
  } as unknown as PluginContext;
  const loop = new PluginAnimationLoop(plugin);
  /** Runs the queued frame, `ms` later; false when none is queued. */
  const step = (ms = 16) => {
    clock.t += ms;
    const f = frames.shift();
    if (!f) return false;
    f(clock.t);
    return true;
  };
  return { loop, canvas, animation, isAnimating, frames, step };
}

describe('PluginAnimationLoop, paused when idle', () => {
  it('sleeps after the canvas has been idle for a while, and wakes on the canvas wake', () => {
    const { loop, canvas, frames, step } = setup(true);
    loop.start();
    canvas.isIdle = true;
    let n = 0;
    while (step() && n < 1000) n++;
    expect(n).toBeGreaterThan(5);
    expect(n).toBeLessThan(40);
    expect(loop.isSleeping).toBe(true);
    expect(loop.isAnimating).toBe(true);
    expect(frames.length).toBe(0);

    canvas.isIdle = false;
    canvas.wake.next();
    expect(loop.isSleeping).toBe(false);
    expect(frames.length).toBe(1);
    step();
    expect(frames.length).toBe(1);
    loop.dispose();
  });

  it('a start after a sleep and a stop runs one frame chain, not two', () => {
    const { loop, canvas, frames, step } = setup(true);
    loop.start();
    canvas.isIdle = true;
    while (step());
    expect(loop.isSleeping).toBe(true);
    // Out of view, then in view again (canvas-presence): Canvas3D's resume must not ask for the old callback too.
    loop.stop();
    canvas.isIdle = false;
    loop.start();
    expect(frames.length).toBe(1);
    step();
    expect(frames.length).toBe(1);
    loop.dispose();
  });

  it('keeps going while a state animation runs, and wakes when one starts', () => {
    const { loop, canvas, animation, isAnimating, frames, step } = setup(true);
    loop.start();
    canvas.isIdle = true;
    animation.isAnimating = true;
    for (let i = 0; i < 100; i++) expect(step()).toBe(true);
    animation.isAnimating = false;
    while (step());
    expect(loop.isSleeping).toBe(true);
    isAnimating.next(true);
    expect(frames.length).toBe(1);
    loop.dispose();
  });

  it('a wake while stopped asks for nothing; off by default it never sleeps', () => {
    const asleep = setup(true);
    asleep.loop.start();
    asleep.loop.stop();
    asleep.frames.length = 0;
    asleep.canvas.wake.next();
    expect(asleep.frames.length).toBe(0);
    asleep.loop.dispose();

    const { loop, canvas, step } = setup(false);
    loop.start();
    canvas.isIdle = true;
    for (let i = 0; i < 200; i++) expect(step()).toBe(true);
    expect(loop.isSleeping).toBe(false);
    loop.dispose();
  });
});
