/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { Subject } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import type { PluginContext } from '@molstar/plugin/context';
import { CanvasContexts, followCanvasPresence, unfollowCanvasPresence } from '../canvas-presence.js';

// The canvas' presence (canvas-presence.ts) with fake plugins: a fake WebGL context whose loss and restore fire their
// events a task later (as browsers do), an IntersectionObserver driven by hand.

const observers = new Map<Element, (entries: { isIntersecting: boolean }[]) => void>();
class FakeIntersectionObserver {
  constructor(private cb: (entries: { isIntersecting: boolean }[]) => void) {}
  observe(el: Element) {
    observers.set(el, this.cb);
  }
  disconnect() {}
}

const tick = () => new Promise<void>((r) => setTimeout(r, 0));

function fakeViewer(inView: boolean) {
  const root = document.createElement('div');
  root.getBoundingClientRect = () =>
    (inView
      ? { width: 100, height: 100, top: 0, left: 0, bottom: 100, right: 100 }
      : { width: 100, height: 100, top: 5000, left: 0, bottom: 5100, right: 100 }) as DOMRect;
  let lost = false;
  let lostEventSent = false;
  const contextLost = new Subject<number>();
  const contextRestored = new Subject<number>();
  const ext = {
    loseContext() {
      lost = true;
      setTimeout(() => {
        lostEventSent = true;
        contextLost.next(0);
      });
    },
    restoreContext() {
      if (!lostEventSent) throw new Error('restore not allowed yet');
      setTimeout(() => {
        lost = false;
        lostEventSent = false;
        contextRestored.next(0);
      });
    },
  };
  const loop = { isAnimating: true, starts: 0, stops: 0 };
  const plugin = {
    layout: { root },
    canvas3dContext: { contextLost },
    canvas3d: {
      webgl: {
        get isContextLost() {
          return lost;
        },
        xr: { session: undefined },
        gl: { getExtension: () => ext },
        contextRestored,
      },
      requestDraw() {},
    },
    animationLoop: {
      get isAnimating() {
        return loop.isAnimating;
      },
      start() {
        loop.isAnimating = true;
        loop.starts++;
      },
      stop() {
        loop.isAnimating = false;
        loop.stops++;
      },
    },
  } as unknown as PluginContext;
  return {
    plugin,
    loop,
    root,
    isLost: () => lost,
    show: (on: boolean) => observers.get(root)?.([{ isIntersecting: on }]),
  };
}

describe('canvas presence', () => {
  const made: ReturnType<typeof fakeViewer>[] = [];
  const make = (inView: boolean) => {
    const v = fakeViewer(inView);
    followCanvasPresence(v.plugin);
    made.push(v);
    return v;
  };
  beforeEach(() => {
    (globalThis as any).IntersectionObserver = FakeIntersectionObserver;
    CanvasContexts.budget = 3;
    CanvasContexts.max = 4;
  });
  afterEach(() => {
    for (const v of made.splice(0)) unfollowCanvasPresence(v.plugin);
    observers.clear();
    CanvasContexts.budget = 12;
    CanvasContexts.max = 15;
  });

  it('stops the loop out of view and starts it in view', async () => {
    const v = make(true);
    v.show(false);
    await tick();
    expect(v.loop.isAnimating).toBe(false);
    v.show(true);
    await tick();
    expect(v.loop.isAnimating).toBe(true);
  });

  it('keeps live contexts within the budget, releasing those out of view first', async () => {
    const shown = [make(true), make(true)];
    const hidden = [make(false), make(false), make(false)];
    await tick();
    await tick();
    expect(CanvasContexts.stats().live).toBe(3);
    for (const v of shown) expect(v.isLost()).toBe(false);
    expect(hidden.filter((v) => v.isLost()).length).toBe(2);
  });

  it('restores a context when its viewer comes into view, without passing the budget', async () => {
    const a = make(true);
    const b = make(false);
    const c = make(false);
    const d = make(false);
    await tick();
    await tick();
    const lostOnes = [b, c, d].filter((v) => v.isLost());
    expect(lostOnes.length).toBeGreaterThan(0);
    const back = lostOnes[0];
    a.show(false);
    back.show(true);
    await tick();
    await tick();
    await tick();
    expect(back.isLost()).toBe(false);
    expect(CanvasContexts.stats().live).toBeLessThanOrEqual(CanvasContexts.budget);
    expect(CanvasContexts.stats().dark).toBe(0);
  });

  it('never restores before the browser said the context was lost', async () => {
    const a = make(true);
    const b = make(true);
    const c = make(true);
    const d = make(false);
    await tick();
    await tick();
    // d was released at mount; show it at once, then hide another: the restore waits for d's lost event.
    expect(d.isLost()).toBe(true);
    c.show(false);
    d.show(true);
    await tick();
    await tick();
    await tick();
    expect(d.isLost()).toBe(false);
    expect([a, b].every((v) => !v.isLost())).toBe(true);
  });

  it('lets viewers in view use up to max when more than the budget show', async () => {
    for (let i = 0; i < 5; i++) make(false);
    await tick();
    await tick();
    for (const v of made) v.show(true);
    for (let i = 0; i < 6; i++) await tick();
    const s = CanvasContexts.stats();
    expect(s.live).toBe(CanvasContexts.max);
    expect(s.dark).toBe(1);
  });
});
