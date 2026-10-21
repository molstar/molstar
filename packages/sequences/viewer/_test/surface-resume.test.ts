/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

// A paused surface resumed with an animation running (seqspec:viewer#the-render-backend), in the quick lane: a surface
// mounted in jsdom over a 2D context that draws nothing, its composer animating (a skeleton's shimmer, a region's fade).
// Paused, its polling frame is cancelled; resumed, by the host (`setActive`) or the viewer (`ui.active`), a frame draws
// the overlay and the polling runs on. Chromium's own: `canvas/active.test.ts`.

import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { SectionSurface } from '../src/canvas/index.js';
import { createSeqViewer, DefaultSeqViewerSpec, type Scheduler, type SeqViewerContext } from '../src/index.js';
import { fixtureMeasurer } from '../src/testing/index.js';
import { scaleDoc } from './prepare-fixtures.js';

/** A 2D context whose every call does nothing (jsdom has none): what a frame draws is not this test's. */
const noop = () => undefined;
const paint = (canvas: HTMLCanvasElement): unknown =>
  new Proxy(
    { canvas },
    {
      get: (t, k) =>
        k in t
          ? t[k as keyof typeof t]
          : k === 'measureText'
            ? () => ({ width: 0, actualBoundingBoxAscent: 0, actualBoundingBoxDescent: 0 })
            : k === 'getImageData' || k === 'createImageData'
              ? () => ({ data: new Uint8ClampedArray(4), width: 1, height: 1 })
              : noop,
      set: (t, k, v) => {
        (t as Record<string | symbol, unknown>)[k] = v;
        return true;
      },
    },
  );

/** Frames run by `flush`, in order, a second apart: an animation's polling frame is due at each. */
class Frames implements Scheduler {
  private cbs = new Map<number, (t: number) => void>();
  private next = 1;
  private t = 0;
  now() {
    return this.t;
  }
  requestFrame(cb: (t: number) => void) {
    this.cbs.set(this.next, cb);
    return this.next++;
  }
  cancelFrame(id: number) {
    this.cbs.delete(id);
  }
  setTimer() {
    return 0;
  }
  clearTimer() {}
  get pending() {
    return this.cbs.size;
  }
  flush() {
    const run = [...this.cbs.values()];
    this.cbs.clear();
    this.t += 1000;
    for (const cb of run) cb(this.t);
  }
}

beforeAll(() => {
  HTMLCanvasElement.prototype.getContext = function (this: HTMLCanvasElement) {
    return paint(this);
  } as never;
  globalThis.OffscreenCanvas ??= class {
    width: number;
    height: number;
    constructor(width: number, height: number) {
      this.width = width;
      this.height = height;
    }
    getContext() {
      return paint(this as never);
    }
  } as never;
  // jsdom lays nothing out: the surface's size, as the browser would report it.
  Object.defineProperty(HTMLDivElement.prototype, 'clientWidth', { configurable: true, get: () => 900 });
  Object.defineProperty(HTMLDivElement.prototype, 'clientHeight', { configurable: true, get: () => 420 });
  if (!document.fonts)
    Object.defineProperty(document, 'fonts', {
      configurable: true,
      value: { addEventListener: noop, removeEventListener: noop },
    });
  window.matchMedia ??= (() => ({ matches: false, addEventListener: noop, removeEventListener: noop })) as never;
  globalThis.ResizeObserver ??= class {
    observe() {}
    unobserve() {}
    disconnect() {}
  } as never;
});

let done: (() => void) | null = null;
afterEach(() => {
  done?.();
  done = null;
});

async function mounted() {
  const scheduler = new Frames();
  const ctx: SeqViewerContext = createSeqViewer(DefaultSeqViewerSpec(), { scheduler, text: fixtureMeasurer() });
  const r = await ctx.setDocument(scaleDoc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  const host = document.createElement('div');
  document.body.appendChild(host);
  const section = ctx.viewports.sections.value[0];
  const surface = new SectionSurface(ctx, section, host, { dpr: 1 });
  // The composer animating after each of its frames (a skeleton shimmering), as the stream's progress sets it.
  const composer = surface.composer;
  const draw = composer.frame.bind(composer);
  let overlays = 0;
  composer.frame = ((be: unknown, flags: number) => {
    draw(be as never, flags);
    if (flags & 2) overlays++;
    composer.animating = true;
  }) as never;
  surface.renderNow();
  done = () => {
    surface.dispose();
    ctx.dispose();
    host.remove();
  };
  return { ctx, surface, scheduler, overlays: () => overlays };
}

describe('a paused surface resumed with an animation running', () => {
  it('the host resumes it (setActive): its overlay drawn and polled again', async () => {
    const { surface, scheduler, overlays } = await mounted();
    scheduler.flush();
    expect(surface.stats.frames).toBeGreaterThan(0);
    expect(scheduler.pending).toBe(1); // the animation's polling frame
    surface.setActive(false);
    expect(scheduler.pending).toBe(0);
    surface.setActive(true);
    expect(scheduler.pending).toBe(1);
    const before = overlays();
    scheduler.flush();
    expect(overlays()).toBe(before + 1);
    // Running on: the next polling frame asked for, and drawn.
    expect(scheduler.pending).toBe(1);
    scheduler.flush();
    expect(overlays()).toBe(before + 2);
  });

  it('the viewer resumes it (ui.active): the same', async () => {
    const { ctx, scheduler, overlays } = await mounted();
    ctx.setActive(false);
    // The polling frame finds the viewer paused and lapses.
    scheduler.flush();
    expect(scheduler.pending).toBe(0);
    ctx.setActive(true);
    expect(scheduler.pending).toBeGreaterThan(0);
    const before = overlays();
    scheduler.flush();
    expect(overlays()).toBe(before + 1);
    expect(scheduler.pending).toBe(1);
  });
});
