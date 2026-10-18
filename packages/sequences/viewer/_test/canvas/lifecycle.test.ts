/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A surface's lifecycle (CS2; seqspec:viewer#the-render-backend), in Chromium: a context that completes takes its
// mounted surfaces with it (a host that called `mountSection` directly leaks nothing), a surface cannot mount on a
// disposed context, a host of no size draws nothing until it has one, and a renderer that throws is omitted and
// reported once while the rest of the frame is drawn.

import { afterEach, describe, expect, it } from 'vitest';
import { mountSection } from '../../src/canvas/index.js';
import { KernelRenderers, type LayerRenderer } from '../../src/index.js';
import { syntheticFamily } from '../families.js';
import { allTypes } from '../prepare-fixtures.js';
import { canvasesOf, close, type Opened, opened, pixelsOf, settle } from './open.js';

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
});

describe('disposal', () => {
  it('ctx.dispose() disposes every mounted surface of the context: canvases released, element and listeners gone', async () => {
    t = await opened(syntheticFamily(30, 8), { width: 600, height: 260 });
    const o = t;
    const host = document.createElement('div');
    host.style.cssText = 'width:300px;height:200px';
    document.body.appendChild(host);
    const second = mountSection(o.ctx, o.section, host, { dpr: 1 });
    const mounted = [o.surface, second];
    const all = mounted.map((s) => canvasesOf({ surface: s } as Opened));
    expect(o.surface.element.isConnected && second.element.isConnected).toBe(true);
    o.ctx.dispose();
    for (const s of mounted) expect(s.element.isConnected).toBe(false);
    for (const c of all.flat()) expect([c.width, c.height]).toEqual([0, 0]);
    // Nothing left to run: a draw now and a frame later do nothing, and a second dispose is a no-op.
    const frames = mounted.map((s) => s.stats.frames);
    for (const s of mounted) {
      s.renderNow();
      s.dispose();
    }
    expect(mounted.map((s) => s.stats.frames)).toEqual(frames);
    // The document's listeners are gone: a release on it reaches no dead surface.
    document.dispatchEvent(new PointerEvent('pointerup', { isPrimary: true }));
    host.remove();
  });

  it('a surface cannot mount on a disposed context', async () => {
    t = await opened(syntheticFamily(10, 6));
    const o = t;
    o.ctx.dispose();
    const host = document.createElement('div');
    expect(() => mountSection(o.ctx, o.section, host)).toThrow(/context is disposed/);
  });

  it('disposing a surface leaves its context and its other surfaces drawing', async () => {
    t = await opened(syntheticFamily(30, 8), { width: 600, height: 260 });
    const o = t;
    const host = document.createElement('div');
    host.style.cssText = 'width:300px;height:200px';
    document.body.appendChild(host);
    const second = mountSection(o.ctx, o.section, host, { dpr: 1 });
    o.surface.dispose();
    expect(second.element.isConnected).toBe(true);
    const frames = second.stats.frames;
    o.ctx.viewports.get(o.section)?.panBy(30, 0, 'api');
    o.advance();
    expect(second.stats.frames).toBeGreaterThan(frames);
    // The last frame's main-thread ms, for a host's own numbers.
    expect(second.stats.lastMs).toBeGreaterThanOrEqual(0);
    expect(Number.isFinite(second.stats.lastMs)).toBe(true);
    second.dispose();
    host.remove();
  });
});

it('a host of no size (hidden, collapsed) draws nothing, and its first size draws every pass whole', async () => {
  // A canvas of no size cannot be drawn from: the back target's present threw at the mount.
  t = await opened(syntheticFamily(120, 8), { width: 0, height: 0 });
  const o = t;
  const vp = o.ctx.viewports.get(o.section);
  for (let k = 0; k < 3; k++) {
    vp?.panBy(13, 7, 'pointer');
    settle(o);
  }
  expect(o.surface.stats.frames).toBe(0);
  o.host.style.cssText += ';width:600px;height:260px';
  // The ResizeObserver reports in the next rendering frame, after its animation frames.
  for (let k = 0; k < 2; k++) await new Promise((r) => requestAnimationFrame(r));
  settle(o);
  const [main] = canvasesOf(o);
  expect([main.width, main.height]).toEqual([600, 260]);
  expect(pixelsOf(main).some((v, i) => i % 4 === 3 && v > 0)).toBe(true);
});

describe('a renderer that throws', () => {
  it('is omitted after its first throw and reported, the rest of the frame drawn, frame after frame', async () => {
    let calls = 0;
    const kernel = (KernelRenderers as LayerRenderer[]).find((r) => r.name === 'ranges') as LayerRenderer;
    const bad: LayerRenderer = {
      ...kernel,
      draw: () => {
        calls++;
        throw new Error('boom');
      },
    };
    // The registry entry replaces the kernel's by name: a host's `ranges`.
    t = await opened(allTypes({ length: 400 }), { replace: { ranges: bad }, width: 700, height: 300 });
    const o = t;
    o.ctx.viewports.get(o.section)?.setView(6, 0, 'api');
    settle(o);
    const first = calls;
    expect(first).toBeGreaterThan(0);
    for (let k = 0; k < 4; k++) {
      o.ctx.viewports.get(o.section)?.panBy(30, 0, 'api');
      settle(o);
    }
    const errors = o.ctx.issues.current.value.filter((i) => i.code === 'seqv.renderer.error');
    expect(errors.length).toBeGreaterThan(0);
    // Once per layer, never again: the failed layer is omitted from every later frame.
    expect(calls).toBe(first);
    // The other layers are on the canvas.
    const px = pixelsOf(canvasesOf(o)[0]);
    let painted = 0;
    for (let i = 3; i < px.length; i += 4) if (px[i] > 0) painted++;
    expect(painted).toBeGreaterThan(2000);
  });
});
