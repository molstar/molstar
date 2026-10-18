/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { flushFrame } from '@molstar/prism-ui/dom/scope';
import { pendingFrameCallbacks } from '@molstar/prism-ui/reactive/frame-state';
import { lostReleases, stubPointerCapture } from '@molstar/prism-ui/testing/gestures';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import type { VolumeLevelUnit } from '../../../model/volume-levels.js';
import type {
  TransferChange,
  TransferColors,
  TransferPoint,
  VolumeTransferModel,
} from '../../../model/volume-transfer.js';
import { VolumeTransfer } from '../volume-transfer.js';

// The transfer function's view: a point per control point, stepped with the keys (left and right along the axis, up
// and down in opacity, Delete removes it), a stop on a point moving with it; the strip's stops.

beforeAll(stubPointerCapture);
afterEach(() => {
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

function mount() {
  const points = observableState<readonly TransferPoint[]>([
    { x: 0.2, alpha: 0 },
    { x: 0.3, alpha: 0.05 },
    { x: 0.4, alpha: 0.05 },
    { x: 0.5, alpha: 0 },
  ]);
  const colors = observableState<TransferColors>({
    kind: 'value',
    interpolate: true,
    domain: [0, 1],
    stops: [
      { x: 0.3, color: 0xff0000 },
      { x: 0.4, color: 0xff0000 },
    ],
  });
  const changes: TransferChange[] = [];
  const settled: number[] = [];
  const subscribed = vi.spyOn(colors.asReadonly, 'subscribe');
  const model: VolumeTransferModel = {
    stats: { min: -2, max: 8, mean: 0, sigma: 1 },
    histogram: observableState({ min: -2, max: 8, counts: [1, 4, 9, 4, 1] }).asReadonly,
    points: points.asReadonly,
    colors: colors.asReadonly,
    logMin: 1e-3,
    presets: () => [],
    unit: observableState<VolumeLevelUnit>('sigma').asReadonly,
    setUnit: () => {},
    change: (c) => {
      changes.push(c);
      if (c.points) points.set(c.points);
      if (c.stops && colors.value.kind === 'value') colors.set({ ...colors.value, stops: c.stops });
    },
    settle: () => void settled.push(1),
    colorsAtPeaks: () => {},
    palettes: [],
    applyPreset: () => {},
  };
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <VolumeTransfer model={model} />, container);
  const sliders = () => [...container.querySelectorAll<HTMLElement>('[data-point]')];
  const key = (el: HTMLElement, k: string) => el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  return { container, dispose, points, colors, changes, settled, subscribed, sliders, key };
}

describe('the transfer function view', () => {
  it('shows a point per control point, named by its value and opacity', () => {
    const t = mount();
    expect(t.sliders().map((s) => s.getAttribute('aria-label'))).toEqual([
      'Point at 0.00 σ, opacity 0',
      'Point at 1.00 σ, opacity 0.050',
      'Point at 2.00 σ, opacity 0.050',
      'Point at 3.00 σ, opacity 0',
    ]);
    expect(t.container.querySelectorAll('[data-stop]')).toHaveLength(2);
  });

  it('steps a point with the keys, and a stop on it moves with it', () => {
    const t = mount();
    t.key(t.sliders()[1], 'ArrowUp');
    expect(t.points.value[1].alpha).toBeCloseTo(0.0625, 6);
    t.key(t.sliders()[1], 'ArrowRight');
    expect(t.points.value[1].x).toBeCloseTo(0.305, 6);
    const c = t.colors.value;
    expect(c.kind === 'value' && c.stops.map((s) => s.x)).toEqual([0.305, 0.4]);
  });

  it('keeps a point between its neighbours, and removes one with Delete (two stay)', () => {
    const t = mount();
    for (let i = 0; i < 40; i++) t.key(t.sliders()[1], 'ArrowRight');
    expect(t.points.value[1].x).toBeCloseTo(0.399, 6);
    t.key(t.sliders()[0], 'Delete');
    t.key(t.sliders()[0], 'Delete');
    expect(t.points.value).toHaveLength(2);
    t.key(t.sliders()[0], 'Delete');
    expect(t.points.value).toHaveLength(2);
  });

  it('reads the colors once for all its bindings (the gradient has dozens of stops)', () => {
    const t = mount();
    expect(t.container.querySelectorAll('stop').length).toBeGreaterThan(60);
    expect(t.subscribed.mock.calls.length).toBeLessThan(10);
  });

  it('closes the step of a drag it is disposed in the middle of, and asks for no more frames', () => {
    // The drag's step waits for the frame (a frame state, published in the page's one frame tick); disposed, it is
    // dropped, and no frame callback is left.
    const before = pendingFrameCallbacks();
    const t = mount();
    const point = t.sliders()[1];
    point.dispatchEvent(new MouseEvent('pointerdown', { bubbles: true, buttons: 1, clientX: 5, clientY: 5 }));
    point.dispatchEvent(new MouseEvent('pointermove', { bubbles: true, buttons: 1, clientX: 6, clientY: 6 }));
    expect(t.settled).toEqual([]);
    expect(t.changes).toHaveLength(0);

    t.dispose();
    expect(t.settled).toEqual([1]);
    flushFrame();
    expect(t.changes).toHaveLength(0);
    expect(pendingFrameCallbacks()).toBe(before);
  });

  /** Presses a second apart: none is a double click (which removes what it is on). */
  const slowPresses = () => {
    let now = 0;
    vi.spyOn(performance, 'now').mockImplementation(() => (now += 1000));
  };

  it('a point drag that loses its release leaves no listener, and the next one moves the point', () => {
    slowPresses();
    const t = mount();
    const plot = t.sliders()[0].closest<HTMLElement>('[title^="Drag a point"]')!;
    plot.getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    const r = lostReleases({
      target: () => t.sliders()[1],
      from: { x: 30, y: 50 },
      to: { x: 35, y: 50 },
      reset: () => (t.changes.length = 0),
      worked: () => t.points.value.length === 4 && Math.abs(t.points.value[1].x - 0.35) < 1e-9,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });

  it('a stop drag that loses its release leaves no listener, and the next one moves the stop', () => {
    slowPresses();
    const t = mount();
    const strip = t.container.querySelector<HTMLElement>('[data-stop]')!.parentElement!.parentElement!;
    // Wide enough that 0.45 is not within a snap (6 px) of a point.
    strip.getBoundingClientRect = () => ({ left: 0, top: 0, width: 1000, height: 10 }) as DOMRect;
    const stop = () => t.container.querySelector<HTMLElement>('[data-stop="1"]')!;
    const r = lostReleases({
      target: stop,
      from: { x: 400, y: 5 },
      to: { x: 450, y: 5 },
      reset: () => (t.changes.length = 0),
      worked: () => {
        const c = t.colors.value;
        return c.kind === 'value' && c.stops.length === 2 && Math.abs(c.stops[1].x - 0.45) < 1e-9;
      },
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
