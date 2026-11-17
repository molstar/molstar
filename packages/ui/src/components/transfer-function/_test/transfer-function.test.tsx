/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { flush, flushFrame, pendingFrameCallbacks, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { lostReleases, pointer, stubPointerCapture } from '../../../testing/gestures.js';
import { valueAxis } from '../../histogram/histogram.js';
import {
  applyPalette,
  TransferFunction,
  type TransferPoint,
  type TransferPreset,
  transferColorAt,
} from '../transfer-function.js';

// The transfer function: a point per control point (its value and opacity), stepped with the keys (left and right
// along the axis, up and down in opacity, Delete removes it), dragged on the plot, added by a press; a color handle per
// point under the plot opening its color.

beforeAll(stubPointerCapture);
afterEach(() => {
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

const stats = { min: -2, max: 8, mean: 0, sigma: 1 };

function mount(opts: { presets?: TransferPreset[] } = {}) {
  const points = state<readonly TransferPoint[]>([
    { x: 0, alpha: 0, color: 0xff0000 },
    { x: 1, alpha: 0.05, color: 0xff0000 },
    { x: 2, alpha: 0.05, color: 0x0000ff },
    { x: 3, alpha: 0, color: 0x0000ff },
  ]);
  const changes: (readonly TransferPoint[])[] = [];
  let commits = 0;
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <TransferFunction
        stats={stats}
        histogram={{ min: -2, max: 8, counts: [1, 4, 9, 4, 1] }}
        points={points}
        logMin={1e-3}
        presets={opts.presets}
        palettes={[{ label: 'Greens', colors: [0x00ff00, 0x008800] }]}
        onChange={(p) => {
          changes.push(p);
          points.set(p);
        }}
        onCommit={() => commits++}
      />
    ),
    container,
  );
  flush();
  const sliders = () => [...container.querySelectorAll<HTMLElement>('[data-point]')];
  const handles = () => [...container.querySelectorAll<HTMLButtonElement>('button[aria-haspopup="dialog"]')];
  const plot = () => sliders()[0].closest<HTMLElement>('[data-tip^="Drag a point"]')!;
  const key = (el: HTMLElement, k: string) => {
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
    flush();
  };
  return { container, dispose, points, changes, sliders, handles, plot, key, commits: () => commits };
}

/** Presses a second apart: none is a double click (which removes what it is on). */
const slowPresses = () => {
  let now = 0;
  vi.spyOn(performance, 'now').mockImplementation(() => (now += 1000));
};

describe('the transfer function', () => {
  it('shows a point per control point, read by its value in the unit and its opacity', () => {
    const t = mount();
    expect(t.sliders().map((s) => s.getAttribute('aria-valuetext'))).toEqual([
      '0.00 σ, opacity 0',
      '1.00 σ, opacity 0.050',
      '2.00 σ, opacity 0.050',
      '3.00 σ, opacity 0',
    ]);
    // At 20% of the plot (−2 … 8), opacity 0 at its bottom.
    expect(t.sliders()[0].parentElement!.style.transform).toBe('translate(20%, 100%)');
    expect(t.handles().map((h) => h.getAttribute('aria-label'))).toEqual([
      'Color at 0.00 σ',
      'Color at 1.00 σ',
      'Color at 2.00 σ',
      'Color at 3.00 σ',
    ]);
  });

  it('draws the curve, the area in the points colors and the color bar', () => {
    const t = mount();
    const stops = [...t.container.querySelectorAll('linearGradient stop')];
    expect(stops.map((s) => [s.getAttribute('offset'), s.getAttribute('stop-color')])).toEqual([
      ['0.2000', '#ff0000'],
      ['0.3000', '#ff0000'],
      ['0.4000', '#0000ff'],
      ['0.5000', '#0000ff'],
    ]);
    const bar = t.handles()[0].parentElement!.parentElement!.firstElementChild as HTMLElement;
    expect(bar.style.background).toContain('linear-gradient');
  });

  it('steps a point with the keys (opacity on the log axis), kept between its neighbours', () => {
    const t = mount();
    t.key(t.sliders()[1], 'ArrowUp');
    expect(t.points.value[1].alpha).toBeCloseTo(0.0625, 6);
    t.key(t.sliders()[1], 'ArrowRight');
    expect(t.points.value[1].x).toBeCloseTo(1.05, 6);
    for (let i = 0; i < 40; i++) t.key(t.sliders()[1], 'ArrowRight');
    expect(t.points.value[1].x).toBeCloseTo(1.99, 6);
    expect(t.commits()).toBe(0);
    t.sliders()[1].dispatchEvent(new KeyboardEvent('keyup', { key: 'ArrowRight', bubbles: true }));
    expect(t.commits()).toBe(1);
  });

  it('removes a point with Delete, two staying', () => {
    const t = mount();
    t.key(t.sliders()[0], 'Delete');
    t.key(t.sliders()[0], 'Delete');
    expect(t.points.value).toHaveLength(2);
    t.key(t.sliders()[0], 'Delete');
    expect(t.points.value).toHaveLength(2);
    expect(t.sliders()).toHaveLength(2);
  });

  it('adds a point where the plot is pressed, in the color the curve has there, and drags it', () => {
    slowPresses();
    const t = mount();
    t.plot().getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    // 35 of 100 px: 1.5, between the second and the third point; at the top: opacity 1.
    t.plot().dispatchEvent(pointer('pointerdown', { x: 35, y: 0 }));
    expect(t.points.value).toHaveLength(5);
    expect(t.points.value[2]).toEqual({ x: 1.5, alpha: 1, color: transferColorAt(t.points.value, 1.5) });
    expect(t.points.value[2].color).toBe(0x7f007f);
    t.plot().dispatchEvent(pointer('pointermove', { x: 38, y: 0 }));
    t.plot().dispatchEvent(pointer('pointerup', { x: 38, y: 0 }));
    expect(t.points.value[2].x).toBeCloseTo(1.8, 6);
    expect(t.commits()).toBe(1);
  });

  it('removes a point pressed twice in a row (a double click)', () => {
    const t = mount();
    const at = { x: 0, y: 0 };
    const p = t.sliders()[1];
    p.dispatchEvent(pointer('pointerdown', at));
    p.dispatchEvent(pointer('pointerup', at));
    t.sliders()[1].dispatchEvent(pointer('pointerdown', at));
    expect(t.points.value.map((q) => q.x)).toEqual([0, 2, 3]);
  });

  it("opens a point's color from its handle: a pick recolors that point only, closing ends the step", async () => {
    const t = mount();
    t.handles()[1].click();
    flush();
    await new Promise((r) => setTimeout(r, 0));
    expect(t.handles()[1].getAttribute('aria-expanded')).toBe('true');
    const dialog = document.querySelector('[role="dialog"]')!;
    expect(dialog.getAttribute('aria-label')).toBe('Color at 1.00 σ');
    const input = dialog.querySelector<HTMLInputElement>('input[aria-label="Color"]')!;
    input.value = '#00ff00';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    expect(t.points.value.map((p) => p.color)).toEqual([0xff0000, 0x00ff00, 0x0000ff, 0x0000ff]);
    expect(t.commits()).toBe(0);
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    flush();
    expect(t.commits()).toBe(1);
    expect(document.querySelector('[role="dialog"]')).toBeNull();
  });

  it('removes the point from its color popover', async () => {
    const t = mount();
    t.handles()[2].click();
    flush();
    await new Promise((r) => setTimeout(r, 0));
    [...document.querySelectorAll<HTMLButtonElement>('[role="dialog"] button')]
      .find((b) => b.textContent?.includes('Remove point'))!
      .click();
    flush();
    expect(t.points.value.map((q) => q.x)).toEqual([0, 1, 3]);
    expect(document.querySelector('[role="dialog"]')).toBeNull();
  });

  it('applies a preset, the colors kept where it has none, and marks the one in effect', async () => {
    const preset: TransferPreset = {
      label: 'Wide',
      points: [
        { x: 0, alpha: 0 },
        { x: 3, alpha: 0.1 },
      ],
    };
    const t = mount({ presets: [preset] });
    t.container.querySelector<HTMLButtonElement>('button[aria-label="Curve presets"]')!.click();
    flush();
    await new Promise((r) => setTimeout(r, 0));
    const option = [...document.querySelectorAll<HTMLElement>('[role="option"]')].find((o) =>
      o.textContent?.includes('Wide'),
    )!;
    option.click();
    flush();
    expect(t.points.value).toEqual([
      { x: 0, alpha: 0, color: 0xff0000 },
      { x: 3, alpha: 0.1, color: 0x0000ff },
    ]);
    expect(t.commits()).toBe(1);
  });

  it('closes the step of a drag it is disposed in the middle of, and asks for no more frames', () => {
    const before = pendingFrameCallbacks();
    const t = mount();
    const point = t.sliders()[1];
    point.dispatchEvent(pointer('pointerdown', { x: 5, y: 5 }));
    point.dispatchEvent(pointer('pointermove', { x: 6, y: 6 }));
    expect(t.commits()).toBe(0);
    expect(t.changes).toHaveLength(0);
    t.dispose();
    expect(t.commits()).toBe(1);
    flushFrame();
    expect(t.changes).toHaveLength(0);
    expect(pendingFrameCallbacks()).toBe(before);
  });

  it('a point drag that loses its release leaves no listener, and the next one moves the point', () => {
    slowPresses();
    const t = mount();
    t.plot().getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    const r = lostReleases({
      target: () => t.sliders()[1],
      from: { x: 30, y: 50 },
      to: { x: 35, y: 50 },
      reset: () => {
        t.changes.length = 0;
      },
      worked: () => t.points.value.length === 4 && Math.abs(t.points.value[1].x - 1.5) < 1e-9,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});

describe('applyPalette', () => {
  const axis = valueAxis({ min: 0, max: 10, mean: 0, sigma: 1 });
  const at = (x: number, alpha: number): TransferPoint => ({ x, alpha, color: 0 });

  it('colors each peak in turn, a point at zero taking the next peak (the last after the last)', () => {
    const points = [at(1, 0), at(2, 0.1), at(3, 0), at(4, 0.2), at(5, 0.2), at(6, 0)];
    const colored = applyPalette(points, { label: 'c', colors: [1, 2] }, axis);
    expect(colored.map((p) => p.color)).toEqual([1, 1, 2, 2, 2, 2]);
  });

  it('spreads a sequential palette along the axis', () => {
    const colored = applyPalette(
      [at(0, 0), at(5, 1), at(10, 0)],
      { label: 's', colors: [0x000000, 0xfefefe], sequential: true },
      axis,
    );
    expect(colored.map((p) => p.color)).toEqual([0x000000, 0x7f7f7f, 0xfefefe]);
  });
});

describe('the transfer function, against its host', () => {
  const four = (): TransferPoint[] => [
    { x: 0, alpha: 0, color: 0xff0000 },
    { x: 1, alpha: 0.05, color: 0xff0000 },
    { x: 2, alpha: 0.05, color: 0x0000ff },
    { x: 3, alpha: 0, color: 0x0000ff },
  ];

  it('drags the point a press added even when the host applies the changes later', async () => {
    slowPresses();
    const points = state<readonly TransferPoint[]>(four());
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <TransferFunction
          stats={stats}
          points={points}
          // The host's own commit: the change lands a microtask later.
          onChange={(p) => queueMicrotask(() => points.set(p))}
        />
      ),
      container,
    );
    flush();
    const plot = container.querySelector<HTMLElement>('[data-tip^="Drag a point"]')!;
    plot.getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    plot.dispatchEvent(pointer('pointerdown', { x: 35, y: 50 }));
    // Before the host applied the insert: a move moves the new point, not the third one at its index.
    plot.dispatchEvent(pointer('pointerup', { x: 38, y: 50 }));
    await Promise.resolve();
    await Promise.resolve();
    expect(points.value.map((q) => q.x)).toEqual([0, 1, expect.closeTo(1.8, 6), 2, 3]);
    dispose();
  });

  it('reads the plot where it is at each step: a panel that scrolls under a drag does not shift the value', () => {
    slowPresses();
    const t = mount();
    let top = 0;
    t.plot().getBoundingClientRect = () => ({ left: 0, top, width: 100, height: 100 }) as DOMRect;
    const p = t.sliders()[1];
    p.dispatchEvent(pointer('pointerdown', { x: 10, y: 50 }));
    // The panel scrolls 40px up under the pointer; the pointer follows the point.
    top = -40;
    p.dispatchEvent(pointer('pointerup', { x: 20, y: 10 }));
    // 20 of 100: value 0, kept after the first point (Gap); its height 0.5 of the plot (y 10 at top -40).
    expect(t.points.value[1].x).toBeCloseTo(0 + 0.001 * 10, 6);
    // Mid-height on the log axis (1e-3 to 1, its bottom band zero): 10^(−3 + 3 × 0.46 / 0.96) ≈ 0.0274; the box read at
    // the press (top 0) would have put it at 0.9 of the height, ≈ 0.49.
    expect(t.points.value[1].alpha).toBeCloseTo(0.0274, 3);
    t.dispose();
  });

  it('is drawn again for new statistics; its unit radios and points say what they are', () => {
    const s = state({ ...stats });
    const points = state<readonly TransferPoint[]>(four());
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => <TransferFunction stats={s} points={points} onChange={(p) => points.set(p)} />,
      container,
    );
    flush();
    const slider = () => container.querySelector<HTMLElement>('[data-point]')!;
    expect(slider().getAttribute('aria-valuemax')).toBe('8');
    s.set({ min: -4, max: 16, mean: 0, sigma: 2 });
    flush();
    expect(slider().getAttribute('aria-valuemax')).toBe('8');
    // In σ: 16 is 8 σ with σ 2; the axis' end in absolute values follows too.
    const unit = container.querySelector<HTMLElement>('[role=radio][aria-label="Absolute values"]')!;
    unit.click();
    flush();
    expect(slider().getAttribute('aria-valuemax')).toBe('16');
    expect(container.querySelector('[role=radio][aria-label="σ from the mean"]')).not.toBeNull();
    const help = document.getElementById(slider().getAttribute('aria-describedby')!)!;
    expect(help.textContent).toContain('up and down change its opacity');
    dispose();
  });
});
