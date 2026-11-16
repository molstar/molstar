/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { flush, flushFrame, pendingFrameCallbacks, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { lostReleases, pointer, stubPointerCapture } from '../../../testing/gestures.js';
import type { ValueUnit } from '../../histogram/histogram.js';
import { type IsoLevel, IsoLevelControl } from '../iso-level.js';

// The iso levels: a marker per level on the values' histogram, stepped with the keys in σ, read in the unit chosen;
// dragged on the plot (the marker under the pointer, else the nearest), once a frame.

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const levels = state<readonly IsoLevel[]>([
    { id: 'plus', label: 'Fo-Fc(+)', value: 3, color: 0x00aa00 },
    { id: 'minus', label: 'Fo-Fc(−)', value: -3, color: 0xaa0000, hidden: true },
  ]);
  const unit = state<ValueUnit>('sigma');
  const changes: [string, number][] = [];
  let commits = 0;
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <IsoLevelControl
        title="Fo-Fc"
        stats={{ min: -10, max: 10, mean: 0, sigma: 1 }}
        histogram={{ min: -10, max: 10, counts: [1, 5, 20, 5, 1] }}
        levels={levels}
        unit={unit}
        onChange={(id, value) => {
          changes.push([id, value]);
          levels.set(levels.value.map((l) => (l.id === id ? { ...l, value } : l)));
        }}
        onCommit={() => commits++}
      />
    ),
    container,
  );
  flush();
  const markers = () => [...container.querySelectorAll<HTMLElement>('[role="slider"]')];
  const plot = () => markers()[0].closest<HTMLElement>('section')!.children[1] as HTMLElement;
  const key = (el: HTMLElement, k: string, shiftKey = false) => {
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, shiftKey, bubbles: true }));
    flush();
  };
  return { container, dispose, levels, changes, unit, markers, plot, key, commits: () => commits };
}

describe('the iso level control', () => {
  it('has a header with its title and the unit, the histogram, the mean and the axis in the unit', () => {
    const t = mount();
    expect(t.container.querySelector('section')!.getAttribute('aria-label')).toBe('Fo-Fc');
    expect(t.container.textContent).toContain('Fo-Fc');
    expect(t.container.querySelector('path')!.getAttribute('d')).toMatch(/^M0\.0000 1L/);
    expect(t.container.textContent).toContain('−10.00 σ');
    t.unit.set('absolute');
    flush();
    expect(t.container.textContent).toContain('−10');
    expect(t.container.textContent).not.toContain('−10.00 σ');
  });

  it('shows a marker per level where it is, in its color, a hidden one faded', () => {
    const t = mount();
    const [plus, minus] = t.markers();
    expect(plus.getAttribute('aria-label')).toBe('Fo-Fc(+) iso value');
    // At 65% and 35% of the plot (−10 … 10).
    expect(plus.parentElement!.style.transform).toBe('translateX(65%)');
    expect(minus.parentElement!.style.transform).toBe('translateX(35%)');
    expect(plus.parentElement!.style.getPropertyValue('--msui-level')).toBe('#00aa00');
    expect(minus.className).toMatch(/off/);
    expect(plus.className).not.toMatch(/off/);
  });

  it('steps a marker with the keys in σ, reads it in the unit, and a key released ends the step', () => {
    const t = mount();
    const minus = t.markers()[1];
    expect(minus.getAttribute('aria-valuetext')).toBe('−3.00 σ');
    expect(minus.getAttribute('aria-valuenow')).toBe('-3');
    t.key(minus, 'ArrowLeft', true);
    t.key(minus, 'ArrowRight');
    expect(t.changes).toEqual([
      ['minus', -4],
      ['minus', -3.9],
    ]);
    expect(t.commits()).toBe(0);
    minus.dispatchEvent(new KeyboardEvent('keyup', { key: 'ArrowRight', bubbles: true }));
    expect(t.commits()).toBe(1);
    t.key(minus, 'Home');
    expect(t.changes.at(-1)).toEqual(['minus', -10]);
    t.key(minus, 'PageUp');
    expect(t.changes.at(-1)).toEqual(['minus', -9]);
    t.unit.set('absolute');
    flush();
    expect(minus.getAttribute('aria-valuetext')).toBe('−9');
    expect(minus.getAttribute('aria-valuemax')).toBe('10');
  });

  it('drags the nearest marker from a press anywhere on the plot, once a frame', () => {
    const t = mount();
    const plot = t.plot();
    plot.getBoundingClientRect = () => ({ left: 0, top: 0, width: 200, height: 48 }) as DOMRect;
    // 50 of 200 px: −5, nearest the − level (at 70 px).
    plot.dispatchEvent(pointer('pointerdown', { x: 50, y: 10 }));
    expect(t.changes).toEqual([['minus', -5]]);
    plot.dispatchEvent(pointer('pointermove', { x: 40, y: 10 }));
    plot.dispatchEvent(pointer('pointermove', { x: 30, y: 10 }));
    expect(t.changes).toHaveLength(1);
    flushFrame();
    expect(t.changes).toEqual([
      ['minus', -5],
      ['minus', -7],
    ]);
    plot.dispatchEvent(pointer('pointerup', { x: 20, y: 10 }));
    expect(t.changes.at(-1)).toEqual(['minus', -8]);
    expect(t.commits()).toBe(1);
  });

  it('closes the step of a drag it is disposed in the middle of, and asks for no more frames', () => {
    const before = pendingFrameCallbacks();
    const t = mount();
    const marker = t.markers()[0];
    marker.dispatchEvent(pointer('pointerdown', { x: 5, y: 0 }));
    const pressed = t.changes.length;
    marker.dispatchEvent(pointer('pointermove', { x: 6, y: 0 }));
    expect(t.commits()).toBe(0);
    t.dispose();
    expect(t.commits()).toBe(1);
    flushFrame();
    expect(t.changes).toHaveLength(pressed);
    expect(pendingFrameCallbacks()).toBe(before);
  });

  it('a marker drag that loses its release leaves no listener, and the next one moves the level', () => {
    const t = mount();
    t.plot().getBoundingClientRect = () => ({ left: 0, width: 200 }) as DOMRect;
    const r = lostReleases({
      target: () => t.markers()[0],
      from: { x: 130, y: 0 },
      to: { x: 180, y: 0 },
      reset: () => {
        t.changes.length = 0;
      },
      // 180 of 200 px: 8 on −10 … 10.
      worked: () => t.changes.at(-1)?.[0] === 'plus' && t.changes.at(-1)?.[1] === 8,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
