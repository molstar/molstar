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
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import type {
  VolumeLevel,
  VolumeLevelChange,
  VolumeLevelsModel,
  VolumeLevelUnit,
} from '../../../model/volume-levels.js';
import { VolumeLevels } from '../volume-levels.js';

// The levels view: a marker per level on the map's histogram, stepped with the keys in σ; a row per level whose value
// is typed in the unit chosen, the other unit beside it.

beforeAll(() => {
  stubPointerCapture();
  (globalThis as { CSS?: unknown }).CSS ??= { escape: (s: string) => s };
});
afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const levels = observableState<readonly VolumeLevel[]>([
    { id: 'plus', label: 'Fo-Fc(+ve)', value: 3, color: 0x00aa00, opacity: 1, hidden: false },
    { id: 'minus', label: 'Fo-Fc(-ve)', value: -3, color: 0xaa0000, opacity: 1, hidden: false },
  ]);
  const unit = observableState<VolumeLevelUnit>('sigma');
  const changes: [string, VolumeLevelChange][] = [];
  let settled = 0;
  const model: VolumeLevelsModel = {
    title: 'Fo-Fc',
    stats: { min: -10, max: 10, mean: 0, sigma: 1 },
    histogram: observableState({ min: -10, max: 10, counts: [1, 5, 20, 5, 1] }).asReadonly,
    levels: levels.asReadonly,
    unit: unit.asReadonly,
    setUnit: (u) => unit.set(u),
    change: (id, c) => {
      changes.push([id, c]);
      levels.set(levels.value.map((l) => (l.id === id ? { ...l, ...c } : l)));
    },
    settle: () => settled++,
  };
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <VolumeLevels model={model} />, container);
  return { container, dispose, changes, unit, settled: () => settled };
}

describe('the levels view', () => {
  it('shows a marker per level where it is, and steps it with the arrow keys in σ', () => {
    const t = mount();
    const markers = [...t.container.querySelectorAll<HTMLElement>('[role="slider"]')];
    expect(markers.map((m) => m.getAttribute('aria-valuetext'))).toEqual(['3.00 σ', '−3.00 σ']);
    // At 65% and 35% of the plot (−10 … 10).
    expect((markers[0].parentElement as HTMLElement).style.transform).toBe('translateX(65%)');
    markers[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowLeft', shiftKey: true, bubbles: true }));
    expect(t.changes).toEqual([['minus', { value: -4 }]]);
    expect(markers[1].getAttribute('aria-valuetext')).toBe('−4.00 σ');
  });

  it('takes a typed value in the unit chosen, showing the other beside it', () => {
    const t = mount();
    const input = t.container.querySelector<HTMLInputElement>('input[aria-label="Fo-Fc(+ve) iso value"]')!;
    expect(input.value).toBe('3.00');
    t.unit.set('absolute');
    expect(input.value).toBe('3');
    expect(t.container.textContent).toContain('3.00 σ');
    input.value = '2.5';
    input.dispatchEvent(new Event('change', { bubbles: true }));
    expect(t.changes.at(-1)).toEqual(['plus', { value: 2.5 }]);
    expect(t.settled()).toBe(1);
  });

  it('opens a quick editor from the swatch: a row of colors applied with one click, the whole picker behind More', async () => {
    const t = mount();
    const swatch = t.container.querySelector<HTMLButtonElement>('button[aria-expanded]')!;
    swatch.click();
    await Promise.resolve();
    const quick = t.container.querySelector('[role="listbox"][aria-label="Fo-Fc(+ve) color"]')!;
    const colors = [...quick.querySelectorAll<HTMLButtonElement>('[role="option"]')];
    expect(colors.length).toBeGreaterThan(5);
    colors[0].click();
    expect(t.changes.at(-1)).toEqual(['plus', { color: 0x3362b2 }]);
    expect(t.settled()).toBe(1);
    expect(colors[0].getAttribute('aria-selected')).toBe('true');
    // The field picker only once More is opened.
    expect(t.container.querySelector('[aria-label="Palette"]')).toBeNull();
    t.container.querySelector<HTMLButtonElement>('button[title="More colors"]')!.click();
    await Promise.resolve();
    expect(t.container.querySelector('[aria-label="Palette"]')).not.toBeNull();
  });

  it('closes the step of a drag it is disposed in the middle of, and asks for no more frames', () => {
    // The drag's step waits for the frame (a frame state, published in the page's one frame tick); disposed, it is
    // dropped, and no frame callback is left.
    const before = pendingFrameCallbacks();
    const t = mount();
    const marker = t.container.querySelector<HTMLElement>('[data-level]')!;
    marker.dispatchEvent(new MouseEvent('pointerdown', { bubbles: true, buttons: 1, clientX: 5 }));
    const pressed = t.changes.length;
    marker.dispatchEvent(new MouseEvent('pointermove', { bubbles: true, buttons: 1, clientX: 6 }));
    expect(t.settled()).toBe(0);
    expect(t.changes).toHaveLength(pressed);

    t.dispose();
    expect(t.settled()).toBe(1);
    flushFrame();
    expect(t.changes).toHaveLength(pressed);
    expect(pendingFrameCallbacks()).toBe(before);
  });

  it('a marker drag that loses its release leaves no listener, and the next one moves the level', () => {
    const t = mount();
    const plot = t.container.querySelector<HTMLElement>('[data-level]')!.closest<HTMLElement>('section')!
      .children[1] as HTMLElement;
    plot.getBoundingClientRect = () => ({ left: 0, width: 200 }) as DOMRect;
    const marker = () => t.container.querySelector<HTMLElement>('[data-level="plus"]')!;
    const r = lostReleases({
      target: marker,
      from: { x: 130, y: 0 },
      to: { x: 180, y: 0 },
      reset: () => (t.changes.length = 0),
      // 180 of 200 px: 8 on −10 … 10.
      worked: () => t.changes.at(-1)?.[0] === 'plus' && t.changes.at(-1)?.[1].value === 8,
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
