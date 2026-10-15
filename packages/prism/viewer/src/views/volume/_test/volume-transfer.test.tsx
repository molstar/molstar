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
import { afterEach, describe, expect, it } from 'vitest';
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

afterEach(() => {
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
    settle: () => {},
    colorsAtPeaks: () => {},
    palettes: [],
    applyPreset: () => {},
  };
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <VolumeTransfer model={model} />, container);
  const sliders = () => [...container.querySelectorAll<HTMLElement>('[data-point]')];
  const key = (el: HTMLElement, k: string) => el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  return { container, points, colors, changes, sliders, key };
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
});
