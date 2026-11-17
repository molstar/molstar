/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import {
  formatValue,
  fromSigma,
  histogramPath,
  sigmaOf,
  UnitSwitch,
  type ValueUnit,
  valueAxis,
  valueNumber,
  valueText,
} from '../histogram.js';

afterEach(() => {
  document.body.innerHTML = '';
});

const stats = { min: -2, max: 8, mean: 1, sigma: 2 };

describe('the value axis', () => {
  it('maps values to [0, 1] and back, clamped', () => {
    const axis = valueAxis(stats);
    expect(axis.at(-2)).toBe(0);
    expect(axis.at(3)).toBe(0.5);
    expect(axis.at(100)).toBe(1);
    expect(axis.valueAt(0.25)).toBe(0.5);
    expect(axis.valueAt(-1)).toBe(-2);
    expect(axis.clamp(9)).toBe(8);
  });

  it('keeps an empty range from dividing by zero', () => {
    const axis = valueAxis({ min: 1, max: 1, mean: 1, sigma: 0 });
    expect(axis.span).toBe(1);
    expect(axis.at(1)).toBe(0);
  });
});

describe('values in a unit', () => {
  it('writes σ from the mean (two decimals) and absolute values (three digits), with a typographic minus', () => {
    expect(sigmaOf(4, stats)).toBe(1.5);
    expect(fromSigma(1.5, stats)).toBe(4);
    expect(formatValue('sigma', 4, stats)).toBe('1.50');
    expect(formatValue('sigma', -2, stats)).toBe('−1.50');
    expect(formatValue('absolute', 0.073421, stats)).toBe('0.0734');
    expect(formatValue('absolute', -1234.5, stats)).toBe('−1230');
    expect(valueText('sigma', 4, stats)).toBe('1.50 σ');
    expect(valueText('absolute', 4, stats)).toBe('4');
    expect(valueNumber('sigma', -2, stats)).toBe(-1.5);
  });
});

describe('histogramPath', () => {
  it('draws a step per bin, placed by value on the axis, log counts', () => {
    const axis = valueAxis({ min: 0, max: 4, mean: 2, sigma: 1 });
    // Two bins over [2, 4]: the right half of the axis.
    const d = histogramPath({ min: 2, max: 4, counts: [0, 9] }, axis);
    expect(d).toBe('M0.5000 1L0.5000 1.0000L0.7500 1.0000L0.7500 0.0000L1.0000 0.0000L1.0000 1Z');
    expect(histogramPath(undefined, axis)).toBe('');
    expect(histogramPath({ min: 0, max: 1, counts: [] }, axis)).toBe('');
  });
});

describe('UnitSwitch', () => {
  it('is a radio group of σ and abs: a click or the arrows pick, one Tab stop', () => {
    const unit = state<ValueUnit>('sigma');
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <UnitSwitch unit={unit} />, container);
    const radios = [...container.querySelectorAll<HTMLButtonElement>('[role="radio"]')];
    expect(radios.map((r) => r.textContent)).toEqual(['σ', 'abs']);
    expect(radios.map((r) => r.getAttribute('aria-checked'))).toEqual(['true', 'false']);
    expect(radios.map((r) => r.tabIndex)).toEqual([0, -1]);
    radios[1].click();
    flush();
    expect(unit.value).toBe('absolute');
    expect(radios.map((r) => r.getAttribute('aria-checked'))).toEqual(['false', 'true']);
    radios[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    flush();
    expect(unit.value).toBe('sigma');
    expect(radios.map((r) => r.tabIndex)).toEqual([0, -1]);
    dispose();
  });
});
