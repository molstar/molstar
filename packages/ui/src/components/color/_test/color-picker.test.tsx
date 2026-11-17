/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { lostReleases, stubPointerCapture } from '../../../testing/gestures.js';
import {
  ColorPalette,
  ColorPicker,
  colorToHex,
  colorToRgb,
  colorWeights,
  fieldColor,
  parseColor,
  shade,
} from '../color-picker.js';

afterEach(() => {
  document.body.innerHTML = '';
});

describe('parseColor', () => {
  it.each([
    ['#fff', 0xffffff],
    ['0090ff', 0x0090ff],
    ['#0090FF', 0x0090ff],
    ['0x3362b2', 0x3362b2],
    // `0x` is a number: never the short form.
    ['0x123', 0x000123],
    ['0xff', 0x0000ff],
    ['123', 0x112233],
    ['22,44,122', 0x162c7a],
    ['22, 44, 122', 0x162c7a],
    ['22 44 122', 0x162c7a],
    ['rgb(11,12,13)', 0x0b0c0d],
    ['RGB( 11 12 13 )', 0x0b0c0d],
    ['rgba(11, 12, 13, 0.5)', 0x0b0c0d],
    ['rgb(11 12 13 / 50%)', 0x0b0c0d],
    ['rgb(11,12,13', 0x0b0c0d],
    ['(11, 12, 13)', 0x0b0c0d],
    ['rgb(100%, 0%, 50%)', 0xff0080],
    ['1.0, 0.5, 0', 0xff8000],
    ['1, 0, 0', 0x010000],
    ['hsl(0, 100%, 50%)', 0xff0000],
    ['hsl(120 100% 25%)', 0x008000],
    ['steelblue', 0x4682b4],
    ['Steel Blue', 0x4682b4],
  ])('%s', (text, color) => {
    expect(parseColor(text)).toBe(color);
  });

  it.each(['', 'nope', '#12345', '300, 0, 0', '1, 2', 'rgb(a, b, c)', '1,2,3,4,5', 'constructor'])(
    'rejects %j',
    (text) => {
      expect(parseColor(text)).toBeUndefined();
    },
  );
});

describe('ColorPicker', () => {
  it('shades, writes channels, and the smooth palette ends in its hues', () => {
    expect(shade(0x000000, 0.5)).toBe(0x808080);
    expect(shade(0xffffff, -1)).toBe(0x000000);
    expect(colorToRgb(0x3362b2)).toBe('51, 98, 178');
    expect(colorToHex(0x0090ff)).toBe('#0090ff');
    expect(colorToHex(0)).toBe('#000000');
    expect(fieldColor(0, 0.5)).toBe(0xff0000);
    expect(fieldColor(1 / 3, 0.5)).toBe(0x00ff00);
    expect(fieldColor(0, 1)).toBe(shade(0xff0000, -0.9));
  });

  it('spreads the weights over the room a color has: black first lightens, white darkens, a mid color both ways', () => {
    const black = colorWeights(0x000000);
    expect(black[7]).toBe(0x000000);
    expect(new Set(black).size).toBe(8);
    const white = colorWeights(0xffffff);
    expect(white[0]).toBe(0xffffff);
    expect(new Set(white).size).toBe(8);
    const blue = colorWeights(0x123456);
    expect(blue).toHaveLength(8);
    expect(new Set(blue).size).toBe(8);
    expect(blue.indexOf(0x123456)).toBeGreaterThan(3);
  });

  it('picks a hue, then its weights; takes typed colors; the current swatch is selected', () => {
    const color = state(0x0090ff);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ColorPicker value={color} onChange={(c) => color.set(c)} />, container);
    const selected = () => container.querySelector('[role=radio][aria-checked=true]')?.getAttribute('aria-label');
    const text = container.querySelector('input') as HTMLInputElement;
    const preview = text.previousElementSibling as HTMLElement;
    expect(preview.style.background).toBe('rgb(0, 144, 255)');
    expect(selected()).toBe('0, 144, 255');
    expect(text.value).toBe('0, 144, 255');
    (container.querySelector('[aria-label="229, 72, 77"]') as HTMLButtonElement).click();
    expect(color.value).toBe(0xe5484d);
    flush();
    expect(selected()).toBe('229, 72, 77');
    expect(text.value).toBe('229, 72, 77');
    text.focus();
    text.value = 'rgb(18, 52, 86';
    text.dispatchEvent(new Event('input'));
    // Being typed in: left as typed, applied and written out when left.
    expect(color.value).toBe(0xe5484d);
    expect(text.value).toBe('rgb(18, 52, 86');
    text.blur();
    text.dispatchEvent(new Event('blur'));
    expect(color.value).toBe(0x123456);
    expect(text.value).toBe('18, 52, 86');
    flush();
    // Not a hue of the palette: only its own weight is marked.
    expect(selected()).toBe('18, 52, 86');
    // The weights are of the color picked; picking one keeps them.
    const weights = () =>
      [...container.querySelectorAll('[aria-label=Weights] [role=radio]')].map((b) => b.getAttribute('aria-label'));
    const row = colorWeights(0x123456);
    expect(weights()).toEqual(row.map(colorToRgb));
    (container.querySelectorAll('[aria-label=Weights] [role=radio]')[1] as HTMLButtonElement).click();
    expect(color.value).toBe(row[1]);
    flush();
    expect(weights()).toEqual(row.map(colorToRgb));
    expect(selected()).toBe(colorToRgb(row[1]));
  });

  it('applies typed text on Enter, not per keystroke', () => {
    const color = state(0x0090ff);
    const onChange = vi.fn((c: number) => color.set(c));
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ColorPicker value={color} onChange={onChange} />, container);
    const text = container.querySelector('input') as HTMLInputElement;
    text.focus();
    for (const typed of ['2', '25', '255', '255,', '255, 0, 0']) {
      text.value = typed;
      text.dispatchEvent(new Event('input'));
    }
    expect(onChange).not.toHaveBeenCalled();
    // An IME's Enter is the composition's.
    text.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', isComposing: true, bubbles: true }));
    expect(onChange).not.toHaveBeenCalled();
    text.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    expect(onChange).toHaveBeenCalledTimes(1);
    expect(color.value).toBe(0xff0000);
    expect(text.value).toBe('255, 0, 0');
    // Left after Enter: nothing new to apply.
    text.dispatchEvent(new Event('blur'));
    expect(onChange).toHaveBeenCalledTimes(1);
    // Not a color: the current one written back.
    text.focus();
    text.value = 'nope';
    text.dispatchEvent(new Event('input'));
    text.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    expect(onChange).toHaveBeenCalledTimes(1);
    expect(text.value).toBe('255, 0, 0');
  });

  it('keeps the weights row after a weight is picked, also when the host applies it later', async () => {
    const color = state(0x123456);
    const container = document.createElement('div');
    document.body.append(container);
    // The host applies a change in a later task.
    render(() => <ColorPicker value={color} onChange={(c) => setTimeout(() => color.set(c), 0)} />, container);
    const swatches = () => [...container.querySelectorAll<HTMLButtonElement>('[aria-label=Weights] [role=radio]')];
    const row = colorWeights(0x123456).map(colorToRgb);
    swatches()[2].click();
    await new Promise((r) => setTimeout(r, 0));
    flush();
    expect(color.value).toBe(colorWeights(0x123456)[2]);
    expect(swatches().map((b) => b.getAttribute('aria-label'))).toEqual(row);
    expect(swatches()[2].getAttribute('aria-checked')).toBe('true');
    // A change from elsewhere moves the row.
    color.set(0xe5484d);
    flush();
    expect(swatches().map((b) => b.getAttribute('aria-label'))).toEqual(colorWeights(0xe5484d).map(colorToRgb));
  });

  it('each row of swatches is one Tab stop; the arrows move and pick', () => {
    const color = state(ColorPalette[0]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ColorPicker value={color} onChange={(c) => color.set(c)} />, container);
    const palette = [...container.querySelectorAll<HTMLButtonElement>('[aria-label=Palette] [role=radio]')];
    expect(container.querySelector('[aria-label=Palette]')?.getAttribute('role')).toBe('radiogroup');
    expect(palette.map((b) => b.tabIndex)).toEqual([0, -1, -1, -1, -1, -1, -1, -1]);
    palette[0].focus();
    palette[0].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true, cancelable: true }));
    expect(color.value).toBe(ColorPalette[1]);
    expect(document.activeElement).toBe(palette[1]);
    flush();
    expect(palette.map((b) => b.tabIndex)).toEqual([-1, 0, -1, -1, -1, -1, -1, -1]);
    const weights = [...container.querySelectorAll<HTMLButtonElement>('[aria-label=Weights] [role=radio]')];
    expect(weights.filter((b) => b.tabIndex === 0)).toHaveLength(1);
    const at = colorWeights(ColorPalette[1]).indexOf(ColorPalette[1]);
    weights[at].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true, cancelable: true }));
    expect(color.value).toBe(colorWeights(ColorPalette[1])[at + 1]);
    flush();
    // The row stays: walked with the arrows.
    expect(weights[at + 1].getAttribute('aria-label')).toBe(colorToRgb(colorWeights(ColorPalette[1])[at + 1]));
    expect(weights[at + 1].tabIndex).toBe(0);
  });

  describe('drags (gestures attach once: `drag`, dom/scoped.ts)', () => {
    beforeAll(stubPointerCapture);

    // jsdom lays nothing out: each strip 100 px wide and tall at the origin.
    const box = (el: Element) => {
      el.getBoundingClientRect = () => ({ left: 0, top: 0, width: 100, height: 100 }) as DOMRect;
    };
    const mount = () => {
      const color = state(0x0090ff);
      const container = document.createElement('div');
      document.body.append(container);
      render(() => <ColorPicker value={color} onChange={(c) => color.set(c)} />, container);
      const [field, greys] = container.querySelectorAll<HTMLElement>(':scope > div > [aria-hidden=true]');
      box(field);
      box(greys);
      return { color, field, greys };
    };

    it('the smooth palette: no listener left by lost releases, and the next drag picks', () => {
      const { color, field } = mount();
      const r = lostReleases({
        target: () => field,
        from: { x: 10, y: 50 },
        to: { x: 60, y: 30 },
        reset: () => color.set(0),
        worked: () => color.value === fieldColor(0.6, 0.3),
      });
      expect(r.grown).toBe(0);
      expect(r.worked).toBe(true);
    });

    it('the greys: no listener left by lost releases, and the next drag picks', () => {
      const { color, greys } = mount();
      const r = lostReleases({
        target: () => greys,
        from: { x: 10, y: 50 },
        to: { x: 99, y: 50 },
        reset: () => color.set(0x123456),
        worked: () => color.value === 0xffffff,
      });
      expect(r.grown).toBe(0);
      expect(r.worked).toBe(true);
    });
  });
});
