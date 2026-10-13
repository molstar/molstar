/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { ColorPicker, colorToRgb, colorWeights, fieldColor, parseColor, shade } from '../color-picker.js';

describe('parseColor', () => {
  it.each([
    ['#fff', 0xffffff],
    ['0090ff', 0x0090ff],
    ['#0090FF', 0x0090ff],
    ['0x3362b2', 0x3362b2],
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
    const color = observableState(0x0090ff);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ColorPicker value={color.asReadonly} onChange={(c) => color.set(c)} />, container);
    const selected = () => container.querySelector('[role=option][aria-selected=true]')?.getAttribute('aria-label');
    const text = container.querySelector('input') as HTMLInputElement;
    expect(selected()).toBe('0, 144, 255');
    expect(text.value).toBe('0, 144, 255');
    (container.querySelector('[aria-label="229, 72, 77"]') as HTMLButtonElement).click();
    expect(color.value).toBe(0xe5484d);
    expect(selected()).toBe('229, 72, 77');
    expect(text.value).toBe('229, 72, 77');
    text.focus();
    text.value = 'rgb(18, 52, 86';
    text.dispatchEvent(new Event('input'));
    expect(color.value).toBe(0x123456);
    // Being typed in: left as typed, written out when left.
    expect(text.value).toBe('rgb(18, 52, 86');
    text.blur();
    text.dispatchEvent(new Event('blur'));
    expect(text.value).toBe('18, 52, 86');
    // Not a hue of the palette: only its own weight is marked.
    expect(selected()).toBe('18, 52, 86');
    // The weights are of the color picked; picking one keeps them.
    const weights = () =>
      [...container.querySelectorAll('[aria-label=Weights] [role=option]')].map((b) => b.getAttribute('aria-label'));
    const row = colorWeights(0x123456);
    expect(weights()).toEqual(row.map(colorToRgb));
    (container.querySelectorAll('[aria-label=Weights] [role=option]')[1] as HTMLButtonElement).click();
    expect(color.value).toBe(row[1]);
    expect(weights()).toEqual(row.map(colorToRgb));
    expect(selected()).toBe(colorToRgb(row[1]));
  });
});
