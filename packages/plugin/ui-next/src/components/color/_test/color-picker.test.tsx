/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { ColorPicker, parseHex, shade } from '../color-picker.js';

describe('ColorPicker', () => {
  it('parses hex and shades', () => {
    expect(parseHex('#fff')).toBe(0xffffff);
    expect(parseHex('0090ff')).toBe(0x0090ff);
    expect(parseHex('nope')).toBeUndefined();
    expect(shade(0x000000, 0.5)).toBe(0x808080);
    expect(shade(0xffffff, -1)).toBe(0x000000);
  });

  it('picks from the palette, types hex, lightens; the current swatch is selected', () => {
    const color = observableState(0x0090ff);
    const container = document.createElement('div');
    render(() => <ColorPicker value={color.asReadonly} onChange={(c) => color.set(c)} />, container);
    const selected = () => container.querySelector('[role=option][aria-selected=true]')?.getAttribute('aria-label');
    expect(selected()).toBe('#0090ff');
    (container.querySelector('[aria-label="#e5484d"]') as HTMLButtonElement).click();
    expect(color.value).toBe(0xe5484d);
    expect(selected()).toBe('#e5484d');
    const hex = container.querySelector('input') as HTMLInputElement;
    hex.value = '#123456';
    hex.dispatchEvent(new Event('input'));
    expect(color.value).toBe(0x123456);
    expect(selected()).toBeUndefined();
    const lighten = [...container.querySelectorAll('button')].find((b) => b.textContent === 'Lighten')!;
    lighten.click();
    expect(color.value).toBe(shade(0x123456, 0.15));
  });
});
