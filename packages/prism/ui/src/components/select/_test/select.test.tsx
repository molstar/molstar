/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { Floating, Select } from '../select.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

describe('Select', () => {
  it('expands inline on the current option, moves with arrows, picks, closes on Escape and outside', async () => {
    const value = observableState('b');
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <Select
          aria-label="Letter"
          options={[
            ['a', 'Alpha'],
            ['b', 'Beta'],
            ['c', 'Gamma'],
          ]}
          value={value}
          onChange={(v) => value.set(v)}
        />
      ),
      container,
    );
    const button = container.querySelector('button[aria-haspopup=listbox]') as HTMLButtonElement;
    expect(button.textContent).toBe('Beta');
    button.click();
    await flush();
    const options = () => [...container.querySelectorAll('[role=option]')] as HTMLButtonElement[];
    expect(document.activeElement).toBe(options()[1]);
    options()[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
    expect(document.activeElement).toBe(options()[2]);
    options()[2].click();
    expect(value.value).toBe('c');
    expect(button.textContent).toBe('Gamma');
    expect(container.querySelector('[role=listbox]')).toBeNull();

    button.click();
    await flush();
    options()[0].dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(container.querySelector('[role=listbox]')).toBeNull();

    button.click();
    document.body.dispatchEvent(new Event('pointerdown', { bubbles: true }));
    expect(container.querySelector('[role=listbox]')).toBeNull();
    container.remove();
  });
});

describe('Floating', () => {
  it('places itself on page scrolls, not on scrolls of its own content, and writes its min width once', async () => {
    const anchor = document.createElement('button');
    document.body.append(anchor);
    const rect = vi.spyOn(anchor, 'getBoundingClientRect');
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <Floating anchor={anchor}>
          <div class="list" />
        </Floating>
      ),
      container,
    );
    await flush();
    expect(rect).toHaveBeenCalledTimes(1);
    const floating = container.firstElementChild as HTMLElement;
    const minWidth = vi.spyOn(floating.style, 'minWidth', 'set');

    floating.querySelector('.list')!.dispatchEvent(new Event('scroll'));
    expect(rect).toHaveBeenCalledTimes(1);

    document.dispatchEvent(new Event('scroll'));
    expect(rect).toHaveBeenCalledTimes(2);
    expect(minWidth).not.toHaveBeenCalled();

    dispose();
    document.dispatchEvent(new Event('scroll'));
    expect(rect).toHaveBeenCalledTimes(2);
    anchor.remove();
    container.remove();
  });
});
