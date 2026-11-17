/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush as flushDom, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Floating, OptionList, Select } from '../select.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

describe('Select', () => {
  it('inline (`over={false}`): expands on the current option, moves with arrows, picks, closes on Escape and outside', async () => {
    const value = state('b');
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <Select
          over={false}
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
    expect(button.getAttribute('aria-expanded')).toBe('false');
    button.click();
    await flush();
    expect(button.getAttribute('aria-expanded')).toBe('true');
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
  });

  it('opens over its field by default, in the top layer, and the layout keeps its place', async () => {
    const value = state('b');
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <Select
          aria-label="Letter"
          options={[
            ['a', 'Alpha'],
            ['b', 'Beta'],
          ]}
          value={value}
          onChange={(v) => value.set(v)}
        />
      ),
      container,
    );
    const select = container.firstElementChild as HTMLElement;
    const button = select.querySelector('button') as HTMLButtonElement;
    button.click();
    await flush();
    const floating = select.querySelector<HTMLElement>('[popover]')!;
    expect(floating.dataset.placement).toBe('over');
    expect(floating.querySelector('[role=listbox]')).not.toBeNull();
    // Only the floating layer was added: nothing in the flow under the field.
    expect([...select.children].filter((c) => c !== button && c !== floating && c.nodeType === 1)).toEqual([]);
    floating.querySelectorAll<HTMLButtonElement>('[role=option]')[0].click();
    expect(value.value).toBe('a');
    expect(select.querySelector('[popover]')).toBeNull();
  });

  it('takes a plain value', () => {
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <Select options={[['x', 'Ex']]} value="x" />, container);
    expect(container.querySelector('button')?.textContent).toBe('Ex');
  });
});

describe('OptionList', () => {
  it('marks the current option through a keyed selector: a change rewrites the two options it concerns', async () => {
    const value = state('a');
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <OptionList
          options={[
            ['a', 'A'],
            ['b', 'B'],
            ['c', 'C'],
          ]}
          value={value}
          onPick={() => {}}
          onClose={() => {}}
        />
      ),
      container,
    );
    const options = [...container.querySelectorAll('[role=option]')];
    const selected = () => options.map((o) => o.getAttribute('aria-selected'));
    expect(selected()).toEqual(['true', 'false', 'false']);
    const writes = vi.fn();
    const observer = new MutationObserver((records) => writes(records.map((r) => r.target)));
    for (const o of options) observer.observe(o, { attributes: true, attributeFilter: ['aria-selected'] });
    value.set('c');
    flushDom();
    expect(selected()).toEqual(['false', 'false', 'true']);
    await Promise.resolve();
    expect(new Set(writes.mock.calls.flat(2))).toEqual(new Set([options[0], options[2]]));
    observer.disconnect();
  });
});

describe('Floating', () => {
  it('over: centers the current option on its anchor, measured once', async () => {
    const anchor = document.createElement('button');
    document.body.append(anchor);
    const rect = (top: number, height: number) =>
      ({ top, bottom: top + height, left: 10, right: 110, width: 100, height, x: 10, y: top }) as DOMRect;
    anchor.getBoundingClientRect = () => rect(100, 20);
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <Floating anchor={anchor} placement="over">
          <div role="option" aria-selected="false" />
          <div role="option" aria-selected="true" />
        </Floating>
      ),
      container,
    );
    const floating = container.firstElementChild as HTMLElement;
    floating.getBoundingClientRect = () => rect(0, 100);
    const [, current] = floating.querySelectorAll<HTMLElement>('[role=option]');
    const read = vi.fn(() => rect(52, 24));
    current.getBoundingClientRect = read;
    await flush();
    // The option's middle (52 + 12 into the list) on the anchor's (100 + 10): the list's top at 46.
    expect(floating.style.top).toBe('46px');
    expect(floating.style.left).toBe('10px');
    document.dispatchEvent(new Event('scroll'));
    expect(read).toHaveBeenCalledTimes(1);
    expect(floating.style.top).toBe('46px');
  });

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
  });

  it('records the side it took, a flip included, for its entrance; its rise offsets are its own', async () => {
    const anchor = document.createElement('button');
    document.body.append(anchor);
    const rect = (top: number, height: number) =>
      ({ top, bottom: top + height, left: 10, right: 110, width: 100, height, x: 10, y: top }) as DOMRect;
    // Near the bottom of the window: a list asked below opens above.
    anchor.getBoundingClientRect = () => rect(window.innerHeight - 30, 20);
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
    const floating = container.firstElementChild as HTMLElement;
    Object.defineProperty(floating, 'offsetHeight', { value: 200 });
    await flush();
    expect(floating.dataset.side).toBe('above');
    anchor.getBoundingClientRect = () => rect(40, 20);
    document.dispatchEvent(new Event('scroll'));
    expect(floating.dataset.side).toBe('below');
    dispose();
    // Custom properties inherit: the floating layer sets both offsets itself before its side's, so a layer inside
    // another (a select in a popover) never rises by its parent's.
    // The source (the generated module keeps its text in the registry); tests run from the repository's root.
    const css = readFileSync(resolve(process.cwd(), 'packages/ui/src/components/select/select.module.css'), 'utf8');
    const base = /\n\.floating \{[^}]*?--msui-rise-x: 0px;\s*--msui-rise-y: 0px;/;
    expect(css).toMatch(base);
  });
});
