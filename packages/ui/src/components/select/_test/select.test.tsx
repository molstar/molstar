/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { flush as flushDom, render, Show } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Floating, MenuPopover, OptionList, Select } from '../select.js';

beforeAll(() => {
  // jsdom has no popover API: a manual popover (`Floating`) stays `display: none`, and nothing in it takes the focus.
  const style = document.createElement('style');
  style.textContent = '[popover] { display: block !important; }';
  document.head.append(style);
});

const key = (el: Element, k: string) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, composed: true, cancelable: true }));

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

  it('uncontrolled (a plain value or none): a pick shows on the button', async () => {
    const container = document.createElement('div');
    document.body.append(container);
    const onChange = vi.fn();
    render(
      () => (
        <Select
          options={[
            ['x', 'Ex'],
            ['y', 'Why'],
          ]}
          value="x"
          onChange={onChange}
        />
      ),
      container,
    );
    const button = container.querySelector('button')!;
    button.click();
    await flush();
    container.querySelectorAll<HTMLButtonElement>('[role=option]')[1].click();
    flushDom();
    expect(onChange).toHaveBeenCalledWith('y');
    expect(button.textContent).toBe('Why');
  });

  it('closed from the list (a pick, Escape, Tab), the focus goes back to its button; its rows are no Tab stops', async () => {
    for (const over of [true, false]) {
      const container = document.createElement('div');
      document.body.append(container);
      render(
        () => (
          <Select
            over={over}
            options={[
              ['a', 'Alpha'],
              ['b', 'Beta'],
            ]}
            value="a"
          />
        ),
        container,
      );
      const button = container.querySelector<HTMLButtonElement>('[aria-haspopup=listbox]')!;
      const options = () => [...container.querySelectorAll<HTMLButtonElement>('[role=option]')];
      button.click();
      await flush();
      expect(options().map((o) => o.tabIndex)).toEqual([-1, -1]);
      expect(document.activeElement).toBe(options()[0]);
      options()[1].focus();
      options()[1].click();
      flushDom();
      expect(container.querySelector('[role=listbox]')).toBeNull();
      expect(document.activeElement).toBe(button);

      for (const k of ['Escape', 'Tab']) {
        button.click();
        await flush();
        key(options()[0], k);
        flushDom();
        expect(container.querySelector('[role=listbox]')).toBeNull();
        expect(document.activeElement).toBe(button);
      }
      container.remove();
    }
  });

  it("a press on another select's open button closes this one's list (no global query for open buttons)", async () => {
    const container = document.createElement('div');
    document.body.append(container);
    const opts = [['a', 'A'] as const];
    render(
      () => (
        <>
          <Select over={false} options={opts} value="a" />
          <Select over={false} options={opts} value="a" />
        </>
      ),
      container,
    );
    const [one, two] = container.querySelectorAll<HTMLButtonElement>('[aria-haspopup=listbox]');
    one.click();
    await flush();
    // The other's button, marked open (as a form row's own select would be).
    two.setAttribute('aria-expanded', 'true');
    two.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));
    flushDom();
    expect(container.querySelector('[role=listbox]')).toBeNull();
  });
});

describe('MenuPopover', () => {
  function mount() {
    const open = state(true);
    const onClose = vi.fn(() => open.set(false));
    const anchor = document.createElement('button');
    const elsewhere = document.createElement('input');
    const container = document.createElement('div');
    document.body.append(anchor, elsewhere, container);
    render(
      () => (
        <Show when={open}>
          {() => (
            <MenuPopover anchor={anchor} onClose={onClose} aria-label="Row">
              <button type="button" role="menuitem">
                Rename
              </button>
            </MenuPopover>
          )}
        </Show>
      ),
      container,
    );
    return { onClose, elsewhere, container };
  }

  it('takes its own Escape (in it), not one pressed elsewhere on the page; Tab closes it', async () => {
    const a = mount();
    await flush();
    key(a.elsewhere, 'Escape');
    expect(a.onClose).not.toHaveBeenCalled();
    key(a.container.querySelector('[role=menuitem]')!, 'Escape');
    expect(a.onClose).toHaveBeenCalledWith(true);

    document.body.innerHTML = '';
    const b = mount();
    await flush();
    key(b.container.querySelector('[role=menuitem]')!, 'Tab');
    expect(b.onClose).toHaveBeenCalledWith(true);
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
  it('places itself again when its anchor or its column resizes (never wider than its column)', async () => {
    const observed: Element[] = [];
    let callback: (() => void) | undefined;
    vi.stubGlobal(
      'ResizeObserver',
      class {
        constructor(fn: () => void) {
          callback = fn;
        }
        observe(el: Element) {
          observed.push(el);
        }
        disconnect() {}
      },
    );
    const anchor = document.createElement('button');
    const column = document.createElement('div');
    document.body.append(column, anchor);
    let width = 200;
    column.getBoundingClientRect = () =>
      ({ top: 0, bottom: 100, left: 5, right: 5 + width, width, height: 100, x: 5, y: 0 }) as DOMRect;
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <Floating anchor={anchor} column={column}>
          <div />
        </Floating>
      ),
      container,
    );
    const floating = container.firstElementChild as HTMLElement;
    await flush();
    expect(observed).toEqual([floating, anchor, column]);
    expect(floating.style.width).toBe('200px');
    width = 150;
    callback!();
    expect(floating.style.width).toBe('150px');
  });

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
