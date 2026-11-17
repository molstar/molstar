// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { focusFirstItem, menuItems, menuKeys, rovingFocus, rovingTabs } from '../focus.js';
import { runInScope, ViewScope } from '../scope.js';

afterEach(() => {
  document.body.innerHTML = '';
});

const key = (el: Element, k: string) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));

describe('rovingFocus', () => {
  it('rovingTabs: one Tab stop on the selected tab; both axes move; right to left swaps ← and →', () => {
    const list = document.createElement('div');
    list.innerHTML =
      '<button role="tab">a</button><button role="tab" aria-selected="true">b</button><button role="tab">c</button>';
    document.body.append(list);
    const scope = new ViewScope();
    runInScope(scope, () => rovingTabs(list));
    const [a, b, c] = list.querySelectorAll<HTMLElement>('[role="tab"]');
    expect([a.tabIndex, b.tabIndex, c.tabIndex]).toEqual([-1, 0, -1]);
    b.focus();
    key(b, 'ArrowRight');
    expect(document.activeElement).toBe(c);
    key(c, 'ArrowUp');
    expect(document.activeElement).toBe(b);
    list.setAttribute('dir', 'rtl');
    key(b, 'ArrowRight');
    expect(document.activeElement).toBe(a);
    key(a, 'ArrowLeft');
    expect(document.activeElement).toBe(b);
    key(b, 'ArrowDown');
    expect(document.activeElement).toBe(c);
    scope.dispose();
  });

  it('one axis: the other arrows move nothing; the Tab stop is the last focused item', () => {
    const bar = document.createElement('div');
    bar.innerHTML = '<button class="i">1</button><button class="i">2</button><button class="i">3</button>';
    document.body.append(bar);
    const scope = new ViewScope();
    runInScope(scope, () => rovingFocus(bar, { items: '.i', axis: 'vertical' }));
    const [one, two, three] = bar.querySelectorAll<HTMLElement>('.i');
    one.focus();
    key(one, 'ArrowRight');
    expect(document.activeElement).toBe(one);
    key(one, 'ArrowDown');
    expect(document.activeElement).toBe(two);
    key(two, 'End');
    expect(document.activeElement).toBe(three);
    three.blur();
    expect([one.tabIndex, two.tabIndex, three.tabIndex]).toEqual([-1, -1, 0]);
    scope.dispose();
  });
});

describe('menu items', () => {
  it('a disabled (aria-disabled) item is reached by the arrows; the menu opens on the first usable one', () => {
    const menu = document.createElement('div');
    menu.innerHTML =
      '<button role="menuitem" tabindex="-1" aria-disabled="true">off</button>' +
      '<button role="menuitem" tabindex="-1">on</button><button role="menuitem" disabled>gone</button>';
    document.body.append(menu);
    const [off, on] = menu.querySelectorAll<HTMLElement>('[role="menuitem"]');
    expect(menuItems(menu)).toEqual([off, on]);
    focusFirstItem(menu);
    expect(document.activeElement).toBe(on);
    const up = new KeyboardEvent('keydown', { key: 'ArrowUp', cancelable: true });
    Object.defineProperty(up, 'currentTarget', { value: menu });
    menuKeys(up);
    expect(document.activeElement).toBe(off);
  });
});
