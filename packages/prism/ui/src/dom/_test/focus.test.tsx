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
import { openContextMenu } from '../../components/context-menu/context-menu.js';
import { menuKeys, rememberFocus, rovingTabs, trapTab } from '../focus.js';
import { render } from '../render.js';

const key = (el: Element, k: string, init: KeyboardEventInit = {}) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true, ...init }));

function mount(fn: () => any) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(fn, container);
  return {
    container,
    done: () => {
      dispose();
      container.remove();
    },
  };
}

describe('menuKeys', () => {
  it('moves between the usable items by the arrows, wrapping; Home and End jump', () => {
    const { container, done } = mount(() => (
      <div role="menu" onKeyDown={menuKeys}>
        <button type="button" role="menuitem">
          A
        </button>
        <button type="button" role="menuitem" disabled>
          B
        </button>
        <button type="button" role="menuitemcheckbox" aria-checked="true">
          C
        </button>
        <button type="button" role="menuitem">
          D
        </button>
      </div>
    ));
    const [a, , c, d] = [...container.querySelectorAll('button')];
    a.focus();
    key(a, 'ArrowDown');
    expect(document.activeElement).toBe(c);
    key(c, 'ArrowDown');
    expect(document.activeElement).toBe(d);
    key(d, 'ArrowDown');
    expect(document.activeElement).toBe(a);
    key(a, 'ArrowUp');
    expect(document.activeElement).toBe(d);
    key(d, 'Home');
    expect(document.activeElement).toBe(a);
    key(a, 'End');
    expect(document.activeElement).toBe(d);
    done();
  });
});

describe('rovingTabs', () => {
  it('makes a tab list one Tab stop (the selected tab) and moves between tabs by the arrows', async () => {
    const { container, done } = mount(() => (
      <div role="tablist">
        <button type="button" role="tab" aria-selected="false">
          Home
        </button>
        <button type="button" role="tab" aria-selected="true">
          Structure
        </button>
        <button type="button" role="tab" aria-selected="false">
          Settings
        </button>
      </div>
    ));
    const list = container.querySelector<HTMLElement>('[role=tablist]')!;
    const stop = rovingTabs(list);
    const tabs = [...list.querySelectorAll<HTMLElement>('[role=tab]')];
    expect(tabs.map((t) => t.tabIndex)).toEqual([-1, 0, -1]);
    tabs[1].focus();
    key(tabs[1], 'ArrowDown');
    expect(document.activeElement).toBe(tabs[2]);
    expect(tabs.map((t) => t.tabIndex)).toEqual([-1, -1, 0]);
    key(tabs[2], 'ArrowDown');
    expect(document.activeElement).toBe(tabs[0]);
    key(tabs[0], 'End');
    expect(document.activeElement).toBe(tabs[2]);
    // Another tab selected while the focus is elsewhere: it becomes the stop.
    tabs[2].blur();
    tabs[0].setAttribute('aria-selected', 'true');
    tabs[1].setAttribute('aria-selected', 'false');
    await new Promise((r) => setTimeout(r, 0));
    expect(tabs.map((t) => t.tabIndex)).toEqual([0, -1, -1]);
    stop();
    done();
  });
});

describe('trapTab and rememberFocus', () => {
  it('keeps Tab inside a dialog and gives the focus back to what opened it', () => {
    const opener = document.createElement('button');
    document.body.append(opener);
    opener.focus();
    const restore = rememberFocus();
    const { container, done } = mount(() => (
      <div role="dialog" onKeyDown={(e) => trapTab(e, e.currentTarget)}>
        <button type="button">first</button>
        <input />
        <button type="button">last</button>
      </div>
    ));
    const dialog = container.querySelector<HTMLElement>('[role=dialog]')!;
    const [first, last] = [...dialog.querySelectorAll('button')];
    last.focus();
    expect(key(last, 'Tab')).toBe(false);
    expect(document.activeElement).toBe(first);
    expect(key(first, 'Tab', { shiftKey: true })).toBe(false);
    expect(document.activeElement).toBe(last);
    // A Tab between the ends is the browser's.
    first.focus();
    expect(key(first, 'Tab')).toBe(true);
    // Closed with the focus inside: back to the opener.
    restore(dialog);
    expect(document.activeElement).toBe(opener);
    done();
    opener.remove();
  });

  it('leaves a focus that moved on to something else', () => {
    const opener = document.createElement('button');
    const other = document.createElement('input');
    document.body.append(opener, other);
    opener.focus();
    const restore = rememberFocus();
    other.focus();
    restore(null);
    expect(document.activeElement).toBe(other);
    opener.remove();
    other.remove();
  });
});

describe('openContextMenu', () => {
  it('takes the arrows and gives the focus back on Escape', () => {
    const root = document.createElement('div');
    const opener = document.createElement('button');
    document.body.append(root, opener);
    opener.focus();
    const ran: string[] = [];
    openContextMenu(root, { x: 10, y: 10 }, [
      { id: 'a', label: 'Hide', run: () => ran.push('a') },
      { id: 'b', label: 'Focus', run: () => ran.push('b') },
    ]);
    const items = [...root.querySelectorAll<HTMLElement>('[role=menuitem]')];
    expect(document.activeElement).toBe(items[0]);
    key(items[0], 'ArrowDown');
    expect(document.activeElement).toBe(items[1]);
    key(items[1], 'Escape');
    expect(root.querySelector('[role=menu]')).toBeNull();
    expect(document.activeElement).toBe(opener);
    root.remove();
    opener.remove();
  });
});
