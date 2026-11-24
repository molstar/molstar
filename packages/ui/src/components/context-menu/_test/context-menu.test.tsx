/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { openContextMenu, openMenuBeside } from '../context-menu.js';

afterEach(() => {
  // A menu left open keeps its listeners on the document past its DOM (the files share a worker): closed first.
  for (const menu of document.querySelectorAll('[role="menu"]')) {
    menu.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  }
  document.body.innerHTML = '';
});

describe('openContextMenu', () => {
  const setup = () => {
    const root = document.createElement('div');
    document.body.append(root);
    const run = vi.fn();
    const closed = vi.fn();
    const close = openContextMenu(root, { x: 10, y: 20 }, [{ id: 'a', label: 'Color by pLDDT', run }], closed);
    return { root, run, closed, close };
  };
  const menu = (root: HTMLElement) => root.querySelector('[role="menu"]');

  it('opens with its items, focuses the first, runs one and closes on the choice', () => {
    const { root, run, closed } = setup();
    expect(menu(root)?.textContent).toContain('Color by pLDDT');
    expect(document.activeElement).toBe(root.querySelector('[role="menuitem"]'));
    (root.querySelector('[role="menuitem"]') as HTMLElement).click();
    expect(run).toHaveBeenCalledTimes(1);
    expect(menu(root)).toBeNull();
    expect(closed).toHaveBeenCalledTimes(1);
  });

  it('keeps an icon slot for every item when one has an icon, and none when none has', () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 0, y: 0 }, [
      { id: 'a', label: 'With', icon: 'eye', run: () => {} },
      { id: 'b', label: 'Without', run: () => {} },
    ]);
    const rows = root.querySelectorAll('[role="menuitem"]');
    expect(rows[0].querySelector('svg')).not.toBeNull();
    expect(rows[1].querySelector('svg')).toBeNull();
    expect(rows[0].children).toHaveLength(3);
    expect(rows[1].children).toHaveLength(3);
    const plain = document.createElement('div');
    document.body.append(plain);
    openContextMenu(plain, { x: 0, y: 0 }, [{ id: 'a', label: 'Plain', run: () => {} }]);
    expect(plain.querySelector('[role="menuitem"]')!.children).toHaveLength(2);
  });

  it('greys a disabled item with its reason (its tooltip and description), opens past it, reaches it by the arrows, and does not run it', () => {
    const root = document.createElement('div');
    document.body.append(root);
    const run = vi.fn();
    openContextMenu(root, { x: 0, y: 0 }, [
      { id: 'a', label: 'Off', disabled: 'Nothing selected', run },
      { id: 'b', label: 'On', run: () => {} },
    ]);
    const [off, on] = root.querySelectorAll<HTMLElement>('[role="menuitem"]');
    expect(off.getAttribute('aria-disabled')).toBe('true');
    expect(off.dataset.tip).toBe('Nothing selected');
    expect(off.getAttribute('aria-description')).toBe('Nothing selected');
    expect(on.hasAttribute('aria-description')).toBe(false);
    expect(document.activeElement).toBe(on);
    // APG menus: a disabled item takes the focus (its reason is read), it is not activated.
    on.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowUp', bubbles: true, cancelable: true }));
    expect(document.activeElement).toBe(off);
    off.click();
    expect(run).not.toHaveBeenCalled();
    expect(menu(root)).not.toBeNull();
  });

  it('closes on Escape and on a press elsewhere, and not on a press inside', () => {
    const a = setup();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(menu(a.root)).toBeNull();
    const b = setup();
    menu(b.root)?.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(menu(b.root)).not.toBeNull();
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(menu(b.root)).toBeNull();
  });

  it("does not take Escape pressed in another UI root's shadow root, and does not stop it", () => {
    const hostA = document.createElement('div');
    const hostB = document.createElement('div');
    document.body.append(hostA, hostB);
    const rootA = document.createElement('div');
    hostA.attachShadow({ mode: 'open' }).append(rootA);
    const inputB = document.createElement('input');
    hostB.attachShadow({ mode: 'open' }).append(inputB);
    openContextMenu(rootA, { x: 0, y: 0 }, [{ id: 'a', label: 'A', run: () => {} }]);
    const seen = vi.fn();
    document.addEventListener('keydown', seen);
    inputB.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, composed: true }));
    expect(menu(rootA)).not.toBeNull();
    expect(seen).toHaveBeenCalledTimes(1);
    // inside A it closes
    rootA
      .querySelector('button')
      ?.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, composed: true }));
    expect(menu(rootA)).toBeNull();
    document.removeEventListener('keydown', seen);
  });

  it('its rows are no Tab stops; Tab closes it, the focus back where it was', () => {
    const before = document.createElement('button');
    document.body.append(before);
    before.focus();
    const { root, closed } = setup();
    const row = root.querySelector<HTMLElement>('[role="menuitem"]')!;
    expect(row.tabIndex).toBe(-1);
    expect(document.activeElement).toBe(row);
    row.dispatchEvent(new KeyboardEvent('keydown', { key: 'Tab', bubbles: true, cancelable: true }));
    expect(menu(root)).toBeNull();
    expect(closed).toHaveBeenCalledTimes(1);
    expect(document.activeElement).toBe(before);
  });

  it('closes by its own return, once, and takes its listeners with it', () => {
    const { closed, close } = setup();
    close();
    close();
    expect(closed).toHaveBeenCalledTimes(1);
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
    expect(closed).toHaveBeenCalledTimes(1);
  });
});

describe('openMenuBeside', () => {
  /** A button at a place of the viewport (jsdom has no layout: its rect is given). */
  const button = (rect: { left: number; top: number; width: number; height: number }) => {
    const el = document.createElement('button');
    document.body.append(el);
    el.getBoundingClientRect = () =>
      ({
        ...rect,
        right: rect.left + rect.width,
        bottom: rect.top + rect.height,
        x: rect.left,
        y: rect.top,
      }) as DOMRect;
    return el;
  };
  const choice = (chosen: string) =>
    ['light', 'dark', 'auto'].map((id) => ({ id, label: id, checked: id === chosen, run: vi.fn() }));

  it('a choice: radio rows, the chosen one checked and focused', () => {
    const root = document.createElement('div');
    document.body.append(root);
    const close = openMenuBeside(root, button({ left: 900, top: 100, width: 32, height: 32 }), choice('dark'));
    const rows = [...root.querySelectorAll('[role="menuitemradio"]')];
    expect(rows.map((r) => r.getAttribute('aria-checked'))).toEqual(['false', 'true', 'false']);
    expect(document.activeElement).toBe(rows[1]);
    // The check's slot is kept on every row: the labels hold their column.
    expect(rows.every((r) => r.lastElementChild?.tagName === 'SPAN')).toBe(true);
    // Closed: an open menu's listeners outlive the DOM (the files share a worker).
    close();
  });

  it('opens to the left of the button, its top at the button’s; low on the window, its bottom at the button’s', () => {
    const root = document.createElement('div');
    document.body.append(root);
    const place = (top: number) => {
      root.innerHTML = '';
      const close = openMenuBeside(root, button({ left: 900, top, width: 32, height: 32 }), choice('light'));
      const transform = root.querySelector<HTMLElement>('[role="menu"]')!.style.transform;
      close();
      return transform;
    };
    // jsdom's own getters, put back after (the files share a worker's modules and globals).
    const width = Object.getOwnPropertyDescriptor(HTMLElement.prototype, 'offsetWidth');
    const height = Object.getOwnPropertyDescriptor(HTMLElement.prototype, 'offsetHeight');
    Object.defineProperty(HTMLElement.prototype, 'offsetWidth', { configurable: true, get: () => 160 });
    Object.defineProperty(HTMLElement.prototype, 'offsetHeight', { configurable: true, get: () => 100 });
    try {
      // window.innerHeight is 768 in jsdom: the right edge 4 px left of the button.
      expect(place(100)).toBe('translate(736px, 100px)');
      expect(place(700)).toBe(`translate(736px, ${700 + 32 - 100}px)`);
    } finally {
      if (width) Object.defineProperty(HTMLElement.prototype, 'offsetWidth', width);
      if (height) Object.defineProperty(HTMLElement.prototype, 'offsetHeight', height);
    }
  });
});

describe('a submenu (a flyout)', () => {
  const setup = () => {
    const root = document.createElement('div');
    document.body.append(root);
    const pick = vi.fn();
    const closed = vi.fn();
    const subClosed = vi.fn();
    const close = openContextMenu(
      root,
      { x: 10, y: 20 },
      [
        { id: 'hide', label: 'Hide', run: vi.fn() },
        {
          id: 'color',
          label: 'Color',
          run: vi.fn(),
          items: [
            { id: 'a', label: 'By pLDDT', run: pick },
            { id: 'b', label: 'By chain', run: vi.fn() },
          ],
        },
      ],
      closed,
      { onSubmenuClose: subClosed },
    );
    const menus = () => [...root.querySelectorAll<HTMLElement>('[role="menu"]')];
    const color = root.querySelector<HTMLElement>('[aria-haspopup="menu"]')!;
    return { root, pick, closed, subClosed, close, menus, color };
  };
  const key = (el: Element, k: string) =>
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));

  it('a row with a submenu ends in its arrow (the label never carries one), the same last cell as a key or a check', () => {
    const { color } = setup();
    expect(color.textContent).toBe('Color');
    const cells = [...color.children];
    const last = cells[cells.length - 1];
    expect(last.querySelector('svg')).not.toBeNull();
    expect(last.hasAttribute('data-submenu-arrow')).toBe(true);
    // the label is the row's first text cell: the arrow comes after it, at the row's end
    expect(cells.indexOf(color.querySelector('span:not([data-submenu-arrow])')!)).toBeLessThan(cells.indexOf(last));
  });

  it('a row with a key and a submenu shows both: the key, then the arrow, in one trailing cell', () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 0, y: 0 }, [
      { id: 'a', label: 'Both', shortcut: 'K', run: () => {}, items: [{ id: 'b', label: 'Sub', run: () => {} }] },
    ]);
    const row = root.querySelector<HTMLElement>('[aria-haspopup="menu"]')!;
    expect(row.children).toHaveLength(2);
    const end = row.lastElementChild!;
    expect(end.querySelector('kbd')?.textContent).toBe('K');
    const arrow = end.querySelector('[data-submenu-arrow]')!;
    expect(arrow.querySelector('svg')).not.toBeNull();
    expect(end.lastElementChild).toBe(arrow);
  });

  it('opens beside its row on a click, the menu staying, the row marked; a choice runs it and closes both', () => {
    const { pick, closed, menus, color } = setup();
    expect(color.getAttribute('aria-expanded')).toBe('false');
    color.click();
    expect(menus()).toHaveLength(2);
    expect(menus()[0].textContent).toContain('Hide');
    expect(color.getAttribute('aria-expanded')).toBe('true');
    (menus()[1].querySelector('[role="menuitem"]') as HTMLElement).click();
    expect(pick).toHaveBeenCalledTimes(1);
    expect(menus()).toHaveLength(0);
    expect(closed).toHaveBeenCalledTimes(1);
  });

  it('→ opens it with the focus in; ← and Escape close it alone, the focus back on the row', () => {
    const { menus, color, subClosed, closed } = setup();
    color.focus();
    key(color, 'ArrowRight');
    expect(menus()[1].contains(document.activeElement)).toBe(true);
    key(document.activeElement!, 'ArrowLeft');
    expect(menus()).toHaveLength(1);
    expect(document.activeElement).toBe(color);
    key(color, 'ArrowRight');
    key(document.activeElement!, 'Escape');
    expect(menus()).toHaveLength(1);
    expect(document.activeElement).toBe(color);
    expect(subClosed).toHaveBeenCalledTimes(2);
    expect(closed).not.toHaveBeenCalled();
    key(color, 'Escape');
    expect(menus()).toHaveLength(0);
    expect(closed).toHaveBeenCalledTimes(1);
  });

  it('a press in the flyout is not outside; a press outside closes all', () => {
    const { menus, color, closed } = setup();
    color.click();
    menus()[1].dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(menus()).toHaveLength(2);
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(menus()).toHaveLength(0);
    expect(closed).toHaveBeenCalledTimes(1);
  });

  it('Tab in a flyout closes the whole menu, the focus back where it was', () => {
    const before = document.createElement('button');
    document.body.append(before);
    before.focus();
    const { menus, color, closed } = setup();
    color.focus();
    key(color, 'ArrowRight');
    expect(menus()[1].contains(document.activeElement)).toBe(true);
    key(document.activeElement!, 'Tab');
    expect(menus()).toHaveLength(0);
    expect(closed).toHaveBeenCalledTimes(1);
    expect(document.activeElement).toBe(before);
  });

  it('a choice in a flyout puts the focus back too', () => {
    const before = document.createElement('button');
    document.body.append(before);
    before.focus();
    const { menus, color } = setup();
    color.click();
    (menus()[1].querySelector('[role="menuitem"]') as HTMLElement).focus();
    (menus()[1].querySelector('[role="menuitem"]') as HTMLElement).click();
    expect(document.activeElement).toBe(before);
  });

  it('names its pair for assistive technology: the row controls the flyout, the flyout is labelled by the row; it scrolls inside the window', () => {
    const { menus, color } = setup();
    color.click();
    const flyout = menus()[1];
    expect(color.getAttribute('aria-controls')).toBe(flyout.id);
    expect(flyout.getAttribute('aria-labelledby')).toBe(color.id);
    expect(color.id).not.toBe('');
    expect(flyout.style.overflowY).toBe('auto');
    expect(flyout.style.maxHeight).toBe(`${window.innerHeight - 8}px`);
  });

  it('leaving the menu stops a pending rest: nothing opens after the pointer is gone', () => {
    vi.useFakeTimers();
    try {
      const { menus, color } = setup();
      color.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
      menus()[0].dispatchEvent(new PointerEvent('pointerleave'));
      vi.advanceTimersByTime(400);
      expect(menus()).toHaveLength(1);
    } finally {
      vi.useRealTimers();
    }
  });

  it('a menu that opened under a still pointer: its first move on the row (no pointerover came) starts the rest', () => {
    vi.useFakeTimers();
    try {
      const { menus, color } = setup();
      color.dispatchEvent(new PointerEvent('pointermove', { bubbles: true }));
      color.dispatchEvent(new PointerEvent('pointermove', { bubbles: true }));
      vi.advanceTimersByTime(300);
      expect(menus()).toHaveLength(2);
    } finally {
      vi.useRealTimers();
    }
  });

  it('draws a separator element above an item that starts a group (`separated`), never one above the first; the arrows skip it', () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 0, y: 0 }, [
      { id: 'a', label: 'A', separated: true, run: () => {} },
      { id: 'b', label: 'B', run: () => {} },
      { id: 'c', label: 'C', separated: true, run: () => {} },
    ]);
    const menu = root.querySelector('[role="menu"]')!;
    expect([...menu.children].map((c) => c.getAttribute('role'))).toEqual([
      'menuitem',
      'menuitem',
      'separator',
      'menuitem',
    ]);
    const rows = [...menu.querySelectorAll<HTMLElement>('[role="menuitem"]')];
    rows[1].focus();
    rows[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
    expect(document.activeElement).toBe(rows[2]);
  });

  it.each(['transform: translate(10px, 10px)', 'filter: blur(0px)', 'contain: paint', 'will-change: transform'])(
    'inside an ancestor that would offset its fixed position (%s) it is put in the body, in a wrapper with the UI root’s class and theme',
    (style) => {
      const ui = document.createElement('div');
      ui.className = 'msui-root';
      ui.setAttribute('data-msui-theme', 'dark');
      ui.style.setProperty('--msui-accent', 'red');
      const skew = document.createElement('div');
      skew.setAttribute('style', style);
      const root = document.createElement('div');
      skew.append(root);
      ui.append(skew);
      document.body.append(ui);
      const close = openContextMenu(root, { x: 0, y: 0 }, [
        { id: 'a', label: 'A', run: () => {} },
        { id: 'b', label: 'B', run: () => {} },
      ]);
      const menu = document.querySelector<HTMLElement>('[role="menu"]')!;
      expect(root.contains(menu)).toBe(false);
      const wrap = menu.parentElement!;
      expect(wrap.parentElement).toBe(document.body);
      expect(wrap.classList.contains('msui-root')).toBe(true);
      expect(wrap.getAttribute('data-msui-theme')).toBe('dark');
      expect(wrap.style.getPropertyValue('--msui-accent')).toBe('red');
      // no ancestor of the menu is the one that offsets it
      for (let el: Element | null = menu; el; el = el.parentElement) expect(el.getAttribute('style')).not.toBe(style);
      close();
      expect(wrap.isConnected).toBe(false);
    },
  );

  it('stays in its root where no ancestor offsets it', () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 0, y: 0 }, [{ id: 'a', label: 'A', run: () => {} }]);
    expect(root.querySelector('[role="menu"]')).not.toBeNull();
  });

  it('opens when the pointer rests on the row', () => {
    vi.useFakeTimers();
    try {
      const { menus, color } = setup();
      color.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
      expect(menus()).toHaveLength(1);
      vi.advanceTimersByTime(300);
      expect(menus()).toHaveLength(2);
    } finally {
      vi.useRealTimers();
    }
  });
});
