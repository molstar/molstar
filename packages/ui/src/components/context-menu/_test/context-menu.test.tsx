/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { openContextMenu } from '../context-menu.js';

afterEach(() => {
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

  it('greys a disabled item with its reason, skips it for the focus, and does not run it', () => {
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
    expect(document.activeElement).toBe(on);
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

  it('closes by its own return, once, and takes its listeners with it', () => {
    const { closed, close } = setup();
    close();
    close();
    expect(closed).toHaveBeenCalledTimes(1);
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
    expect(closed).toHaveBeenCalledTimes(1);
  });
});
