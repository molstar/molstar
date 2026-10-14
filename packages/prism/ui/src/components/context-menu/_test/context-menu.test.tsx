// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { openContextMenu } from '../context-menu.js';

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

  it('opens with its items, runs one and closes on the choice', () => {
    const { root, run, closed } = setup();
    expect(menu(root)?.textContent).toContain('Color by pLDDT');
    (root.querySelector('[role="menuitem"]') as HTMLElement).click();
    expect(run).toHaveBeenCalledTimes(1);
    expect(menu(root)).toBeNull();
    expect(closed).toHaveBeenCalledTimes(1);
    root.remove();
  });

  it("shows an item's icon, and keeps the label column where an item has none", () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 0, y: 0 }, [
      { id: 'a', label: 'With', icon: 'gauge', run: () => {} },
      { id: 'b', label: 'Without', run: () => {} },
    ]);
    const rows = root.querySelectorAll('[role="menuitem"]');
    expect(rows[0].querySelector('svg')).not.toBeNull();
    expect(rows[1].querySelector('svg')).toBeNull();
    root.remove();
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
    a.root.remove();
    b.root.remove();
  });

  it('closes by its own return, once, and takes its listeners with it', () => {
    const { root, closed, close } = setup();
    close();
    close();
    expect(closed).toHaveBeenCalledTimes(1);
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
    expect(closed).toHaveBeenCalledTimes(1);
    root.remove();
  });
});
