/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, of } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { type CanvasMenuItem, ShellModel } from '../../../shell/shell.js';
import { CanvasStack, ToolMenuRows } from '../canvas-chrome.js';

// A menu row with `items` is a submenu: beside the canvas stack's menu (a click, →, hover; ← back to its row), its rows
// following the list while open; under the row where a menu has no room beside it.

beforeEach(() => {
  vi.stubGlobal('matchMedia', (q: string) => ({ matches: false, media: q }) as MediaQueryList);
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      disconnect() {}
      unobserve() {}
    },
  );
});
afterEach(() => {
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

function rowsOf(list: string[], ran: string[]) {
  return list.map((id): CanvasMenuItem => ({ id, label: id, pressed: of(id !== 'b'), run: () => void ran.push(id) }));
}

describe('menu submenus', () => {
  it('opens beside the canvas stack menu, follows its list, runs a row and goes back with ←', () => {
    const shell = new ShellModel();
    const ran: string[] = [];
    const list = new BehaviorSubject(rowsOf(['a', 'b'], ran));
    shell.registerToolbar([
      {
        id: 'seq',
        place: 'canvas',
        icon: 'sequence',
        title: 'Sequence',
        menu: [
          { id: 'one', label: 'One', run: () => ran.push('one') },
          { id: 'sub', label: 'Tracks', run: () => ran.push('never'), items: list },
        ],
      },
    ]);
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <CanvasStack shell={shell} />, host);
    host
      .querySelector<HTMLButtonElement>('button[aria-label="Sequence"]')!
      .dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
    const menus = () => [...host.querySelectorAll<HTMLElement>('[role=menu]')];
    expect(menus()).toHaveLength(1);
    const sub = menus()[0].querySelector<HTMLButtonElement>('[aria-haspopup=menu]')!;
    expect(sub.textContent).toBe('Tracks');
    sub.click();
    expect(ran).toEqual([]);
    expect(menus()).toHaveLength(2);
    const subRows = () => [...menus()[1].querySelectorAll('button')].map((b) => b.textContent);
    expect(subRows()).toEqual(['a', 'b']);
    expect([...menus()[1].querySelectorAll('button')].map((b) => b.getAttribute('aria-checked'))).toEqual([
      'true',
      'false',
    ]);
    list.next(rowsOf(['a', 'b', 'c'], ran));
    expect(subRows()).toEqual(['a', 'b', 'c']);
    // ← in the submenu: back to its row, the menu still open.
    menus()[1]
      .querySelector('button')!
      .dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowLeft', bubbles: true }));
    expect(menus()).toHaveLength(1);
    // → on the row opens it again; a row runs and the menu closes.
    sub.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    expect(menus()).toHaveLength(2);
    (menus()[1].querySelectorAll('button')[2] as HTMLButtonElement).click();
    expect(ran).toEqual(['c']);
    expect(menus()).toHaveLength(0);
  });

  it('shows its rows under it where the menu has no room beside it', () => {
    const ran: string[] = [];
    const host = document.createElement('div');
    document.body.append(host);
    render(
      () => (
        <ToolMenuRows
          items={[{ id: 'sub', label: 'Tracks', run: () => {}, items: of(rowsOf(['a'], ran)) }]}
          onRun={() => {}}
        />
      ),
      host,
    );
    expect(host.textContent).toBe('Tracksa');
    host.querySelector('button')!.click();
    expect(ran).toEqual(['a']);
  });
});
