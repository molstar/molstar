/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { domView } from '../../../dom/view.js';
import { ShellModel } from '../../../shell/shell.js';
import { FloatingDock } from '../dock-view.js';

// A rail tab's menu (shell.md §3): Hide this tab, Show each hidden one, Reset layout; from the pointer's context menu or
// the keyboard's menu key.

beforeEach(() => {
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

const tab = (id: string, order: number) => ({
  id,
  title: id[0].toUpperCase() + id.slice(1),
  icon: 'square' as const,
  order,
  render: domView(() => null),
});

describe('a rail tab’s menu', () => {
  it('hides a tab, brings it back, and resets the layout', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('home', 1), tab('history', 2)]);
    const host = document.createElement('div');
    host.className = 'msp-root';
    document.body.append(host);
    render(() => <FloatingDock shell={shell} dock={shell.dock('right')} />, host);
    const button = (id: string) => host.querySelector<HTMLButtonElement>(`[role=tab][data-tab="${id}"]`);
    const rows = () => [...document.querySelectorAll('[role=menu] button')].map((b) => b.textContent);
    const choose = (label: string) =>
      [...document.querySelectorAll<HTMLButtonElement>('[role=menu] button')]
        .find((b) => b.textContent === label)!
        .click();

    button('history')!.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, clientX: 5, clientY: 5 }));
    expect(rows()).toEqual(['Hide History', 'Reset layout']);
    choose('Hide History');
    expect(button('history')).toBeNull();

    // From the keyboard: the menu key on a tab.
    button('home')!.dispatchEvent(new KeyboardEvent('keydown', { key: 'ContextMenu', bubbles: true }));
    expect(rows()).toEqual(['Hide Home', 'Show History', 'Reset layout']);
    choose('Show History');
    expect(button('history')).not.toBeNull();

    shell.hideTab('home');
    button('history')!.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true }));
    choose('Reset layout');
    expect(button('home')).not.toBeNull();
    expect(shell.hiddenTabs.value).toEqual([]);
  });
});
