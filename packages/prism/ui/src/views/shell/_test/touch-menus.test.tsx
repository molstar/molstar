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
import { ShellModel } from '../../../shell/shell.js';
import { CanvasStack } from '../canvas-chrome.js';

// A touch screen has no hover: a canvas button with an action and a menu (the camera) opens its menu on a tap, its own
// action the first row; a tap outside closes it.

beforeEach(() => {
  vi.stubGlobal('matchMedia', (q: string) => ({ matches: q === '(hover: none)', media: q }) as MediaQueryList);
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

describe('canvas menus on a touch screen', () => {
  it('opens on a tap with the button’s own action first, and closes on a tap outside', () => {
    const shell = new ShellModel();
    const ran: string[] = [];
    shell.registerToolbar([
      {
        id: 'reset',
        place: 'canvas',
        icon: 'focus',
        title: 'Reset camera',
        run: () => ran.push('reset'),
        menu: [{ id: 'zoom', label: 'Reset zoom', run: () => ran.push('zoom') }],
      },
    ]);
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <CanvasStack shell={shell} />, host);
    const button = host.querySelector<HTMLButtonElement>('button[aria-label="Reset camera"]')!;
    button.click();
    expect(ran).toEqual([]);
    const rows = () => [...host.querySelectorAll('[role=menu] button')].map((b) => b.textContent);
    expect(rows()).toEqual(['Reset camera', 'Reset zoom']);
    (host.querySelector('[role=menu] button') as HTMLButtonElement).click();
    expect(ran).toEqual(['reset']);
    button.click();
    expect(rows()).toHaveLength(2);
    document.body.dispatchEvent(new Event('pointerdown', { bubbles: true, composed: true }));
    expect(rows()).toHaveLength(0);
  });
});
