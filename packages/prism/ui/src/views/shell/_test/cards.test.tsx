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
import { flush, onCleanup } from '../../../dom/scope.js';
import { domView } from '../../../dom/view.js';
import { render } from '../../../dom/render.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { ShellModel, type ShellOptions, type ShellViewContext } from '../../../shell/shell.js';
import { ShellView } from '../shell-view.js';

// Where a card opens (shell.md, Cards and Compact): over the canvas in the panel's place, or, compact, in the bottom
// sheet in the panel's place; one card element, moved between the two as the layout changes.

let mounts = 0;
let unmounts = 0;

function mount(options?: ShellOptions) {
  const shell = new ShellModel(options);
  shell.registerTabs([
    { id: 'structure', title: 'Structure', icon: 'house', render: domView<ShellViewContext>(() => <i>panel</i>) },
  ]);
  shell.registerOverlays([
    {
      id: 'card',
      title: 'Card',
      render: domView<ShellViewContext>(() => {
        mounts++;
        onCleanup(() => unmounts++);
        return <b>card</b>;
      }),
    },
  ]);
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ShellView shell={shell} status={new NotificationsManager()} canvas={<canvas />} />, container);
  const card = () => container.querySelector<HTMLElement>('[role="dialog"]');
  const sheet = () => container.querySelector<HTMLElement>('main ~ [data-covered], main ~ div:not([data-covered])');
  return { shell, container, card, sheet };
}

describe('overlay cards', () => {
  beforeEach(() => {
    mounts = unmounts = 0;
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
    document.body.replaceChildren();
    vi.unstubAllGlobals();
  });

  it('compact: opens in the sheet in the panel’s place, the panel kept out of it, and back when the card closes', () => {
    const { shell, container, card } = mount({ mode: 'compact' });
    shell.show('structure');
    shell.openOverlay('card');
    flush();
    const main = container.querySelector('main') as HTMLElement;
    const dock = card()?.closest('[data-covered]') as HTMLElement;
    expect(dock).not.toBeNull();
    expect(main.contains(card())).toBe(false);
    // Not a dialog over the layout: part of it, as the panel is.
    expect(card()?.getAttribute('aria-modal')).toBe('false');
    const panel = dock.querySelector('[role="tabpanel"]') as HTMLElement;
    expect(panel.hasAttribute('inert')).toBe(true);

    shell.closeOverlay();
    flush();
    expect(card()).toBeNull();
    expect(dock.hasAttribute('data-covered')).toBe(false);
    expect(panel.hasAttribute('inert')).toBe(false);
  });

  it('moves the one card between the canvas and the sheet as the mode changes, never making a second', () => {
    const { shell, container, card } = mount({ mode: 'desktop', docks: { right: { pinned: false } } });
    shell.openOverlay('card');
    flush();
    const main = container.querySelector('main') as HTMLElement;
    const first = card();
    expect(main.contains(first)).toBe(true);
    expect(first?.getAttribute('aria-modal')).toBe('true');

    shell.setMode('compact');
    flush();
    expect(card()).toBe(first);
    expect(main.contains(first)).toBe(false);
    expect(first?.getAttribute('aria-modal')).toBe('false');

    shell.setMode('desktop');
    flush();
    expect(card()).toBe(first);
    expect(main.contains(first)).toBe(true);
    // The card's view was made once, and is still there.
    expect([mounts, unmounts]).toEqual([1, 0]);
    expect(container.querySelectorAll('[role="dialog"]')).toHaveLength(1);
  });

  it('compact and expanded (no sheet): stays over the canvas, where it can be seen', () => {
    const { shell, container, card } = mount({ mode: 'compact' });
    shell.setExpanded(true);
    flush();
    shell.openOverlay('card');
    flush();
    const main = container.querySelector('main') as HTMLElement;
    expect(card()).not.toBeNull();
    expect(main.contains(card())).toBe(true);
  });

  it('desktop, docked: the panel beside the card stays live (not covered)', () => {
    const { shell, container, card } = mount({ mode: 'desktop', docks: { right: { floating: false } } });
    shell.show('structure');
    shell.openOverlay('card');
    flush();
    expect(card()).not.toBeNull();
    const panel = container.querySelector('[role="tabpanel"]') as HTMLElement;
    expect(panel.hasAttribute('inert')).toBe(false);
  });

  it('offers no pin beside a compact rail, and one with the full rail', () => {
    const { shell, card } = mount({ mode: 'desktop', compactRail: true });
    shell.openOverlay('card');
    flush();
    expect(card()?.querySelector('button[aria-label="Pin beside the canvas"]')).toBeNull();
    shell.setCompactRail(false);
    flush();
    shell.openOverlay('card');
    flush();
    expect(card()?.querySelector('button[aria-label="Pin beside the canvas"]')).not.toBeNull();
  });
});
