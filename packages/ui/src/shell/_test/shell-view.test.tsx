/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { flush } from '../../dom/index.js';
import { stubPointerCapture } from '../../testing/gestures.js';
import { findView, railGroup } from '../layout.js';
import { key, microtasks, mountShell, probeView } from './fixture.js';

beforeAll(() => stubPointerCapture());

afterEach(() => {
  document.body.innerHTML = '';
});

const menu = () => document.querySelector('[role="menu"]');
const menuItem = (label: string) =>
  [...document.querySelectorAll<HTMLElement>('[role="menuitem"]')].find((m) => m.textContent?.includes(label));

describe('ShellView: the rail', () => {
  const views = () => [
    probeView('home', { place: 'head', order: 0 }),
    probeView('layers', { order: 1 }),
    probeView('inspect', { order: 2 }),
    probeView('help', { place: 'foot', order: 3 }),
  ];

  it('shows the rail group’s views as tabs in head, strip, foot order', () => {
    const m = mountShell(views());
    const tabs = m.qa('nav [role="tab"]').map((b) => b.getAttribute('aria-label'));
    expect(tabs).toEqual(['Home', 'Layers', 'Inspect', 'Help']);
    const lists = m.qa('nav [role="tablist"]').map((l) => l.getAttribute('aria-label'));
    expect(lists).toEqual(['head', 'strip', 'foot']);
    m.dispose();
  });

  it('a press opens the group, a press on the open view folds it (the column hides, the rail stays)', () => {
    const m = mountShell(views());
    // Registered views are placed, not opened.
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    expect(m.railButton('Layers')?.getAttribute('aria-selected')).toBe('false');
    m.railButton('Layers')!.click();
    flush();
    expect(m.groupNode('right')?.hasAttribute('data-shown')).toBe(false);
    expect(m.frame('layers')).not.toBeNull();
    expect(m.railButton('Layers')?.getAttribute('aria-selected')).toBe('true');
    expect(m.railButton('Home')?.getAttribute('aria-selected')).toBe('false');
    m.railButton('Inspect')!.click();
    flush();
    expect(m.frame('inspect')).not.toBeNull();
    expect(m.frame('layers')).toBeNull();
    expect(m.railButton('Inspect')?.getAttribute('aria-selected')).toBe('true');
    m.railButton('Inspect')!.click();
    flush();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    expect(m.rail()).not.toBeNull();
    expect(m.railButton('Inspect')?.getAttribute('aria-selected')).toBe('false');
    m.dispose();
  });

  it('Shift+F10 and the context menu open the rail menu: Hide, Show, Reset layout', () => {
    const m = mountShell(views());
    key(m.railButton('Layers')!, 'F10', { shiftKey: true });
    expect(menu()).not.toBeNull();
    menuItem('Hide')!.click();
    flush();
    expect(m.railButton('Layers')).toBeNull();
    expect(m.shell.layout.value.hidden).toEqual(['layers']);
    m.railButton('Home')!.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true }));
    flush();
    expect(menuItem('Show Layers')).toBeDefined();
    menuItem('Show Layers')!.click();
    flush();
    expect(m.railButton('Layers')).not.toBeNull();
    expect(m.frame('layers')).not.toBeNull();
    m.shell.hide('inspect');
    flush();
    key(m.railButton('Home')!, 'ContextMenu');
    menuItem('Reset layout')!.click();
    flush();
    expect(m.railButton('Inspect')).not.toBeNull();
    expect(m.shell.layout.value.hidden).toEqual([]);
    m.dispose();
  });

  it('runs rail items, shows their pressed state, and puts the overflow under More', () => {
    const run = vi.fn();
    const pressed = new State(false);
    const m = mountShell(views());
    m.shell.registerRailItems([{ id: 'theme', icon: 'sun', title: 'Theme', pressed, run }]);
    flush();
    const theme = m.railButton('Theme')!;
    theme.click();
    expect(run).toHaveBeenCalledWith(theme);
    expect(theme.getAttribute('aria-pressed')).toBe('false');
    pressed.set(true);
    flush();
    expect(theme.getAttribute('aria-pressed')).toBe('true');
    m.dispose();
  });

  it('a rail item with a menu opens it on a press (not on hover), read as it opens; a choice picks', () => {
    const theme = new State('light');
    const m = mountShell(views());
    m.shell.registerRailItems([
      {
        id: 'theme',
        icon: 'sun-moon',
        title: 'Theme',
        menu: () =>
          ['light', 'dark'].map((id) => ({ id, label: id, checked: theme.value === id, run: () => theme.set(id) })),
      },
    ]);
    flush();
    const button = m.railButton('Theme')!;
    expect(button.getAttribute('aria-haspopup')).toBe('menu');
    button.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(document.querySelector('[role="menu"]')).toBeNull();
    button.click();
    const rows = [...document.querySelectorAll<HTMLElement>('[role="menuitemradio"]')];
    expect(rows.map((r) => r.getAttribute('aria-checked'))).toEqual(['true', 'false']);
    expect(button.getAttribute('aria-expanded')).toBe('true');
    rows[1].click();
    expect(theme.value).toBe('dark');
    expect(button.getAttribute('aria-expanded')).toBe('false');
    // Opened again: it reads the choice anew.
    button.click();
    expect(document.querySelector('[aria-checked="true"]')?.textContent).toContain('dark');
    // Closed before the shell goes: an open menu's listeners would outlive the test (files share a worker).
    key(document.activeElement!, 'Escape');
    expect(document.querySelector('[role="menu"]')).toBeNull();
    m.dispose();
  });

  it('a visible flag removes and restores a view’s rail button', () => {
    const shown = new State(false);
    const m = mountShell([...views(), probeView('kin', { visible: shown, order: 1.5 })]);
    expect(m.railButton('Kin')).toBeNull();
    shown.set(true);
    flush();
    const tabs = m.qa('nav [role="tab"]').map((b) => b.getAttribute('aria-label'));
    expect(tabs).toEqual(['Home', 'Layers', 'Kin', 'Inspect', 'Help']);
    m.railButton('Kin')!.click();
    flush();
    expect(m.frame('kin')).not.toBeNull();
    shown.set(false);
    flush();
    expect(m.railButton('Kin')).toBeNull();
    expect(m.frame('kin')).toBeNull();
    m.dispose();
  });
});

describe('ShellView: views keep their state', () => {
  it('an input’s value survives float and dock back; the body is paused while not shown; leaving disposes it', async () => {
    const a = probeView('layers');
    const b = probeView('inspect');
    const m = mountShell([a, b]);
    m.shell.show('layers');
    flush();
    const input = m.q<HTMLInputElement>('[data-input="layers"]')!;
    input.value = 'typed';
    expect(a.probe.scope?.active).toBe(true);
    m.shell.float('layers');
    flush();
    expect(m.win('layers')).not.toBeNull();
    expect(m.q<HTMLInputElement>('[data-input="layers"]')).toBe(input);
    expect(input.value).toBe('typed');
    m.shell.dock('layers');
    flush();
    expect(m.win('layers')).toBeNull();
    expect(m.q<HTMLInputElement>('[data-input="layers"]')).toBe(input);
    expect(a.probe.mounts).toBe(1);
    // Another view in its group: kept, paused, out of the document.
    m.shell.show('inspect');
    flush();
    expect(input.isConnected).toBe(false);
    expect(a.probe.scope?.active).toBe(false);
    m.shell.show('layers');
    flush();
    expect(input.isConnected).toBe(true);
    expect(a.probe.scope?.active).toBe(true);
    expect(input.value).toBe('typed');
    // Folded: paused.
    m.shell.select('layers');
    flush();
    expect(a.probe.scope?.active).toBe(false);
    expect(a.probe.cleanups).toBe(0);
    m.dispose();
    expect(a.probe.cleanups).toBe(1);
    await microtasks();
  });

  it('a view leaving the layout (an instance closed) is disposed', () => {
    const plot = probeView('plot', { instances: 'many', home: 'float' });
    const m = mountShell([plot]);
    m.shell.show('plot', { key: 'x' });
    flush();
    expect(m.win('plot#x')).not.toBeNull();
    expect(plot.probe.mounts).toBe(1);
    m.shell.close('plot#x');
    flush();
    expect(m.win('plot#x')).toBeNull();
    expect(plot.probe.cleanups).toBe(1);
    expect(findView(m.shell.layout.value, 'plot#x')).toBeUndefined();
    m.dispose();
  });

  it('an unregistered view’s body goes', () => {
    const a = probeView('layers');
    const extra = probeView('extra');
    const m = mountShell([a]);
    const undo = m.shell.registerViews([extra.def]);
    m.shell.show('extra');
    flush();
    expect(extra.probe.mounts).toBe(1);
    undo();
    flush();
    expect(m.frame('extra')).toBeNull();
    expect(extra.probe.cleanups).toBe(1);
    expect(railGroup(m.shell.layout.value)?.views).toContain('extra');
    m.dispose();
  });
});

describe('ShellView: the defaults (no ⋯, ✕ closes, the rail hides)', () => {
  const views = () => [probeView('home', { place: 'head', order: 0 }), probeView('layers', { order: 1 })];

  it('a docked view has ✕ (it folds the group) and no ⋯; its header’s right-click is the layout menu', () => {
    const m = mountShell(views(), {}, { viewMenu: false });
    m.shell.show('layers');
    flush();
    const frame = m.frame('layers')!;
    expect(frame.querySelector('[aria-haspopup="menu"]')).toBeNull();
    frame.querySelector('header')!.dispatchEvent(new MouseEvent('contextmenu', { bubbles: true, cancelable: true }));
    expect(menuItem('Float')).toBeDefined();
    document
      .querySelector<HTMLElement>('[role="menu"]')
      ?.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    flush();
    const close = [...frame.querySelectorAll('button')].find((b) => b.getAttribute('aria-label') === 'Close Layers');
    close!.click();
    flush();
    expect(railGroup(m.shell.layout.value)?.collapsed).toBe(true);
    m.dispose();
  });

  it('the rail’s Hide hides it and the docked column; the star shows them again', () => {
    const m = mountShell(views(), {}, { viewMenu: false });
    m.shell.show('home');
    flush();
    m.q<HTMLButtonElement>('nav button[aria-label="Hide the views"]')!.click();
    flush();
    expect(m.shell.railShown.value).toBe(false);
    expect(m.rail()).toBeNull();
    m.q<HTMLButtonElement>('[data-msui-shell-main] button')!.click();
    flush();
    expect(m.shell.railShown.value).toBe(true);
    m.dispose();
  });

  it('without `railToggle` the rail has no Hide', () => {
    const m = mountShell(views(), {}, { railToggle: false });
    expect(m.q('nav button[aria-label="Hide the views"]')).toBeNull();
    m.dispose();
  });
});

describe('ShellView: where the focus goes when a layer opens', () => {
  it('a modal focuses the first control of its content, not its header’s ✕', async () => {
    const m = mountShell([probeView('home'), probeView('confirm', { home: 'modal' })], {}, { viewMenu: false });
    m.shell.show('confirm');
    flush();
    await microtasks();
    expect(document.activeElement?.getAttribute('data-input')).toBe('confirm');
    m.dispose();
  });

  it('a card opened by a press (not the keyboard) focuses itself, not a control', async () => {
    const m = mountShell([probeView('home'), probeView('export', { home: 'card' })]);
    m.shell.show('export');
    flush();
    await microtasks();
    expect(document.activeElement?.hasAttribute('data-card')).toBe(true);
    m.dispose();
  });
});

describe('ShellView: a split column never strands its views', () => {
  const views = () => [
    probeView('home', { place: 'head', order: 0 }),
    probeView('layers', { order: 1 }),
    probeView('inspect', { order: 2 }),
    probeView('settings', { order: 3 }),
  ];
  const closeOf = (frame: Element, title: string) =>
    [...frame.querySelectorAll('button')].find((b) => b.getAttribute('aria-label') === `Close ${title}`)!;

  it('✕ on the lower group returns its views to the rail’s group, which keeps what it showed', () => {
    const m = mountShell(views(), {}, { viewMenu: false });
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    m.shell.move('settings', { kind: 'group', id: findView(m.shell.layout.value, 'inspect')!.group.id });
    flush();
    closeOf(m.frame('settings')!, 'Settings').click();
    flush();
    const rail = railGroup(m.shell.layout.value)!;
    expect(rail.views).toEqual(['home', 'layers', 'inspect', 'settings']);
    expect(rail.active).toBe('layers');
    expect(m.frame('layers')).not.toBeNull();
    m.dispose();
  });

  it('with the rail’s group folded, a folded lower group does not show alone: the column closes to the rail', () => {
    const m = mountShell(views(), {}, { viewMenu: false });
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    flush();
    const lower = findView(m.shell.layout.value, 'inspect')!.group.id;
    m.shell.collapse(lower);
    m.shell.collapse('right');
    flush();
    expect(m.groupNode(lower)?.getAttribute('data-shown')).toBe('false');
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    expect(m.rail()).not.toBeNull();
    // A rail button opens its view again.
    m.railButton('Layers')!.click();
    flush();
    expect(m.frame('layers')).not.toBeNull();
    m.dispose();
  });
});

describe('ShellView: a view’s footer', () => {
  it('moves under the body’s scroll, wherever the view writes it', async () => {
    const { ViewFooter, ViewStack } = await import('../view-frame.js');
    const { domView } = await import('../../dom/index.js');
    const m = mountShell([
      probeView('form', {
        render: domView(() => (
          <ViewStack>
            <p>content</p>
            <ViewFooter>
              <button type="button">Apply</button>
            </ViewFooter>
          </ViewStack>
        )),
      }),
    ]);
    m.shell.show('form');
    flush();
    await microtasks();
    const slot = m.frame('form')!.querySelector('[data-view-footer]')!;
    expect(slot.textContent).toBe('Apply');
    expect(
      m
        .frame('form')!
        .querySelector('[data-view-body] .msui-scroll-root, [data-view-body] > :first-child')
        ?.contains(slot),
    ).toBe(false);
    m.dispose();
  });
});

describe('ShellView: tooltips', () => {
  it('installs one tooltip layer; its groups and windows are the edges a row’s card goes outside of', () => {
    const m = mountShell([probeView('home'), probeView('layers')]);
    m.shell.show('home');
    m.shell.float('layers');
    flush();
    expect(m.qa('[role="tooltip"]').length).toBe(1);
    expect(m.q('[data-group="right"]')?.hasAttribute('data-tip-edge')).toBe(true);
    expect(m.q('[data-window]')?.hasAttribute('data-tip-edge')).toBe(true);
    m.dispose();
  });

  it('none of its own with `tooltips: false`', () => {
    const m = mountShell([probeView('home')], {}, { tooltips: false });
    expect(m.q('[role="tooltip"]')).toBeNull();
    m.dispose();
  });
});
