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
import { domView } from '../../../dom/view.js';
import { NavigationManager } from '../../../managers/navigation.js';
import { ViewStateManager } from '../../../managers/view-state.js';
import { type CanvasMenuItem, ShellModel } from '../../../shell/shell.js';
import { expectNoAxeViolations } from '../../../testing/axe.js';
import { CanvasStack } from '../canvas-chrome.js';

// The collapsible canvas toolbar (shell.md §2): a chevron at the bar's edge folds it into one button whose menu lists
// the same tools, with the chevron that opens the bar again as its first row; the state is the layout's.

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

const tick = () => new Promise((r) => setTimeout(r, 0));

function setup(options: ConstructorParameters<typeof ShellModel>[0] = {}) {
  const shell = new ShellModel(options);
  const ran: string[] = [];
  const shown = new BehaviorSubject(true);
  const lit = new BehaviorSubject(false);
  const mode = new BehaviorSubject(true);
  shell.registerOverlays([{ id: 'export', title: 'Export', render: domView(() => <div />) }]);
  shell.registerToolbar([
    {
      id: 'reset',
      place: 'canvas',
      group: 'canvas',
      icon: 'focus',
      title: 'Reset camera',
      run: () => ran.push('reset'),
      menu: [{ id: 'axes', label: 'Orient axes', run: () => ran.push('axes') }],
    },
    {
      id: 'expand',
      place: 'canvas',
      group: 'canvas',
      icon: 'maximize-2',
      title: 'Expand',
      run: () => ran.push('expand'),
      menu: [
        { id: 'expand-viewport', label: 'Expand', run: () => ran.push('expand-viewport') },
        { id: 'fullscreen', label: 'Fullscreen', run: () => ran.push('fullscreen') },
      ],
    },
    {
      id: 'light',
      place: 'canvas',
      group: 'rendering',
      icon: 'lightbulb',
      title: 'Illumination',
      pressed: lit,
      run: () => {
        ran.push('light');
        lit.next(!lit.value);
      },
    },
    { id: 'export', place: 'canvas', group: 'rendering', icon: 'download', title: 'Export', overlay: 'export' },
    {
      id: 'xr',
      place: 'canvas',
      group: 'rendering',
      icon: 'glasses',
      title: 'XR',
      visible: shown,
      run: () => ran.push('xr'),
    },
    // Not on the canvas: not in the bar, not in the menu.
    { id: 'theme', icon: 'sun', title: 'Theme', run: () => ran.push('theme') },
  ]);
  const leadRows: CanvasMenuItem[] = [
    { id: 'interact', label: 'Interact', group: 'mode', choice: true, pressed: mode, run: () => ran.push('interact') },
    {
      id: 'select',
      label: 'Select',
      group: 'mode',
      choice: true,
      pressed: mode.pipe(),
      run: () => ran.push('select'),
    },
    {
      id: 'level',
      label: 'What a click selects',
      group: 'mode',
      run: () => {},
      items: of([{ id: 'residue', label: 'Residue', choice: true, pressed: of(true), run: () => ran.push('residue') }]),
    },
  ];
  const host = document.createElement('div');
  document.body.append(host);
  const dispose = render(
    () => (
      <CanvasStack
        shell={shell}
        lead={
          <div role="group" aria-label="Canvas mode">
            <button type="button" aria-label="Interact" />
            <button type="button" aria-label="Select" />
          </div>
        }
        leadMenu={of(leadRows)}
      />
    ),
    host,
  );
  const q = <T extends Element = HTMLElement>(selector: string) => host.querySelector<T>(selector);
  const qa = (selector: string) => [...host.querySelectorAll<HTMLElement>(selector)];
  return {
    shell,
    ran,
    shown,
    lit,
    host,
    dispose,
    q,
    qa,
    chevron: () => q<HTMLButtonElement>('[data-stack-toggle]')!,
    folded: () => q<HTMLButtonElement>('[data-stack-collapsed]'),
    menu: () => q('[role=menu]'),
    menus: () => qa('[role=menu]'),
    /** The labels of the rows of the open menu. */
    rows: (menu = q('[role=menu]')) =>
      [...(menu?.querySelectorAll('[role^=menuitem]') ?? [])].map((r) => r.textContent?.replace(/\s+/g, ' ').trim()),
  };
}

describe('the collapsible canvas toolbar', () => {
  it('has a chevron button, named and expanded, that folds the bar into one button', () => {
    const t = setup();
    const chevron = t.chevron();
    expect(chevron.getAttribute('aria-label')).toBe('Collapse toolbar');
    expect(chevron.getAttribute('aria-expanded')).toBe('true');
    expect(chevron.tabIndex).toBe(0);
    expect(t.folded()).toBeNull();

    chevron.click();
    expect(t.shell.toolbarCollapsed.value).toBe(true);
    expect(t.shell.viewState.layout.value.toolbar).toBe('collapsed');
    const folded = t.folded()!;
    expect(folded.getAttribute('aria-label')).toBe('Toolbar');
    expect(folded.getAttribute('aria-haspopup')).toBe('menu');
    expect(folded.getAttribute('aria-expanded')).toBe('false');
    // The bar stays, hidden and inert: neither seen, tabbed to nor read.
    const bar = chevron.closest('[hidden]')!;
    expect(bar).not.toBeNull();
    expect(bar.hasAttribute('inert')).toBe(true);
    expect(bar.querySelector('button[aria-label="Reset camera"]')).not.toBeNull();
    t.dispose();
  });

  it('lists the bar’s tools in the menu of the folded button, a checkable row for a toggle, a submenu for a menu', () => {
    const t = setup();
    // What the bar shows, by name (the xr item is visible).
    const barTools = t.qa('button[aria-label]').map((b) => b.getAttribute('aria-label'));
    expect(barTools).toEqual([
      'Collapse toolbar',
      'Interact',
      'Select',
      'Reset camera',
      'Expand',
      'Illumination',
      'Export',
      'XR',
    ]);
    t.chevron().click();
    t.folded()!.click();
    const menu = t.menu()!;
    expect(menu.getAttribute('aria-label')).toBe('Toolbar');
    expect(t.rows(menu)).toEqual([
      'Expand toolbar',
      'Interact',
      'Select',
      'What a click selects',
      'Reset camera',
      'Expand',
      'Illumination',
      'Export',
      'XR',
    ]);
    const row = (label: string) =>
      [...menu.querySelectorAll<HTMLElement>('[role^=menuitem]')].find((r) => r.textContent === label)!;
    // The mode choices are radios, a toggle and an overlay's button checkboxes, a tool with a menu a submenu.
    expect(row('Interact').getAttribute('role')).toBe('menuitemradio');
    expect(row('Illumination').getAttribute('role')).toBe('menuitemcheckbox');
    expect(row('Illumination').getAttribute('aria-checked')).toBe('false');
    expect(row('Export').getAttribute('role')).toBe('menuitemcheckbox');
    expect(row('XR').getAttribute('role')).toBe('menuitem');
    expect(row('Reset camera').getAttribute('aria-haspopup')).toBe('menu');
    expect(row('Reset camera').getAttribute('role')).toBe('menuitem');
    expect(t.rows(menu)).not.toContain('Theme');
    // The expand row is the way back, and says the bar is folded.
    expect(row('Expand toolbar').getAttribute('aria-expanded')).toBe('false');
    t.dispose();
  });

  it('opens a tool’s own menu as a submenu: its action first, unless the menu has it already', () => {
    const t = setup();
    t.chevron().click();
    t.folded()!.click();
    const rows = () => [...t.menu()!.querySelectorAll<HTMLElement>('[role^=menuitem]')];
    const open = (label: string) =>
      rows()
        .find((r) => r.textContent === label)!
        .click();
    open('Reset camera');
    expect(t.menus()).toHaveLength(2);
    expect(t.rows(t.menus()[1])).toEqual(['Reset camera', 'Orient axes']);
    t.menus()[1].querySelectorAll<HTMLElement>('[role^=menuitem]')[0].click();
    expect(t.ran).toEqual(['reset']);
    expect(t.menus()).toHaveLength(0);

    t.folded()!.click();
    open('Expand');
    // Its menu's first row is that action: not listed twice.
    expect(t.rows(t.menus()[1])).toEqual(['Expand', 'Fullscreen']);
    t.menus()[1].querySelectorAll<HTMLElement>('[role^=menuitem]')[1].click();
    expect(t.ran).toEqual(['reset', 'fullscreen']);

    // The lead’s submenu (the picking level) is one too.
    t.folded()!.click();
    open('What a click selects');
    expect(t.rows(t.menus()[1])).toEqual(['Residue']);
    t.dispose();
  });

  it('runs a toggle from its row, follows its state, and an overlay’s opens the card', () => {
    const t = setup();
    t.chevron().click();
    t.folded()!.click();
    const row = (label: string) =>
      [...t.menu()!.querySelectorAll<HTMLElement>('[role^=menuitem]')].find((r) => r.textContent === label)!;
    row('Illumination').click();
    expect(t.ran).toEqual(['light']);
    expect(t.lit.value).toBe(true);
    t.folded()!.click();
    expect(row('Illumination').getAttribute('aria-checked')).toBe('true');
    row('Export').click();
    expect(t.shell.overlay.value?.overlay.id).toBe('export');
    t.folded()!.click();
    expect(row('Export').getAttribute('aria-checked')).toBe('true');
    t.dispose();
  });

  it('follows which tools show', () => {
    const t = setup();
    t.chevron().click();
    t.folded()!.click();
    expect(t.rows()).toContain('XR');
    t.shown.next(false);
    expect(t.rows()).not.toContain('XR');
    t.shown.next(true);
    expect(t.rows()).toContain('XR');
    t.dispose();
  });

  it('opens the bar again from the first row of the menu, the focus on the chevron', async () => {
    const t = setup();
    t.chevron().focus();
    t.chevron().click();
    await tick();
    // The focus followed to the folded button.
    expect(document.activeElement).toBe(t.folded());
    t.folded()!.click();
    t.menu()!.querySelector<HTMLElement>('[role=menuitem]')!.click();
    await tick();
    expect(t.shell.toolbarCollapsed.value).toBe(false);
    expect(t.shell.viewState.layout.value.toolbar).toBe('expanded');
    expect(t.folded()).toBeNull();
    expect(t.menus()).toHaveLength(0);
    expect(t.chevron().closest('[hidden]')).toBeNull();
    expect(t.chevron().closest('[inert]')).toBeNull();
    expect(document.activeElement).toBe(t.chevron());
    t.dispose();
  });

  it('is operated from the keyboard: ↓ opens it at its first row, the arrows move, Escape closes it', async () => {
    const t = setup();
    t.chevron().focus();
    t.chevron().click();
    await tick();
    const folded = t.folded()!;
    folded.focus();
    folded.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }));
    await tick();
    const menu = t.menu()!;
    const items = [...menu.querySelectorAll<HTMLElement>('[role^=menuitem]')];
    expect(document.activeElement).toBe(items[0]);
    expect(items[0].textContent).toBe('Expand toolbar');
    const key = (k: string) =>
      (document.activeElement as HTMLElement).dispatchEvent(
        new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }),
      );
    key('ArrowDown');
    expect(document.activeElement).toBe(items[1]);
    key('End');
    expect(document.activeElement).toBe(items[items.length - 1]);
    key('ArrowDown');
    expect(document.activeElement).toBe(items[0]);
    key('ArrowUp');
    expect(document.activeElement).toBe(items[items.length - 1]);
    // → on a submenu row opens it, ← goes back to the row.
    const sub = items.find((i) => i.textContent === 'What a click selects')!;
    sub.focus();
    key('ArrowRight');
    await tick();
    expect(t.menus()).toHaveLength(2);
    key('ArrowLeft');
    expect(t.menus()).toHaveLength(1);
    expect(document.activeElement).toBe(sub);
    // Escape closes the menu, the focus goes back to the button.
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true }));
    expect(t.menus()).toHaveLength(0);
    expect(document.activeElement).toBe(folded);
    // Enter on the expand row opens the bar again.
    folded.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowUp', bubbles: true, cancelable: true }));
    await tick();
    const first = t.menu()!.querySelector<HTMLElement>('[role=menuitem]')!;
    first.focus();
    first.click();
    await tick();
    expect(document.activeElement).toBe(t.chevron());
    t.dispose();
  });

  it('starts in the state the layout has, with no motion for it, and moves when it changes', () => {
    const t = setup({ toolbar: 'collapsed' });
    expect(t.folded()).not.toBeNull();
    const frame = () => t.q('[data-stack-collapsed]')!.parentElement!;
    expect(frame().hasAttribute('data-moved')).toBe(false);
    expect(t.chevron().closest('[hidden]')).not.toBeNull();
    t.shell.toggleToolbar();
    expect(t.folded()).toBeNull();
    // The bar coming back is a change: it moves in.
    const bar = t.chevron().closest('[data-moved]');
    expect(bar).not.toBeNull();
    t.shell.toggleToolbar();
    expect(frame().hasAttribute('data-moved')).toBe(true);
    t.dispose();
  });

  it('shows no folded button where there is nothing in the stack', () => {
    const shell = new ShellModel({ toolbar: 'collapsed' });
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <CanvasStack shell={shell} />, host);
    expect(host.querySelector('[data-stack-collapsed]')).toBeNull();
    shell.registerToolbar([{ id: 'reset', place: 'canvas', icon: 'focus', title: 'Reset camera', run: () => {} }]);
    expect(host.querySelector('[data-stack-collapsed]')).not.toBeNull();
  });

  it('has no accessibility violations, folded, expanded, and with the menu open', async () => {
    const t = setup();
    await expectNoAxeViolations(t.host);
    t.chevron().click();
    await expectNoAxeViolations(t.host);
    t.folded()!.click();
    await expectNoAxeViolations(t.host);
    t.dispose();
  });
});

describe('menus in a short canvas', () => {
  it('are as tall as the room above the stack, scroll inside it, and keep the keys working', async () => {
    // jsdom has no layout: the canvas's top is at 0, the stack's at 120 (so the room is 120 less the 8px margin), the
    // rows' bottoms at 100.
    vi.spyOn(Element.prototype, 'getBoundingClientRect').mockImplementation(function (this: Element) {
      const at = (top: number, bottom = top) =>
        ({ top, bottom, left: 0, right: 400, width: 400, height: bottom - top }) as DOMRect;
      if (this.hasAttribute('data-msp-viewport')) return at(0, 300);
      if (this.matches('[class*="-stack"]')) return at(120, 154);
      if (this.matches('[role^=menuitem]')) return at(80, 100);
      return at(0);
    });
    const t = setup();
    const viewport = document.createElement('main');
    viewport.setAttribute('data-msp-viewport', '');
    t.host.parentElement!.append(viewport);
    viewport.append(t.host);
    const menuEl = () => t.q('[role=menu][aria-label=Toolbar]')!;

    // The collapsed toolbar's menu: the first row stays, the rows scroll.
    t.chevron().click();
    t.folded()!.click();
    expect(menuEl().style.getPropertyValue('--msp-menu-room')).toBe('112px');
    const scroller = menuEl().querySelector('[class*="scroll-area"]');
    expect(scroller).not.toBeNull();
    expect(scroller!.querySelectorAll('[role^=menuitem]').length).toBeGreaterThan(5);
    expect(scroller!.contains(menuEl().querySelector('[role=menuitem]'))).toBe(false);
    // A submenu has the room above its row.
    const sub = [...menuEl().querySelectorAll<HTMLElement>('[role^=menuitem]')].find(
      (r) => r.textContent === 'Reset camera',
    )!;
    sub.click();
    const submenu = t.menus()[1];
    expect(submenu.style.getPropertyValue('--msp-menu-room')).toBe('92px');
    expect(submenu.querySelector('[class*="scroll-area"]')).not.toBeNull();
    // The arrows still walk the rows, the scrolled ones too.
    submenu.querySelector<HTMLElement>('[role^=menuitem]')!.focus();
    document.activeElement!.dispatchEvent(
      new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }),
    );
    expect(document.activeElement?.textContent).toBe('Orient axes');
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true }));

    // A tool's menu in the bar is bounded the same way.
    t.shell.setToolbarCollapsed(false);
    t.q<HTMLElement>('button[aria-label="Reset camera"]')!.dispatchEvent(
      new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }),
    );
    expect((t.menu() as HTMLElement).style.getPropertyValue('--msp-menu-room')).toBe('112px');
    expect(t.menu()!.querySelector('[class*="scroll-area"]')).not.toBeNull();
    t.dispose();
    vi.restoreAllMocks();
  });
});

describe('the toolbar in the layout', () => {
  it('is saved and restored with the layout, and read as expanded from one that predates it', () => {
    const a = new ShellModel();
    a.setToolbarCollapsed(true);
    const kept = JSON.parse(JSON.stringify(a.viewState.save()));
    expect(kept.layout.toolbar).toBe('collapsed');

    const viewState = new ViewStateManager();
    const b = new ShellModel({ viewState, navigation: new NavigationManager() });
    expect(b.toolbarCollapsed.value).toBe(false);
    expect(viewState.restore(kept)).toBe(true);
    expect(b.toolbarCollapsed.value).toBe(true);

    // A layout kept before the field existed: the starting state (expanded unless the app seeds otherwise).
    delete kept.layout.toolbar;
    const c = new ViewStateManager();
    expect(c.restore(kept)).toBe(true);
    expect(c.layout.value.toolbar).toBe('expanded');
    const seeded = new ViewStateManager({ toolbar: 'collapsed' });
    expect(seeded.restore(kept)).toBe(true);
    expect(seeded.layout.value.toolbar).toBe('collapsed');
    // Another value is not a layout.
    expect(c.restore({ version: 1, layout: { docks: {}, toolbar: 'hidden' } })).toBe(false);
  });

  it('starts as the app says and Reset layout goes back to it', () => {
    const shell = new ShellModel({ toolbar: 'collapsed' });
    expect(shell.toolbarCollapsed.value).toBe(true);
    shell.toggleToolbar();
    expect(shell.toolbarCollapsed.value).toBe(false);
    shell.resetLayout();
    expect(shell.toolbarCollapsed.value).toBe(true);
    const other = new ShellModel();
    other.setToolbarCollapsed(true);
    other.resetLayout();
    expect(other.toolbarCollapsed.value).toBe(false);
  });
});
