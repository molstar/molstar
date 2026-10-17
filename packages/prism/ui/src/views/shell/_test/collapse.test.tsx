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

// The compact canvas toolbar (shell.md §2): `ShellOptions.toolbar` ('auto' by default: compact mode, 'expanded',
// 'collapsed') makes it one button whose menu lists the same tools; the host's, with no control for the person.

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

/** A shell with the toolbar collapsed (the default here) or not. */
function setup(options: ConstructorParameters<typeof ShellModel>[0] = { toolbar: 'collapsed' }) {
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
    folded: () => q<HTMLButtonElement>('[data-stack-collapsed]'),
    menu: () => q('[role=menu]'),
    menus: () => qa('[role=menu]'),
    /** The labels of the rows of the open menu. */
    rows: (menu = q('[role=menu]')) =>
      [...(menu?.querySelectorAll('[role^=menuitem]') ?? [])].map((r) => r.textContent?.replace(/\s+/g, ' ').trim()),
  };
}

describe('the compact canvas toolbar', () => {
  it('is one named button with a menu, and the bar hidden; it has no chevron and no way back', () => {
    const t = setup();
    const folded = t.folded()!;
    expect(folded.getAttribute('aria-label')).toBe('Toolbar');
    expect(folded.getAttribute('aria-haspopup')).toBe('menu');
    expect(folded.getAttribute('aria-expanded')).toBe('false');
    expect(t.shell.toolbarCollapsed.value).toBe(true);
    // The bar is there, hidden and inert: neither seen, tabbed to nor read.
    expect(t.q('[hidden]')).not.toBeNull();
    expect(t.q('[hidden]')!.hasAttribute('inert')).toBe(true);
    expect(t.q('[data-stack-toggle]')).toBeNull();
    expect(t.q('button[aria-label="Collapse toolbar"]')).toBeNull();
    folded.click();
    expect(t.rows()).not.toContain('Expand toolbar');
    t.dispose();
  });

  it('drops the row of an item that left the toolbar, with what it listens to', async () => {
    const t = setup();
    const lit = new BehaviorSubject(false);
    const drop = t.shell.registerToolbar([
      {
        id: 'gone',
        place: 'canvas',
        group: 'rendering',
        icon: 'lightbulb',
        title: 'Gone',
        pressed: lit,
        run: () => {},
      },
    ]);
    await tick();
    t.folded()!.click();
    expect(t.rows()).toContain('Gone');
    expect(lit.observed).toBe(true);
    drop();
    await tick();
    expect(lit.observed).toBe(false);
    t.dispose();
  });

  it('is the full bar when expanded, with no chevron and no folded button', () => {
    const t = setup({ toolbar: 'expanded' });
    expect(t.shell.toolbarCollapsed.value).toBe(false);
    expect(t.folded()).toBeNull();
    expect(t.q('[data-stack-toggle]')).toBeNull();
    expect(t.qa('button[aria-label]').map((b) => b.getAttribute('aria-label'))).toEqual([
      'Interact',
      'Select',
      'Reset camera',
      'Expand',
      'Illumination',
      'Export',
      'XR',
    ]);
    t.dispose();
  });

  it('lists the tools in the menu, a checkable row for a toggle, a submenu for a menu', () => {
    const t = setup();
    t.folded()!.click();
    const menu = t.menu()!;
    expect(menu.getAttribute('aria-label')).toBe('Toolbar');
    expect(t.rows(menu)).toEqual([
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
    t.dispose();
  });

  it('opens a tool’s own menu as a submenu: its action first, unless the menu has it already', () => {
    const t = setup();
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
    t.folded()!.click();
    expect(t.rows()).toContain('XR');
    t.shown.next(false);
    expect(t.rows()).not.toContain('XR');
    t.shown.next(true);
    expect(t.rows()).toContain('XR');
    t.dispose();
  });

  it('is operated from the keyboard: ↓ opens it at its first row, the arrows move, Escape closes it', async () => {
    const t = setup();
    const folded = t.folded()!;
    folded.focus();
    folded.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }));
    await tick();
    const menu = t.menu()!;
    const items = [...menu.querySelectorAll<HTMLElement>('[role^=menuitem]')];
    expect(document.activeElement).toBe(items[0]);
    expect(items[0].textContent).toBe('Interact');
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
    t.dispose();
  });

  it('shows no button where there is nothing in the stack', () => {
    const shell = new ShellModel({ toolbar: 'collapsed' });
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <CanvasStack shell={shell} />, host);
    expect(host.querySelector('[data-stack-collapsed]')).toBeNull();
    shell.registerToolbar([{ id: 'reset', place: 'canvas', icon: 'focus', title: 'Reset camera', run: () => {} }]);
    expect(host.querySelector('[data-stack-collapsed]')).not.toBeNull();
  });

  it('has no accessibility violations, closed and with the menu open', async () => {
    const t = setup();
    await expectNoAxeViolations(t.host);
    t.folded()!.click();
    await expectNoAxeViolations(t.host);
    t.dispose();
  });
});

describe('the toolbar option', () => {
  const widths = (shell: ShellModel, w: number) => shell.setWidth(w);

  it('auto (the default) is the bar in the desktop mode and the button in the compact one', () => {
    const desktop = new ShellModel({ compactWidth: 600 });
    widths(desktop, 1000);
    expect(desktop.mode.value).toBe('desktop');
    expect(desktop.toolbarCollapsed.value).toBe(false);
    const compact = new ShellModel({ toolbar: 'auto', compactWidth: 600 });
    widths(compact, 400);
    expect(compact.mode.value).toBe('compact');
    expect(compact.toolbarCollapsed.value).toBe(true);
    // An embedded layout is compact below the width too, and the bar above it.
    const embedded = new ShellModel({ mode: 'embedded', compactWidth: 600 });
    widths(embedded, 1000);
    expect(embedded.toolbarCollapsed.value).toBe(false);
    widths(embedded, 400);
    expect(embedded.toolbarCollapsed.value).toBe(true);
  });

  it('auto follows the mode live as the shell resizes across compactWidth', () => {
    const t = setup({ compactWidth: 600 });
    t.shell.setWidth(1000);
    expect(t.folded()).toBeNull();
    expect(t.q('button[aria-label="Reset camera"]')!.closest('[hidden]')).toBeNull();
    t.shell.setWidth(400);
    expect(t.folded()).not.toBeNull();
    expect(t.q('button[aria-label="Reset camera"]')!.closest('[hidden]')).not.toBeNull();
    t.shell.setWidth(900);
    expect(t.folded()).toBeNull();
    expect(t.q('button[aria-label="Reset camera"]')!.closest('[hidden]')).toBeNull();
    t.dispose();
  });

  it('expanded and collapsed force one or the other, whatever the mode', () => {
    for (const width of [1000, 400]) {
      const expanded = new ShellModel({ toolbar: 'expanded', compactWidth: 600 });
      expanded.setWidth(width);
      expect(expanded.toolbarCollapsed.value).toBe(false);
      const collapsed = new ShellModel({ toolbar: 'collapsed', compactWidth: 600 });
      collapsed.setWidth(width);
      expect(collapsed.toolbarCollapsed.value).toBe(true);
    }
  });

  it('is not part of the layout: nothing kept, Reset layout leaves it, an old layout’s value is ignored', () => {
    const shell = new ShellModel({ toolbar: 'collapsed' });
    const kept = JSON.parse(JSON.stringify(shell.viewState.save()));
    expect(kept.layout.toolbar).toBeUndefined();
    expect(shell.viewState.layout.value).not.toHaveProperty('toolbar');
    shell.resetLayout();
    expect(shell.toolbarCollapsed.value).toBe(true);
    const viewState = new ViewStateManager();
    expect(viewState.restore({ ...kept, layout: { ...kept.layout, toolbar: 'collapsed' } })).toBe(true);
    expect(viewState.restore({ ...kept, layout: { ...kept.layout, toolbar: 'hidden' } })).toBe(true);
    expect(viewState.layout.value).not.toHaveProperty('toolbar');
    const other = new ShellModel({ viewState, navigation: new NavigationManager() });
    expect(other.toolbarCollapsed.value).toBe(false);
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
    const mountIn = (t: ReturnType<typeof setup>) => {
      const viewport = document.createElement('main');
      viewport.setAttribute('data-msp-viewport', '');
      t.host.parentElement!.append(viewport);
      viewport.append(t.host);
    };
    const t = setup();
    mountIn(t);
    const menuEl = () => t.q('[role=menu][aria-label=Toolbar]')!;

    // The compact toolbar's menu: its rows scroll.
    t.folded()!.click();
    expect(menuEl().style.getPropertyValue('--msp-menu-room')).toBe('112px');
    const scroller = menuEl().querySelector('[class*="scroll-area"]');
    expect(scroller).not.toBeNull();
    expect(scroller!.querySelectorAll('[role^=menuitem]').length).toBeGreaterThan(5);
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
    t.dispose();

    // A tool's menu in the full bar is bounded the same way.
    const bar = setup({ toolbar: 'expanded' });
    mountIn(bar);
    bar
      .q<HTMLElement>('button[aria-label="Reset camera"]')!
      .dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }));
    expect((bar.menu() as HTMLElement).style.getPropertyValue('--msp-menu-room')).toBe('112px');
    expect(bar.menu()!.querySelector('[class*="scroll-area"]')).not.toBeNull();
    bar.dispose();
    vi.restoreAllMocks();
  });
});
