/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { delegateEvents, DomProbe, flush, render, ViewScope } from '../../../dom/index.js';
import { styleText } from '../../../styles/registry.js';
import { state } from '../../../state.js';
import { countListeners, stubPointerCapture } from '../../../testing/gestures.js';
import {
  MenuCloseDelayMs,
  Toolbar,
  type ToolbarFace,
  type ToolbarGroup,
  type ToolbarMenuRow,
  type ToolbarProps,
} from '../toolbar.js';

beforeAll(() => {
  stubPointerCapture();
  // jsdom has no popover API: a manual popover (`Floating`) stays `display: none`, and nothing in it takes the focus.
  const style = document.createElement('style');
  style.textContent = '[popover] { display: block !important; }';
  document.head.append(style);
});

afterEach(() => {
  document.body.innerHTML = '';
  vi.useRealTimers();
  DomProbe.on = undefined;
});

function mount(props: ToolbarProps, parent?: ViewScope) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <Toolbar {...props} />, container, parent);
  const buttons = () => [...container.querySelectorAll<HTMLButtonElement>('[data-toolbar-button]')];
  const button = (label: string) => buttons().find((b) => b.getAttribute('aria-label') === label)!;
  const menu = () => container.querySelector<HTMLElement>('[role="menu"]');
  const menus = () => [...container.querySelectorAll<HTMLElement>('[role="menu"]')];
  const rows = (m: Element | null = menu()) => [...(m?.children ?? [])].filter((r) => r.hasAttribute('data-menu-row'));
  return { container, dispose, buttons, button, menu, menus, rows };
}

const key = (el: Element, k: string) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));
const over = (el: Element, from: Element | null = null) =>
  el.dispatchEvent(new PointerEvent('pointerover', { bubbles: true, pointerType: 'mouse', relatedTarget: from }));
const out = (el: Element, to: Element | null = null) =>
  el.dispatchEvent(new PointerEvent('pointerout', { bubbles: true, pointerType: 'mouse', relatedTarget: to }));
/** The microtasks: a floating menu shows, then the focus goes into it. */
const settle = async () => {
  for (let i = 0; i < 4; i++) await Promise.resolve();
  flush();
};

const levels: ToolbarMenuRow[] = [
  { id: 'atom', label: 'Atom', icon: 'pick-atom', group: 'level', choice: true },
  { id: 'residue', label: 'Residue', icon: 'pick-residue', group: 'level', choice: true },
];

function groups(over: Partial<Record<string, unknown>> = {}): ToolbarGroup[] {
  return [
    {
      id: 'view',
      label: 'View',
      items: [
        { id: 'rotate', icon: 'rotate-3d', title: 'Rotate', run: vi.fn(), ...(over.rotate as object) },
        { id: 'select', icon: 'mouse-pointer-2', title: 'Select', menu: levels, ...(over.select as object) },
      ],
    },
    {
      id: 'camera',
      label: 'Camera',
      items: [
        { id: 'focus', icon: 'focus', title: 'Focus', run: vi.fn() },
        { id: 'expand', icon: 'maximize-2', title: 'Expand', run: vi.fn(), ...(over.expand as object) },
      ],
    },
  ];
}

describe('Toolbar: the toolbar pattern', () => {
  it('is one Tab stop; the arrows along it move between buttons across blocks, wrapping; Home and End jump', () => {
    const t = mount({ groups: groups(), 'aria-label': 'Tools' });
    const bar = t.container.querySelector('[role="toolbar"]')!;
    expect(bar.getAttribute('aria-orientation')).toBe('horizontal');
    expect(t.container.querySelectorAll('[role="group"]')).toHaveLength(2);
    const [rotate, select, focus, expand] = t.buttons();
    expect(t.buttons().map((b) => b.tabIndex)).toEqual([0, -1, -1, -1]);
    rotate.focus();
    key(rotate, 'ArrowRight');
    expect(document.activeElement).toBe(select);
    key(select, 'ArrowRight');
    expect(document.activeElement).toBe(focus);
    expect(t.buttons().map((b) => b.tabIndex)).toEqual([-1, -1, 0, -1]);
    key(focus, 'End');
    expect(document.activeElement).toBe(expand);
    key(expand, 'ArrowRight');
    expect(document.activeElement).toBe(rotate);
    key(rotate, 'ArrowLeft');
    expect(document.activeElement).toBe(expand);
    key(expand, 'Home');
    expect(document.activeElement).toBe(rotate);
    // ↓ is not the bar's axis: it moves nothing here (rotate has no menu).
    key(rotate, 'ArrowDown');
    expect(document.activeElement).toBe(rotate);
  });

  it('groups by blocks by default; by separators on one surface, the arrows the same across groups', () => {
    const blocks = mount({ groups: groups() });
    expect(blocks.container.querySelector('[role="toolbar"]')!.getAttribute('data-grouping')).toBe('blocks');
    const t = mount({ groups: groups(), grouping: 'separators' });
    expect(t.container.querySelector('[role="toolbar"]')!.getAttribute('data-grouping')).toBe('separators');
    expect(t.container.querySelectorAll('[role="group"]')).toHaveLength(2);
    const [, select, focus] = t.buttons();
    select.focus();
    key(select, 'ArrowRight');
    expect(document.activeElement).toBe(focus);
  });

  it('right to left (`dir="rtl"`): ← moves on, → back', () => {
    const t = mount({ groups: groups() });
    t.container.setAttribute('dir', 'rtl');
    const [rotate, select, , expand] = t.buttons();
    rotate.focus();
    key(rotate, 'ArrowLeft');
    expect(document.activeElement).toBe(select);
    key(select, 'ArrowRight');
    expect(document.activeElement).toBe(rotate);
    key(rotate, 'ArrowRight');
    expect(document.activeElement).toBe(expand);
  });

  it('a vertical bar moves by ↓ and ↑, and skips a hidden button', async () => {
    const shown = state(true);
    const t = mount({ groups: groups({ select: { visible: shown } }), orientation: 'vertical' });
    const [rotate, select, focus] = t.buttons();
    rotate.focus();
    key(rotate, 'ArrowDown');
    expect(document.activeElement).toBe(select);
    shown.set(false);
    flush();
    expect(select.hidden).toBe(true);
    rotate.focus();
    key(rotate, 'ArrowDown');
    expect(document.activeElement).toBe(focus);
  });
});

describe('Toolbar: the menu button pattern', () => {
  it('↓ opens the menu at its first item; Esc closes it and the focus goes back to the button', async () => {
    const t = mount({ groups: groups() });
    const select = t.button('Select');
    expect(select.getAttribute('aria-haspopup')).toBe('menu');
    expect(select.getAttribute('aria-expanded')).toBe('false');
    expect(select.querySelector('[class*="mark"]')).not.toBeNull();
    // A button with a menu shows the menu, not a tooltip.
    expect(select.dataset.tip).toBeUndefined();
    select.focus();
    key(select, 'ArrowDown');
    await settle();
    expect(select.getAttribute('aria-expanded')).toBe('true');
    expect(t.menu()?.getAttribute('aria-label')).toBe('Select');
    expect(document.activeElement).toBe(t.rows()[0]);
    key(document.activeElement!, 'ArrowDown');
    expect(document.activeElement).toBe(t.rows()[1]);
    key(document.activeElement!, 'Escape');
    flush();
    expect(t.menu()).toBeNull();
    expect(select.getAttribute('aria-expanded')).toBe('false');
    expect(document.activeElement).toBe(select);
  });

  it('Tab in the menu closes it, the focus on its button first (Tab goes on from the bar); Esc on the button closes it', async () => {
    const t = mount({ groups: groups() });
    const select = t.button('Select');
    select.focus();
    key(select, 'ArrowDown');
    await settle();
    expect(document.activeElement).toBe(t.rows()[0]);
    key(t.rows()[0], 'Tab');
    expect(document.activeElement).toBe(select);
    flush();
    expect(t.menu()).toBeNull();
    over(select);
    await settle();
    expect(t.menu()).not.toBeNull();
    key(select, 'Escape');
    flush();
    expect(t.menu()).toBeNull();
  });

  it('Enter or Space (a click from the keyboard) opens a button that only holds a menu; a button with its own action runs it', async () => {
    const run = vi.fn();
    const t = mount({ groups: groups({ expand: { run, menu: levels } }) });
    const select = t.button('Select');
    select.focus();
    select.click();
    await settle();
    expect(document.activeElement).toBe(t.rows()[0]);
    const expand = t.button('Expand');
    expand.click();
    await settle();
    expect(run).toHaveBeenCalledWith(expand);
    expect(t.menu()).toBeNull();
  });

  it('opens on hover and closes a moment after the pointer leaves, staying while the pointer is over the menu', async () => {
    vi.useFakeTimers();
    const t = mount({ groups: groups() });
    const select = t.button('Select');
    over(select);
    await settle();
    expect(t.menu()).not.toBeNull();
    expect(select.getAttribute('aria-expanded')).toBe('true');
    // Out of the button into the menu: it stays.
    out(select, t.rows()[0]);
    over(t.rows()[0], select);
    vi.advanceTimersByTime(MenuCloseDelayMs * 2);
    expect(t.menu()).not.toBeNull();
    // Out of the menu, and out of the bar: it closes after the delay.
    out(t.rows()[0], document.body);
    vi.advanceTimersByTime(MenuCloseDelayMs - 10);
    expect(t.menu()).not.toBeNull();
    // Back in time: kept.
    over(select, document.body);
    vi.advanceTimersByTime(MenuCloseDelayMs * 2);
    expect(t.menu()).not.toBeNull();
    out(select, document.body);
    vi.advanceTimersByTime(MenuCloseDelayMs);
    flush();
    expect(t.menu()).toBeNull();
    expect(select.getAttribute('aria-expanded')).toBe('false');
  });

  it('a tap opens the menu with the button’s own action first, and a second tap closes it', async () => {
    const run = vi.fn();
    const t = mount({ groups: groups({ select: { run } }) });
    const select = t.button('Select');
    select.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, pointerType: 'touch' }));
    select.dispatchEvent(new MouseEvent('click', { bubbles: true, detail: 1 }));
    await settle();
    expect(run).not.toHaveBeenCalled();
    const rows = t.rows();
    expect(rows.map((r) => r.textContent)).toEqual(['Select', 'Atom', 'Residue']);
    (rows[0] as HTMLElement).click();
    expect(run).toHaveBeenCalledTimes(1);
    expect(t.menu()).toBeNull();
    select.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, pointerType: 'touch' }));
    select.dispatchEvent(new MouseEvent('click', { bubbles: true, detail: 1 }));
    await settle();
    expect(t.menu()).not.toBeNull();
    select.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, pointerType: 'touch' }));
    select.dispatchEvent(new MouseEvent('click', { bubbles: true, detail: 1 }));
    await settle();
    expect(t.menu()).toBeNull();
  });

  it('closes on a press outside it and its button', async () => {
    const t = mount({ groups: groups() });
    over(t.button('Select'));
    await settle();
    t.menu()!.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(t.menu()).not.toBeNull();
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    flush();
    expect(t.menu()).toBeNull();
  });
});

describe('Toolbar: menu rows', () => {
  it('a choice row is a menuitemradio, checked by its state; a toggle a menuitemcheckbox with a check', async () => {
    const level = state('atom');
    const ghost = state(false);
    const pick = vi.fn((id: string) => level.set(id));
    const rows: ToolbarMenuRow[] = ['atom', 'residue', 'chain'].map((id) => ({
      id,
      label: id,
      group: 'level',
      choice: true,
      pressed: derive(level, (v) => v === id),
      run: () => pick(id),
    }));
    rows.push({
      id: 'ghost',
      label: 'Ghost',
      group: 'all',
      heading: 'All instances',
      pressed: ghost,
      run: () => ghost.set(!ghost.value),
    });
    const t = mount({ groups: groups({ select: { menu: rows } }) });
    over(t.button('Select'));
    await settle();
    const [atom, residue, , g] = t.rows();
    expect(atom.getAttribute('role')).toBe('menuitemradio');
    expect(atom.getAttribute('aria-checked')).toBe('true');
    expect(residue.getAttribute('aria-checked')).toBe('false');
    expect(g.getAttribute('role')).toBe('menuitemcheckbox');
    expect(g.getAttribute('aria-checked')).toBe('false');
    expect(g.querySelector('[class*="check"]')).not.toBeNull();
    // The second group is headed, not divided.
    expect(t.menu()!.querySelector('[class*="heading"]')?.textContent).toBe('All instances');
    expect(t.menu()!.querySelector('[role="separator"]')).toBeNull();
    (residue as HTMLElement).click();
    expect(pick).toHaveBeenCalledWith('residue');
    expect(t.menu()).toBeNull();
    over(t.button('Select'));
    await settle();
    expect(t.rows()[0].getAttribute('aria-checked')).toBe('false');
    expect(t.rows()[1].getAttribute('aria-checked')).toBe('true');
  });

  it('a submenu opens beside its row by → at its first item, ← closes it back to its row; hover opens and moves on', async () => {
    const run = vi.fn();
    const rows: ToolbarMenuRow[] = [
      { id: 'a', label: 'Plain', group: 'x' },
      { id: 'more', label: 'More', group: 'y', items: [{ id: 'deep', label: 'Deep', run }] },
    ];
    const t = mount({ groups: groups({ select: { menu: rows } }) });
    const select = t.button('Select');
    select.focus();
    key(select, 'ArrowDown');
    await settle();
    expect(t.menu()!.querySelector('[role="separator"]')).not.toBeNull();
    const [, more] = t.rows();
    expect(more.getAttribute('aria-haspopup')).toBe('menu');
    key(t.rows()[0], 'ArrowDown');
    expect(document.activeElement).toBe(more);
    key(more, 'ArrowRight');
    await settle();
    expect(t.menus()).toHaveLength(2);
    expect(more.getAttribute('aria-expanded')).toBe('true');
    const sub = t.menus()[1];
    expect(sub.getAttribute('aria-label')).toBe('More');
    expect(document.activeElement).toBe(t.rows(sub)[0]);
    key(document.activeElement!, 'ArrowLeft');
    flush();
    expect(t.menus()).toHaveLength(1);
    expect(document.activeElement).toBe(more);
    expect(t.menu()).not.toBeNull();
    // Hover: the submenu row opens it, another row closes it.
    over(more);
    flush();
    expect(t.menus()).toHaveLength(2);
    over(t.rows()[0]);
    flush();
    expect(t.menus()).toHaveLength(1);
    // A sub row runs and closes the whole menu.
    over(more);
    flush();
    (t.rows(t.menus()[1])[0] as HTMLElement).click();
    expect(run).toHaveBeenCalledTimes(1);
    expect(t.menu()).toBeNull();
  });

  it('Esc in a submenu closes only it; Tab there closes the whole menu; ← or → back follows the side it opened on', async () => {
    const rows: ToolbarMenuRow[] = [
      { id: 'more', label: 'More', items: [{ id: 'deep', label: 'Deep', run: () => {} }] },
    ];
    const t = mount({ groups: groups({ select: { menu: rows } }) });
    const select = t.button('Select');
    select.focus();
    key(select, 'ArrowDown');
    await settle();
    const [more] = t.rows();
    key(more, 'ArrowRight');
    await settle();
    expect(t.menus()).toHaveLength(2);
    key(document.activeElement!, 'Escape');
    flush();
    expect(t.menus()).toHaveLength(1);
    expect(document.activeElement).toBe(more);
    key(more, 'ArrowRight');
    await settle();
    key(document.activeElement!, 'Tab');
    expect(document.activeElement).toBe(select);
    flush();
    expect(t.menu()).toBeNull();

    // At the window's right edge the submenu flips to the left: → is back, ← does nothing.
    key(select, 'ArrowDown');
    await settle();
    const row = t.rows()[0] as HTMLElement;
    row.getBoundingClientRect = () =>
      ({ top: 0, bottom: 20, left: 1000, right: 1020, width: 20, height: 20, x: 1000, y: 0 }) as DOMRect;
    key(row, 'ArrowRight');
    await settle();
    const sub = t.menus()[1];
    expect(sub.closest<HTMLElement>('[popover]')!.dataset.side).toBe('left');
    key(document.activeElement!, 'ArrowLeft');
    flush();
    expect(t.menus()).toHaveLength(2);
    key(document.activeElement!, 'ArrowRight');
    flush();
    expect(t.menus()).toHaveLength(1);
    expect(document.activeElement).toBe(row);
  });

  it('a disabled button is dimmed by the toolbar itself (its text token is too near the enabled one)', () => {
    // the registered sheet keeps the nesting: the first `&[aria-disabled="true"]` block after the button's class
    const css = styleText();
    const at = css.indexOf('.msui-toolbar-button {');
    const rule = /&\[aria-disabled="true"\] \{([^}]*)\}/.exec(css.slice(at))?.[1];
    expect(rule).toMatch(/opacity:\s*0?\.[3-6]/);
  });

  it('a disabled button and a disabled row do not run; their reasons are in their tooltips', async () => {
    const run = vi.fn();
    const rowRun = vi.fn();
    const off = state<boolean | string>('Nothing to focus');
    const t = mount({
      groups: [
        {
          id: 'g',
          items: [
            { id: 'f', icon: 'focus', title: 'Focus', disabled: off, run },
            {
              id: 'm',
              icon: 'ellipsis',
              title: 'More',
              menu: [
                { id: 'x', label: 'Export', disabled: 'No structure', run: rowRun },
                { id: 'y', label: 'Other', run: () => {} },
              ],
            },
          ],
        },
      ],
    });
    const focus = t.button('Focus');
    expect(focus.getAttribute('aria-disabled')).toBe('true');
    expect(focus.dataset.tip).toBe('Focus: Nothing to focus');
    // The reason, for assistive technology too.
    expect(focus.getAttribute('aria-description')).toBe('Nothing to focus');
    focus.click();
    expect(run).not.toHaveBeenCalled();
    off.set(false);
    flush();
    expect(focus.hasAttribute('aria-disabled')).toBe(false);
    expect(focus.hasAttribute('aria-description')).toBe(false);
    expect(focus.dataset.tip).toBe('Focus');
    focus.click();
    expect(run).toHaveBeenCalledTimes(1);
    const more = t.button('More');
    more.focus();
    key(more, 'ArrowDown');
    await settle();
    const [x, y] = t.rows();
    expect(x.getAttribute('aria-disabled')).toBe('true');
    expect((x as HTMLElement).dataset.tip).toBe('No structure');
    expect(x.getAttribute('aria-description')).toBe('No structure');
    // It opens on the first usable row; the arrows reach the disabled one (APG menus: its reason is read).
    expect(document.activeElement).toBe(y);
    key(y, 'ArrowUp');
    expect(document.activeElement).toBe(x);
    (x as HTMLElement).click();
    expect(rowRun).not.toHaveBeenCalled();
    expect(t.menu()).not.toBeNull();
  });

  it('a face changes the icon and the name', () => {
    const face = state<ToolbarFace>({ icon: 'sun', title: 'Light theme' });
    const t = mount({ groups: [{ id: 'g', items: [{ id: 't', icon: 'sun', title: 'Theme', face, run: () => {} }] }] });
    const [b] = t.buttons();
    expect(b.getAttribute('aria-label')).toBe('Light theme');
    const svg = b.querySelector('svg')!;
    face.set({ icon: 'moon', title: 'Dark theme' });
    flush();
    expect(b.getAttribute('aria-label')).toBe('Dark theme');
    expect(b.dataset.tip).toBe('Dark theme');
    expect(b.querySelector('svg')).not.toBe(svg);
  });

  it('a toggle is aria-pressed by its state', () => {
    const on = state(false);
    const t = mount({ groups: groups({ rotate: { pressed: on } }) });
    const rotate = t.button('Rotate');
    expect(rotate.getAttribute('aria-pressed')).toBe('false');
    on.set(true);
    flush();
    expect(rotate.getAttribute('aria-pressed')).toBe('true');
  });
});

describe('Toolbar: runtime', () => {
  it('opening a menu renders the button that opens it and the one that closes, no other', async () => {
    const items = ['a', 'b', 'c', 'd', 'e'].map((id) => ({
      id,
      icon: 'ellipsis' as const,
      title: id,
      menu: [{ id: `${id}1`, label: `${id} 1`, run: () => {} }],
    }));
    // Set before rendering: the probe wraps what is made after it.
    let renders = 0;
    DomProbe.on = (kind) => {
      if (kind === 'render') renders++;
    };
    const t = mount({ groups: [{ id: 'g', items }] });
    renders = 0;
    over(t.button('a'));
    await settle();
    expect(renders).toBe(1);
    renders = 0;
    out(t.button('a'), t.button('c'));
    over(t.button('c'), t.button('a'));
    await settle();
    expect(t.button('c').getAttribute('aria-expanded')).toBe('true');
    expect(renders).toBe(2);
  });

  it('leaves no listener after dispose, and none on the document while paused', async () => {
    const container = document.createElement('div');
    document.body.append(container);
    delegateEvents(container);
    delegateEvents(document);
    const listeners = countListeners();
    try {
      const parent = new ViewScope();
      const dispose = render(() => <Toolbar groups={groups()} />, container, parent);
      const select = container.querySelector<HTMLElement>('[aria-haspopup="menu"]')!;
      over(select);
      await settle();
      expect(container.querySelector('[role="menu"]')).not.toBeNull();
      expect(listeners.count).toBeGreaterThan(0);
      parent.setActive(false);
      flush();
      // Paused: the menu closed, and nothing is attached.
      expect(container.querySelector('[role="menu"]')).toBeNull();
      expect(listeners.count).toBe(0);
      parent.setActive(true);
      expect(listeners.count).toBeGreaterThan(0);
      over(select);
      await settle();
      dispose();
      expect(listeners.count).toBe(0);
    } finally {
      listeners.stop();
    }
  });
});
