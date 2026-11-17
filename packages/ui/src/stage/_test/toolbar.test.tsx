/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../dom/index.js';
import type { ToolbarGroup, ToolbarMenuRow } from '../../components/toolbar/toolbar.js';
import { stubPointerCapture } from '../../testing/gestures.js';
import { StageModel } from '../stage-model.js';
import { collapsedRows, StageToolbar, type StageToolbarProps } from '../toolbar.js';

// The stage's toolbar (stage §3.2): menus open upward, bounded by the room above (measured as one opens); folded, one
// button whose menu lists the rows (toggles checkable, actions, menus as submenus led by their own action).

beforeAll(() => {
  stubPointerCapture();
  // jsdom has no popover API: a manual popover (`Floating`) stays `display: none`.
  const style = document.createElement('style');
  style.textContent = '[popover] { display: block !important; }';
  document.head.append(style);
});

// Every mount disposed after its test: an open menu's document listeners would outlive it (files share a worker).
const mounted: (() => void)[] = [];
afterEach(() => {
  for (const dispose of mounted.splice(0)) dispose();
  document.body.innerHTML = '';
});

/** jsdom has no ResizeObserver: a stub keeping the observers by element. */
const observers = new Map<Element, Set<ResizeObserverCallback>>();
beforeEach(() => {
  observers.clear();
  vi.stubGlobal(
    'ResizeObserver',
    class {
      readonly #fn: ResizeObserverCallback;
      constructor(fn: ResizeObserverCallback) {
        this.#fn = fn;
      }
      observe(el: Element) {
        let s = observers.get(el);
        if (!s) observers.set(el, (s = new Set()));
        s.add(this.#fn);
      }
      unobserve() {}
      disconnect() {
        for (const s of observers.values()) s.delete(this.#fn);
      }
    },
  );
});
const resize = (el: Element, width: number) => {
  const entry = { target: el, contentRect: { width, height: 400 } } as unknown as ResizeObserverEntry;
  for (const fn of [...(observers.get(el) ?? [])]) fn([entry], {} as ResizeObserver);
  flush();
};

const settle = async () => {
  for (let i = 0; i < 4; i++) await Promise.resolve();
  flush();
};
const over = (el: Element) =>
  el.dispatchEvent(new PointerEvent('pointerover', { bubbles: true, pointerType: 'mouse' }));

const spin = new State(false);
const reset = vi.fn();
const projection = vi.fn();
function groups(): ToolbarGroup[] {
  return [
    {
      id: 'camera',
      label: 'Camera',
      items: [
        { id: 'reset', icon: 'refresh-cw', title: 'Reset', run: reset },
        { id: 'spin', icon: 'rotate-3d', title: 'Spin', pressed: spin, run: () => spin.set(!spin.value) },
        {
          id: 'projection',
          icon: 'camera',
          title: 'Projection',
          run: projection,
          menu: [
            { id: 'persp', label: 'Perspective', run: () => {} },
            { id: 'ortho', label: 'Orthographic', run: () => {} },
          ],
        },
      ],
    },
    {
      id: 'view',
      label: 'View',
      items: [
        { id: 'hidden', icon: 'eye', title: 'Hidden', visible: false },
        { id: 'expand', icon: 'maximize', title: 'Expand', disabled: 'Not here', run: () => {} },
      ],
    },
  ];
}

function mount(props: Partial<StageToolbarProps>, stageTop = 100, buttonTop = 500) {
  const stage = document.createElement('div');
  stage.setAttribute('data-msui-stage', '');
  document.body.append(stage);
  vi.spyOn(stage, 'getBoundingClientRect').mockReturnValue({ top: stageTop } as DOMRect);
  const dispose = render(() => <StageToolbar groups={groups()} {...props} />, stage);
  mounted.push(dispose);
  for (const b of stage.querySelectorAll('[data-toolbar-button]'))
    vi.spyOn(b, 'getBoundingClientRect').mockReturnValue({ top: buttonTop } as DOMRect);
  const root = stage.querySelector<HTMLElement>('[data-msui-stage-toolbar]')!;
  const button = (label: string) => root.querySelector<HTMLElement>(`[data-toolbar-button][aria-label="${label}"]`)!;
  const rows = (m: Element | null) => [...(m?.children ?? [])].filter((r) => r.hasAttribute('data-menu-row'));
  return { stage, root, dispose, button, rows, menu: () => root.querySelector<HTMLElement>('[role="menu"]') };
}

describe('StageToolbar', () => {
  it('opens menus upward, as tall as the room above in the stage', async () => {
    const t = mount({});
    over(t.button('Projection'));
    await settle();
    expect(t.menu()).not.toBeNull();
    expect(t.menu()!.closest('[popover]')?.getAttribute('data-placement')).toBe('above');
    // The button's top less the stage's, the menu's gap and the edge: 500 - 100 - 16.
    expect(t.root.style.getPropertyValue('--msui-stage-menu-room')).toBe('384px');
  });

  it('expanded: the row shows, scrolling across; no folded button', () => {
    const t = mount({ collapsed: false });
    expect(t.button('Reset')).not.toBeNull();
    expect(t.root.querySelector('[aria-label="Toolbar"][data-toolbar-button]')).toBeNull();
  });

  it('folded: one button; the row stays mounted, hidden and inert', async () => {
    const t = mount({ collapsed: true });
    const bar = t.root.firstElementChild as HTMLElement;
    expect(bar.hidden).toBe(true);
    expect(bar.hasAttribute('inert')).toBe(true);
    const folded = t.root.querySelector<HTMLElement>('[data-toolbar-button][aria-label="Toolbar"]')!;
    vi.spyOn(folded, 'getBoundingClientRect').mockReturnValue({ top: 300 } as DOMRect);
    over(folded);
    await settle();
    const menu = t.menu()!;
    const rows = t.rows(menu);
    expect(rows.map((r) => r.textContent)).toEqual(['Reset', 'Spin', 'Projection', 'Expand']);
    expect(rows[1].getAttribute('role')).toBe('menuitemcheckbox');
    expect(rows[1].getAttribute('aria-checked')).toBe('false');
    expect(rows[2].getAttribute('aria-haspopup')).toBe('menu');
    expect(rows[3].getAttribute('aria-disabled')).toBe('true');
    expect(t.root.style.getPropertyValue('--msui-stage-menu-room')).toBe('184px');
    (rows[1] as HTMLElement).click();
    expect(spin.value).toBe(true);
    spin.set(false);
  });

  it('folded: a submenu opens under the pointer; the open menu is neither rebuilt nor re-measured', async () => {
    const live = new State(groups());
    const t = mount({ collapsed: true, groups: live });
    const folded = t.root.querySelector<HTMLElement>('[data-toolbar-button][aria-label="Toolbar"]')!;
    vi.spyOn(folded, 'getBoundingClientRect').mockReturnValue({ top: 300 } as DOMRect);
    over(folded);
    await settle();
    const menu = t.menu()!;
    const rows = t.rows(menu);
    const projection = rows[2] as HTMLElement;
    expect(projection.getAttribute('aria-haspopup')).toBe('menu');
    expect(projection.getAttribute('aria-expanded')).toBe('false');
    // A submenu's row is a menu button inside the toolbar too: the pointer on it opens its submenu, nothing else.
    vi.spyOn(projection, 'getBoundingClientRect').mockReturnValue({ top: 200 } as DOMRect);
    over(projection);
    await settle();
    expect(projection.getAttribute('aria-expanded')).toBe('true');
    expect(t.root.querySelectorAll('[role="menu"]')).toHaveLength(2);
    expect(t.menu()).toBe(menu);
    expect(t.rows(menu)[2]).toBe(projection);
    expect(t.root.style.getPropertyValue('--msui-stage-menu-room')).toBe('184px');
    // The groups changing while it is open: its rows stay (built again as it next opens).
    live.set(groups().slice(0, 1));
    flush();
    expect(t.rows(t.menu())).toHaveLength(4);
    expect(projection.getAttribute('aria-expanded')).toBe('true');
    t.dispose();
  });

  it('takes its text from `labels` and its folded icon from `foldedIcon`', () => {
    const t = mount({
      collapsed: true,
      foldedIcon: 'settings',
      labels: { folded: 'Werkzeuge', foldedToolbar: 'Leiste', tools: 'Bühne' },
    });
    const folded = t.root.querySelector<HTMLElement>('[data-toolbar-button][aria-label="Werkzeuge"]');
    expect(folded).not.toBeNull();
    expect(folded!.closest('[role="toolbar"]')?.getAttribute('aria-label')).toBe('Leiste');
    expect(t.root.querySelector('[role="toolbar"][aria-label="Bühne"]')).not.toBeNull();
    expect(folded!.querySelector('svg, [data-icon]')).not.toBeNull();
  });

  it('auto without a stage: folds while the element it is placed in is narrow', async () => {
    const t = mount({});
    const bar = t.root.firstElementChild as HTMLElement;
    await settle();
    expect(bar.hidden).toBe(false);
    resize(t.stage, 320);
    expect(bar.hidden).toBe(true);
    resize(t.stage, 900);
    expect(bar.hidden).toBe(false);
  });

  it('auto: folds while the stage is compact', () => {
    const stage = new StageModel({ empty: false });
    const t = mount({ stage });
    const bar = t.root.firstElementChild as HTMLElement;
    expect(bar.hidden).toBe(false);
    stage.setSize(320, 400);
    flush();
    expect(bar.hidden).toBe(true);
    expect(t.root.querySelector('[data-toolbar-button][aria-label="Toolbar"]')).not.toBeNull();
    stage.setSize(900, 400);
    flush();
    expect(bar.hidden).toBe(false);
    stage.dispose();
  });

  it('on touch, a tap opens a menu with the button’s own action first', async () => {
    const t = mount({ collapsed: false });
    const b = t.button('Projection');
    b.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, pointerType: 'touch' }));
    b.dispatchEvent(new MouseEvent('click', { bubbles: true, detail: 1 }));
    await settle();
    expect(projection).not.toHaveBeenCalled();
    expect(t.rows(t.menu()).map((r) => r.textContent)).toEqual(['Projection', 'Perspective', 'Orthographic']);
  });
});

describe('collapsedRows', () => {
  it('lists the lead, then each block’s shown items; a menu is a submenu led by its own action', () => {
    const anchor = document.createElement('button');
    const lead: ToolbarMenuRow[] = [{ id: 'mode', label: 'Select mode', run: () => {} }];
    const rows = collapsedRows(groups(), lead, () => anchor);
    expect(rows.map((r) => r.label)).toEqual(['Select mode', 'Reset', 'Spin', 'Projection', 'Expand']);
    expect(rows[0].group).not.toBe(rows[1].group);
    expect(rows[2].pressed).toBe(spin);
    const sub = rows[3].items!;
    expect(sub.map((r) => r.label)).toEqual(['Projection', 'Perspective', 'Orthographic']);
    sub[0].run!();
    expect(projection).toHaveBeenCalledWith(anchor);
    expect(rows[4].disabled).toBe('Not here');
  });

  it('follows a face and a menu given as states; no own action where a row is it already', () => {
    const face = new State({ icon: 'sun' as const, title: 'Light' });
    const menu = new State<readonly ToolbarMenuRow[]>([{ id: 'l', label: 'Light', run: () => {} }]);
    const rows = collapsedRows(
      [
        {
          id: 'g',
          items: [{ id: 'theme', icon: 'moon', title: 'Theme', face: derive(face, (f) => f), menu, run: () => {} }],
        },
      ],
      [],
      () => document.body,
    );
    expect(rows[0].label).toBe('Light');
    expect(rows[0].items!.map((r) => r.label)).toEqual(['Light']);
  });
});
