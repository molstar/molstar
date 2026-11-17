/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { flush, flushFrame } from '../../dom/index.js';
import { pointer, stubPointerCapture } from '../../testing/gestures.js';
import { findView } from '../layout.js';
import { key, microtasks, mountShell, press, probeView, stubResizeObserver } from './fixture.js';

beforeAll(() => stubPointerCapture());

afterEach(() => {
  document.body.innerHTML = '';
});

const menuItem = (label: string) =>
  [...document.querySelectorAll<HTMLElement>('[role="menuitem"]')].find((m) => m.textContent?.includes(label));
const button = (root: ParentNode, title: string) =>
  root.querySelector<HTMLButtonElement>(`button[aria-label="${title}"]`);
const rectOf = (el: HTMLElement) => ({
  x: Number.parseFloat(el.style.left),
  y: Number.parseFloat(el.style.top),
  width: Number.parseFloat(el.style.width),
  height: Number.parseFloat(el.style.height),
});

describe('ShellView: split below', () => {
  it('a lower group shows labelled tabs and its own frame; folding it hides its body', () => {
    const a = probeView('layers');
    const b = probeView('inspect');
    const c = probeView('plot');
    const m = mountShell([a, b, c]);
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    flush();
    const lower = findView(m.shell.layout.value, 'inspect')!.group;
    expect(lower.id).not.toBe('right');
    m.shell.move('plot', { kind: 'group', id: lower.id });
    m.shell.show('plot');
    flush();
    const g = m.q(`[data-group="${lower.id}"]`)!;
    const tabs = [...g.querySelectorAll('[role="tab"]')];
    expect(tabs.map((t) => t.textContent)).toEqual(['Inspect', 'Plot']);
    expect(tabs.map((t) => t.getAttribute('aria-selected'))).toEqual(['false', 'true']);
    expect(g.querySelector('[data-msui-view="plot"] header')?.getAttribute('data-head')).toBe('header');
    // The upper group keeps its frame (the rail's).
    expect(m.frame('layers')?.querySelector('header')?.getAttribute('data-head')).toBe('rail');
    (tabs[0] as HTMLElement).click();
    flush();
    expect(g.querySelector('[data-msui-view="inspect"]')).not.toBeNull();
    // A press on the open tab folds the group: its header stays, its body is paused.
    g.querySelector<HTMLElement>('[role="tab"][aria-selected="true"]')!.click();
    flush();
    expect(m.groupNode(lower.id)?.hasAttribute('data-folded')).toBe(true);
    expect(m.groupNode(lower.id)?.getAttribute('data-shown')).toBeNull();
    expect(b.probe.scope?.active).toBe(false);
    expect(a.probe.scope?.active).toBe(true);
    g.querySelector<HTMLElement>('[role="tab"][aria-selected="true"]')!.click();
    flush();
    expect(m.groupNode(lower.id)?.hasAttribute('data-folded')).toBe(false);
    expect(b.probe.scope?.active).toBe(true);
    m.dispose();
  });

  it('the rail hidden: the docked column and the rail go, the star in the main area brings them back', () => {
    const a = probeView('layers');
    const m = mountShell([a]);
    m.shell.show('layers');
    flush();
    expect(button(m.main(), 'Show the views')).toBeNull();
    m.shell.toggleRail(false);
    flush();
    expect(m.rail()).toBeNull();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    expect(a.probe.scope?.active).toBe(false);
    const star = button(m.main(), 'Show the views')!;
    expect(star).not.toBeNull();
    star.click();
    flush();
    expect(m.rail()).not.toBeNull();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBeNull();
    expect(a.probe.scope?.active).toBe(true);
    expect(button(m.main(), 'Show the views')).toBeNull();
    m.dispose();
  });
});

describe('ShellView: the docked tree', () => {
  it('resizes by its splitter (arrows, a drag) in place: the splitter and the main area stay, the drag goes on', () => {
    vi.spyOn(HTMLElement.prototype, 'getBoundingClientRect').mockImplementation(function (this: HTMLElement) {
      const w = this.matches('[data-node]')
        ? Number.parseFloat(this.style.getPropertyValue('--msui-shell-size')) || 0
        : 0;
      return { left: 0, top: 0, right: w, bottom: 0, width: w, height: 0, x: 0, y: 0, toJSON() {} } as DOMRect;
    });
    const m = mountShell([probeView('layers')], { size: 300 });
    m.shell.show('layers');
    flush();
    const main = m.main();
    const splitter = m.q('[role="separator"]')!;
    const size = () => Number.parseFloat(m.groupNode('right')!.style.getPropertyValue('--msui-shell-size'));
    expect(size()).toBe(300);
    // The right group: its splitter at its start (a drag left widens it).
    key(splitter, 'ArrowLeft');
    expect(size()).toBe(316);
    key(splitter, 'ArrowRight', { shiftKey: true });
    expect(size()).toBe(252);
    splitter.dispatchEvent(pointer('pointerdown', { x: 500, y: 0 }));
    splitter.dispatchEvent(pointer('pointermove', { x: 450, y: 0 }));
    flushFrame();
    flush();
    expect(size()).toBe(302);
    expect(m.q('[role="separator"]')).toBe(splitter);
    splitter.dispatchEvent(pointer('pointermove', { x: 400, y: 0 }));
    flushFrame();
    flush();
    expect(size()).toBe(352);
    splitter.dispatchEvent(pointer('pointerup', { x: 400, y: 0 }));
    expect(m.main()).toBe(main);
    m.dispose();
  });

  it('keeps the focus in a window when a press or the focus raises it (the tree is not rebuilt)', () => {
    const m = mountShell([probeView('layers'), probeView('inspect')]);
    m.shell.float('layers');
    m.shell.float('inspect');
    flush();
    const input = m.q<HTMLInputElement>('[data-input="layers"]')!;
    input.focus();
    flush();
    expect(document.activeElement).toBe(input);
    expect(Number(m.win('layers')!.style.zIndex)).toBeGreaterThan(Number(m.win('inspect')!.style.zIndex));
    m.dispose();
  });

  it('a group whose views all leave (their flags false) goes, and comes back', () => {
    const shown = new State(true);
    const m = mountShell([probeView('layers'), probeView('inspect'), probeView('kin', { visible: shown })]);
    m.shell.show('layers');
    m.shell.splitBelow('kin');
    flush();
    const lower = findView(m.shell.layout.value, 'kin')!.group.id;
    expect(m.groupNode(lower)?.getAttribute('data-shown')).toBeNull();
    shown.set(false);
    flush();
    expect(m.groupNode(lower)?.getAttribute('data-shown')).toBe('false');
    shown.set(true);
    flush();
    expect(m.groupNode(lower)?.getAttribute('data-shown')).toBeNull();
    m.dispose();
  });

  it('keeps the card open across layout changes', () => {
    const m = mountShell([probeView('layers'), probeView('export', { home: 'card' })]);
    m.shell.show('export');
    flush();
    expect(m.q('[data-card="export"]')).not.toBeNull();
    m.shell.float('layers');
    m.shell.toggleRail(false);
    m.shell.toggleRail(true);
    flush();
    expect(m.q('[data-card="export"]')).not.toBeNull();
    m.dispose();
  });
});

describe('ShellView: the rail’s height', () => {
  it('puts what does not fit under More, which shows the rest', () => {
    const ro = stubResizeObserver();
    const views = ['a', 'b', 'c', 'd', 'e', 'f'].map((id) => probeView(id));
    const m = mountShell(views);
    expect(m.railButton('More')).toBeNull();
    ro.resize(m.rail()!, { height: 120 });
    const more = m.railButton('More')!;
    expect(more).not.toBeNull();
    const hidden = m.qa('nav [role="tab"]').filter((b) => b.hidden);
    expect(hidden.length).toBeGreaterThan(0);
    more.click();
    flush();
    const label = hidden[0].getAttribute('aria-label')!;
    menuItem(label)!.click();
    flush();
    expect(m.frame(label.toLowerCase())).not.toBeNull();
    ro.resize(m.rail()!, { height: 1000 });
    expect(m.railButton('More')).toBeNull();
    expect(m.qa('nav [role="tab"]').filter((b) => b.hidden)).toHaveLength(0);
    m.dispose();
  });
});

describe('ShellView: windows', () => {
  const setup = () => {
    const a = probeView('layers');
    const b = probeView('inspect');
    const m = mountShell([a, b]);
    m.shell.float('layers');
    m.shell.float('inspect');
    flush();
    return { m, a, b };
  };

  it('float makes a window (a non-modal dialog) placed at its rect', () => {
    const { m } = setup();
    const w = m.win('layers')!;
    expect(w.getAttribute('role')).toBe('dialog');
    expect(w.getAttribute('aria-modal')).toBe('false');
    expect(w.getAttribute('aria-label')).toBe('Layers');
    const place = findView(m.shell.layout.value, 'layers')!.window!;
    expect(rectOf(w)).toEqual(place.rect);
    // Windows are in the main area.
    expect(m.main().contains(w)).toBe(true);
    m.dispose();
  });

  it('a press raises a window', () => {
    const { m } = setup();
    const z = (ref: string) => Number(m.win(ref)!.style.zIndex);
    expect(z('inspect')).toBeGreaterThan(z('layers'));
    m.q('[data-probe="layers"]')!.dispatchEvent(pointer('pointerdown', { x: 0, y: 0 }));
    flush();
    expect(z('layers')).toBeGreaterThan(z('inspect'));
    m.dispose();
  });

  it('Escape closes the window with the focus in it, not one elsewhere; ✕ closes; a press outside does not', () => {
    const { m } = setup();
    m.main().focus();
    key(m.main(), 'Escape');
    expect(m.win('layers')).not.toBeNull();
    expect(m.win('inspect')).not.toBeNull();
    const input = m.q<HTMLInputElement>('[data-input="layers"]')!;
    input.focus();
    key(input, 'Escape');
    expect(m.win('layers')).toBeNull();
    expect(m.win('inspect')).not.toBeNull();
    press(m.q('[data-app-main]')!);
    expect(m.win('inspect')).not.toBeNull();
    button(m.win('inspect')!, 'Close Inspect')!.click();
    flush();
    expect(m.win('inspect')).toBeNull();
    // A docked view closed from its window goes back home, hidden.
    expect(m.shell.layout.value.hidden).toContain('inspect');
    m.dispose();
  });

  it('Escape closes a menu opened in a window before the window', () => {
    const { m } = setup();
    // A real press focuses the button: the menu gives the focus back there.
    button(m.win('layers')!, 'Layers: layout')!.focus();
    button(m.win('layers')!, 'Layers: layout')!.click();
    flush();
    const item = menuItem('Dock')!;
    expect(item).toBeDefined();
    key(item, 'Escape');
    expect(menuItem('Dock')).toBeUndefined();
    expect(m.win('layers')).not.toBeNull();
    key(document.activeElement!, 'Escape');
    expect(m.win('layers')).toBeNull();
    m.dispose();
  });

  it('moves and resizes from the keyboard (the menu’s Move and Resize; arrows, Shift, Enter ends)', () => {
    const { m } = setup();
    const r0 = findView(m.shell.layout.value, 'layers')!.window!.rect;
    button(m.win('layers')!, 'Layers: layout')!.click();
    flush();
    menuItem('Move with the keyboard')!.click();
    flush();
    const w = m.win('layers')!;
    expect(w.hasAttribute('data-adjusting')).toBe(true);
    key(w, 'ArrowRight');
    key(w, 'ArrowDown', { shiftKey: true });
    expect(rectOf(w)).toEqual({ ...r0, x: r0.x + 4, y: r0.y + 16 });
    key(w, 'Enter');
    expect(w.hasAttribute('data-adjusting')).toBe(false);
    key(w, 'ArrowRight');
    expect(rectOf(w).x).toBe(r0.x + 4);
    button(w, 'Layers: layout')!.click();
    flush();
    menuItem('Resize with the keyboard')!.click();
    flush();
    key(w, 'ArrowLeft');
    key(w, 'ArrowUp');
    expect(rectOf(w)).toEqual({ x: r0.x + 4, y: r0.y + 16, width: r0.width - 4, height: r0.height - 4 });
    // Escape ends the adjusting, the window stays.
    key(w, 'Escape');
    expect(m.win('layers')).not.toBeNull();
    expect(w.hasAttribute('data-adjusting')).toBe(false);
    m.dispose();
  });

  it('drags by its header, resizes by its grips (in the frame lane)', () => {
    const { m } = setup();
    const r0 = findView(m.shell.layout.value, 'layers')!.window!.rect;
    const w = m.win('layers')!;
    const header = w.querySelector<HTMLElement>('header')!;
    header.dispatchEvent(pointer('pointerdown', { x: 100, y: 100 }));
    header.dispatchEvent(pointer('pointermove', { x: 130, y: 120 }));
    flush();
    expect(rectOf(w)).toEqual(r0);
    flushFrame();
    flush();
    expect(rectOf(w)).toEqual({ ...r0, x: r0.x + 30, y: r0.y + 20 });
    header.dispatchEvent(pointer('pointerup', { x: 140, y: 130 }));
    // The release's place goes through the frame lane (a commit beside it would be overwritten by its pending value).
    flushFrame();
    flush();
    expect(rectOf(w)).toEqual({ ...r0, x: r0.x + 40, y: r0.y + 30 });
    // A press on a header button is the button's.
    const menuButton = button(w, 'Layers: layout')!;
    menuButton.dispatchEvent(pointer('pointerdown', { x: 0, y: 0 }));
    menuButton.dispatchEvent(pointer('pointermove', { x: 50, y: 50 }));
    menuButton.dispatchEvent(pointer('pointerup', { x: 50, y: 50 }));
    flushFrame();
    flush();
    expect(rectOf(w)).toEqual({ ...r0, x: r0.x + 40, y: r0.y + 30 });
    const corner = w.querySelector<HTMLElement>('[data-edge="corner"]')!;
    corner.dispatchEvent(pointer('pointerdown', { x: 0, y: 0 }));
    corner.dispatchEvent(pointer('pointermove', { x: 10, y: 20 }));
    corner.dispatchEvent(pointer('pointerup', { x: 10, y: 20 }));
    flushFrame();
    flush();
    expect(rectOf(w)).toEqual({ x: r0.x + 40, y: r0.y + 30, width: r0.width + 10, height: r0.height + 20 });
    const right = w.querySelector<HTMLElement>('[data-edge="right"]')!;
    right.dispatchEvent(pointer('pointerdown', { x: 0, y: 0 }));
    right.dispatchEvent(pointer('pointerup', { x: 10, y: 50 }));
    flushFrame();
    flush();
    expect(rectOf(w)).toEqual({ x: r0.x + 40, y: r0.y + 30, width: r0.width + 20, height: r0.height + 20 });
    m.dispose();
  });

  it('two views in one window show as tabs; Dock brings a view home', () => {
    const { m } = setup();
    const target = findView(m.shell.layout.value, 'layers')!.group.id;
    m.shell.move('inspect', { kind: 'group', id: target });
    flush();
    const w = m.win('layers') ?? m.win('inspect')!;
    expect([...w.querySelectorAll('[role="tab"]')].map((t) => t.textContent)).toEqual(['Layers', 'Inspect']);
    button(w, 'Inspect: layout')!.click();
    flush();
    menuItem('Dock')!.click();
    flush();
    expect(m.win('inspect')).toBeNull();
    expect(m.groupNode('right')?.contains(m.frame('inspect'))).toBe(true);
    m.dispose();
  });
});

describe('ShellView: card, page, modal', () => {
  const setup = (extra: Parameters<typeof probeView>[1] = {}) => {
    const a = probeView('layers');
    const card = probeView('export', { home: 'card', ...extra });
    const page = probeView('help', { home: 'page' });
    const modal = probeView('confirm', { home: 'modal' });
    const m = mountShell([a, card, page, modal]);
    const anchor = (<button type="button">open</button>) as HTMLButtonElement;
    m.main().append(anchor);
    return { m, anchor, card };
  };

  it('a card opens in place, takes the focus, Escape closes it and the focus goes back to its anchor', async () => {
    const { m, anchor } = setup();
    anchor.focus();
    m.shell.show('export', { anchor });
    flush();
    const card = m.q('[data-card="export"]')!;
    expect(card).not.toBeNull();
    await microtasks();
    expect(card.contains(document.activeElement)).toBe(true);
    key(document.activeElement!, 'Escape');
    expect(m.q('[data-card]')).toBeNull();
    expect(document.activeElement).toBe(anchor);
    m.dispose();
  });

  it('a press outside closes a card (the main area too), not a persistent one; ✕ closes', () => {
    const { m } = setup();
    m.shell.show('export');
    flush();
    press(m.q('[data-probe="export"]')!);
    expect(m.q('[data-card]')).not.toBeNull();
    press(m.q('[data-app-main]')!);
    expect(m.q('[data-card]')).toBeNull();
    m.dispose();
    const p = setup({ persistent: true });
    p.m.shell.show('export');
    flush();
    press(p.m.q('[data-app-main]')!);
    expect(p.m.q('[data-card]')).not.toBeNull();
    button(p.m.q('[data-card]')!, 'Close Export')!.click();
    flush();
    expect(p.m.q('[data-card]')).toBeNull();
    p.m.dispose();
  });

  it('Escape with the focus in a window closes the window, not an open card', () => {
    const { m } = setup();
    m.shell.float('layers');
    m.shell.show('export');
    flush();
    const input = m.q<HTMLInputElement>('[data-input="layers"]')!;
    input.focus();
    key(input, 'Escape');
    expect(m.win('layers')).toBeNull();
    expect(m.q('[data-card]')).not.toBeNull();
    m.dispose();
  });

  it('a page covers the main area and Escape closes it', async () => {
    const { m, anchor } = setup();
    anchor.focus();
    m.shell.show('help', { anchor });
    flush();
    const page = m.q('[role="region"][aria-label="Help"]')!;
    expect(m.main().contains(page)).toBe(true);
    await microtasks();
    expect(document.activeElement).toBe(page);
    key(page, 'Escape');
    expect(m.q('[aria-label="Help"]')).toBeNull();
    expect(document.activeElement).toBe(anchor);
    m.dispose();
  });

  it('a modal makes the rest inert, keeps Tab in it, closes on Escape and lifts the inert', async () => {
    const { m, anchor } = setup();
    anchor.focus();
    m.shell.show('confirm', { anchor });
    flush();
    await microtasks();
    const dialog = m.q('[role="dialog"][aria-modal="true"]')!;
    expect(dialog.getAttribute('aria-label')).toBe('Confirm');
    // The shell's layout (everything but the layers that hold the app) is inert: decided as the modal opens, whatever is
    // in it then (a window opened under it too).
    const layout = m.q('[data-msui-shell-main]')!.closest('[class*="layout"]')!;
    expect(layout.hasAttribute('inert')).toBe(true);
    expect(layout.contains(dialog)).toBe(false);
    expect(dialog.contains(document.activeElement)).toBe(true);
    const all = [...dialog.querySelectorAll<HTMLElement>('button, input')];
    all.at(-1)!.focus();
    const tab = key(all.at(-1)!, 'Tab');
    expect(tab.defaultPrevented).toBe(true);
    expect(document.activeElement).toBe(all[0]);
    key(document.activeElement!, 'Escape');
    expect(m.q('[aria-modal="true"]')).toBeNull();
    expect(layout.hasAttribute('inert')).toBe(false);
    expect(document.activeElement).toBe(anchor);
    m.dispose();
  });
});

describe('ShellView: compact', () => {
  it('follows the shell’s width (a ResizeObserver): the sheet shows the docked view, the card, the top window', () => {
    const ro = stubResizeObserver();
    const a = probeView('layers', { compactBar: true });
    const b = probeView('inspect', { compactBar: true });
    const c = probeView('settings');
    const card = probeView('export', { home: 'card' });
    const m = mountShell([a, b, c, card]);
    m.shell.show('layers');
    flush();
    const input = m.q<HTMLInputElement>('[data-input="layers"]')!;
    input.value = 'kept';
    const shellEl = m.q('[role="application"]')!;
    ro.resize(shellEl, { width: 500 });
    expect(m.shell.mode.value).toBe('compact');
    expect(shellEl.getAttribute('data-mode')).toBe('compact');
    expect(m.rail()).toBeNull();
    // The main area is kept (the same element).
    expect(m.q('[data-app-main]')).not.toBeNull();
    // The bar: the tab list of views, then More beside it (not a tab).
    const sheet = m.q('[role="tablist"][aria-label="Views"]')!.parentElement!;
    const bar = [...sheet.querySelectorAll('button')].map((x) => x.textContent);
    expect(bar).toEqual(['Layers', 'Inspect', 'More']);
    expect(m.frame('layers')?.getAttribute('data-placement')).toBe('sheet');
    expect(m.q<HTMLInputElement>('[data-input="layers"]')).toBe(input);
    expect(input.value).toBe('kept');
    expect(sheet.querySelector('[aria-selected="true"]')?.textContent).toBe('Layers');
    (sheet.querySelector('button:nth-child(2)') as HTMLElement).click();
    flush();
    expect(m.frame('inspect')).not.toBeNull();
    // More lists the rest.
    (sheet.querySelector('[aria-haspopup="menu"]') as HTMLElement).click();
    flush();
    menuItem('Settings')!.click();
    flush();
    expect(m.frame('settings')).not.toBeNull();
    // A window shows in the sheet, over the docked view.
    m.shell.float('layers');
    flush();
    expect(m.frame('layers')?.getAttribute('data-placement')).toBe('sheet');
    expect(m.q('[data-window]')).toBeNull();
    // A card over both.
    m.shell.show('export');
    flush();
    expect(m.frame('export')).not.toBeNull();
    button(m.frame('export')!, 'Close Export')!.click();
    flush();
    expect(m.frame('layers')).not.toBeNull();
    // Wide again: the window is back, with the same body.
    ro.resize(shellEl, { width: 1200 });
    expect(m.shell.mode.value).toBe('desktop');
    expect(m.win('layers')).not.toBeNull();
    expect(m.q<HTMLInputElement>('[data-input="layers"]')).toBe(input);
    expect(input.value).toBe('kept');
    m.dispose();
  });

  it('starts compact by its option', () => {
    const m = mountShell([probeView('layers', { compactBar: true })], { mode: 'compact' });
    expect(m.q('[role="tablist"][aria-label="Views"]')).not.toBeNull();
    // No docked tree beside the main area (the main area's own region is in it, in every mode).
    expect([...m.root.querySelectorAll('[data-node]')].filter((n) => !n.closest('[data-msui-shell-main]'))).toEqual([]);
    m.dispose();
  });
});

describe('ShellView: F6', () => {
  it('cycles the focus between the main area, the groups, the rail and the windows', () => {
    // jsdom lays nothing out: every connected element has a box.
    vi.spyOn(Element.prototype, 'getClientRects').mockImplementation(function (this: Element) {
      return (this.isConnected ? [{}] : []) as unknown as DOMRectList;
    });
    const m = mountShell([probeView('layers'), probeView('inspect')]);
    m.shell.show('layers');
    m.shell.float('inspect');
    flush();
    const regionOf = (el: Element | null) => {
      if (!el) return 'none';
      if (el.closest('[data-window]')) return 'window';
      if (el.closest('[data-group]')) return 'group';
      if (el.closest('nav')) return 'rail';
      if (el.closest('[data-msui-shell-main]')) return 'main';
      return 'other';
    };
    m.main().focus();
    const seen: string[] = [];
    for (let i = 0; i < 4; i++) {
      key(document.activeElement!, 'F6');
      seen.push(regionOf(document.activeElement));
    }
    expect(new Set(seen)).toEqual(new Set(['window', 'group', 'rail', 'main']));
    expect(seen.at(-1)).toBe('main');
    const back: string[] = [];
    for (let i = 0; i < 4; i++) {
      key(document.activeElement!, 'F6', { shiftKey: true });
      back.push(regionOf(document.activeElement));
    }
    expect(back).toEqual([...seen.slice(0, 3).reverse(), 'main']);
    m.dispose();
  });
});
