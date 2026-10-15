/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Menus and keyboard focus in the vanilla UI (UX2; UI2's cases; seqspec:interaction#slots-and-chrome,
// seqspec:interaction#accessibility): a menu opens on its trigger or on the keyboard cursor's Shift+F10, focuses its
// first item, moves with the arrows, closes on Escape, on a choice and on an outside press and gives focus back to its
// trigger; a host's Menu gets `open` through `update`; a node a part rebuilds keeps focus by `data-seqv-key`, and a
// node that goes lands focus on the next of its group, the previous, the section's surface.

import { afterEach, describe, expect, it, vi } from 'vitest';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import type { SeqViewerContext } from '../../src/index.js';
import { at, flush, keys, openedContext, pointer, stubSurface, WIDTH, xOf } from './fixtures.js';

const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
});

async function mounted(spec: D.DomUISpec = {}) {
  const t = await openedContext();
  cleanups.push(() => t.ctx.dispose());
  const host = document.createElement('div');
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  cleanups.push(() => host.remove());
  const m = D.mountSeqViewer(host, t.ctx, { surface: stubSurface, ...spec });
  cleanups.unshift(() => m.dispose());
  await flush();
  return { ...t, host };
}
const S = 'chains.main';
const active = () => document.activeElement as HTMLElement | null;
const press = (el: Element, key: string, init: KeyboardEventInit = {}) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key, bubbles: true, cancelable: true, composed: true, ...init }));
const surfaceOf = (host: HTMLElement) => host.querySelector('.seqv-surface [tabindex]') as HTMLElement;
const trigger = (host: HTMLElement) => host.querySelector('.seqv-track-actions button') as HTMLButtonElement;
const items = (host: HTMLElement) => [...host.querySelectorAll<HTMLButtonElement>('[role^=menuitem]')];

/** The menu at the Domains header, opened by its trigger. */
async function menuAtDomains(spec: D.DomUISpec = {}) {
  const t = await mounted(spec);
  const h = at(t.ctx, S, 't.dom');
  pointer(t.ctx)('move', S, h.header, h.y);
  return t;
}
/** The keyboard cursor on the Domains track: the surface focused, Down until it is there. */
function cursorOnDomains(ctx: SeqViewerContext) {
  const key = keys(ctx);
  ctx.input.focus(ctx.layout.sectionIndex(S), true);
  for (let i = 0; i < 8 && ctx.ui.keyCursor.value?.track !== 't.dom'; i++) key(S, 'ArrowDown');
  expect(ctx.ui.keyCursor.value?.track).toBe('t.dom');
  return key;
}

describe('the menu', () => {
  it('opens on its trigger with the first item focused, the arrows move between items, Home and End go to the ends', async () => {
    const { host } = await menuAtDomains();
    const t = trigger(host);
    expect([t.getAttribute('aria-haspopup'), t.getAttribute('aria-expanded')]).toEqual(['menu', 'false']);
    t.click();
    const all = items(host);
    expect(all.length).toBe(4);
    expect([t.getAttribute('aria-expanded'), t.getAttribute('aria-controls')]).toEqual([
      'true',
      host.querySelector('[role=menu]')?.id,
    ]);
    expect(active()).toBe(all[0]);
    press(all[0], 'ArrowDown');
    expect(active()).toBe(all[1]);
    press(all[1], 'ArrowDown');
    press(all[2], 'ArrowDown');
    expect(active()).toBe(all[3]);
    press(all[3], 'ArrowDown');
    expect(active()).toBe(all[0]);
    press(all[0], 'ArrowUp');
    expect(active()).toBe(all[3]);
    press(all[3], 'Home');
    expect(active()).toBe(all[0]);
    press(all[0], 'End');
    expect(active()).toBe(all[3]);
    // Closed, ArrowDown on the focused trigger opens it at the first item.
    press(active() as HTMLElement, 'Escape');
    expect(active()).toBe(t);
    press(t, 'ArrowDown');
    expect([t.getAttribute('aria-expanded'), active()]).toEqual(['true', items(host)[0]]);
  });

  it('Escape closes it and focus goes back to the trigger; so does choosing an item, which still runs', async () => {
    const { host, ctx } = await menuAtDomains();
    const t = trigger(host);
    t.click();
    expect(press(active() as HTMLElement, 'Escape')).toBe(false);
    expect([host.querySelector('[role=menu]'), active(), t.getAttribute('aria-expanded')]).toEqual([null, t, 'false']);
    // A choice that leaves the trigger in place: "Focus" moves the viewport.
    t.click();
    const focus = items(host)[1];
    expect(focus.textContent).toBe('Focus');
    focus.focus();
    focus.click();
    await flush();
    expect([host.querySelector('[role=menu]'), active()]).toEqual([null, t]);
    expect(ctx.focus.current.value.size).toBeGreaterThan(0);
  });

  it("an open menu reads its items again on the streams' version: a failed region adds Retry, focus stays on its item", async () => {
    const { host, ctx } = await menuAtDomains();
    trigger(host).click();
    const labels = () => items(host).map((i) => i.textContent);
    expect(labels()).not.toContain('Retry');
    const first = active()?.dataset.seqvKey;
    vi.spyOn(ctx.streams, 'hasFailed').mockReturnValue(true);
    (ctx.streams.version as unknown as { set(v: number): void }).set(ctx.streams.version.value + 1);
    expect(labels().at(-1)).toBe('Retry');
    expect(active()?.dataset.seqvKey).toBe(first);
  });

  it('Tab closes it without trapping focus', async () => {
    const { host } = await menuAtDomains();
    const t = trigger(host);
    t.click();
    press(active() as HTMLElement, 'Tab');
    expect([host.querySelector('[role=menu]'), active()]).toEqual([null, t]);
  });

  it('focus moving out of it closes it and stays where it went (an item rebuilt under focus does not: Retry above)', async () => {
    const { host } = await menuAtDomains();
    trigger(host).click();
    expect(host.querySelector('[role=menu]')).not.toBeNull();
    const surface = surfaceOf(host);
    surface.focus();
    expect([host.querySelector('[role=menu]'), active()]).toEqual([null, surface]);
  });

  it('a press outside closes it (a capture-phase pointerdown, which Safari needs), a press inside does not', async () => {
    const { host } = await menuAtDomains();
    const t = trigger(host);
    t.click();
    items(host)[0].dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));
    expect(host.querySelector('[role=menu]')).not.toBeNull();
    expect(active()).toBe(items(host)[0]);
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect([host.querySelector('[role=menu]'), t.getAttribute('aria-expanded')]).toEqual([null, 'false']);
    // Focus is the press's to move, not taken back to the trigger.
    expect(active()).not.toBe(t);
    // And once closed, a press outside asks nothing more of it.
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true }));
    expect(host.querySelector('[role=menu]')).toBeNull();
  });

  it("Shift+F10 and ContextMenu on the keyboard cursor's track open its menu, the first item focused", async () => {
    const { host, ctx } = await mounted();
    const key = cursorOnDomains(ctx);
    expect(host.querySelector('[role=menu]')).toBeNull();
    // F10 alone is not the key.
    expect(key(S, 'F10')).toBe(false);
    expect(host.querySelector('[role=menu]')).toBeNull();
    expect(key(S, 'F10', true)).toBe(true);
    expect(items(host).map((i) => i.textContent)).toEqual([
      'Show only what this track covers',
      'Focus',
      'Colors ▸',
      'Hide track',
    ]);
    expect(active()).toBe(items(host)[0]);
    press(active() as HTMLElement, 'Escape');
    expect([host.querySelector('[role=menu]'), active()?.getAttribute('aria-label')]).toEqual([
      null,
      'Actions for Domains',
    ]);
    expect(key(S, 'ContextMenu')).toBe(true);
    expect(host.querySelector('[role=menu]')).not.toBeNull();
  });

  it("a host's Menu receives `open` through `update` on Shift+F10, and its `onOpenChange` closes it", async () => {
    const seen: boolean[] = [];
    let last: D.MenuProps | null = null;
    const HostMenu: D.DomPart<D.MenuProps> = (host, _ctx, p) => {
      const root = D.h('span', { class: 'host-menu' }, p.trigger);
      host.append(root);
      [last] = [p];
      seen.push(p.open);
      return {
        update(next) {
          last = next;
          seen.push(next.open);
        },
        dispose: () => root.remove(),
      };
    };
    const { host, ctx } = await mounted({ slots: { Menu: HostMenu } });
    const key = cursorOnDomains(ctx);
    expect(host.querySelector('.host-menu')).not.toBeNull();
    expect(seen.at(-1)).toBe(false);
    key(S, 'F10', true);
    expect(seen.at(-1)).toBe(true);
    (last as D.MenuProps | null)?.onOpenChange(false);
    expect(seen.at(-1)).toBe(false);
  });

  it('an open menu keeps its box while the pointer is elsewhere and the cursor is gone, and focus stays in it', async () => {
    const { host, ctx } = await menuAtDomains();
    trigger(host).click();
    const first = items(host)[0];
    expect(active()).toBe(first);
    ctx.input.leave(ctx.layout.sectionIndex(S));
    await flush();
    const res = at(ctx, S, 't.res');
    pointer(ctx)('move', S, xOf(ctx, S, 12), res.y);
    await flush();
    expect([host.querySelector('[role=menu]'), active()]).toEqual([expect.anything(), first]);
    // Closed and focus moved off the box: the box goes, and focus lands on the surface.
    press(first, 'Escape');
    surfaceOf(host).focus();
    ctx.input.leave(ctx.layout.sectionIndex(S));
    await flush();
    expect(host.querySelector('.seqv-track-actions')).toBeNull();
    expect(active()).toBe(surfaceOf(host));
  });
});

describe("the track's box", () => {
  it('stays while its closed trigger holds focus and the pointer is elsewhere, and goes when focus leaves it', async () => {
    const { host, ctx } = await menuAtDomains();
    const t = trigger(host);
    t.click();
    press(active() as HTMLElement, 'Escape');
    expect(active()).toBe(t);
    // Nothing hovered (the pointer on the data area, no keyboard cursor): focus inside keeps it.
    pointer(ctx)('move', S, xOf(ctx, S, 12), at(ctx, S, 't.res').y);
    await flush();
    expect([trigger(host), active()]).toEqual([t, t]);
    surfaceOf(host).focus();
    ctx.input.leave(ctx.layout.sectionIndex(S));
    await flush();
    expect(host.querySelector('.seqv-track-actions')).toBeNull();
  });

  it("stays while a host's Menu is open (its trigger `aria-expanded`), with focus outside the box", async () => {
    let last: D.MenuProps | null = null;
    const HostMenu: D.DomPart<D.MenuProps> = (host, _ctx, p) => {
      host.append(p.trigger);
      const show = (q: D.MenuProps) => {
        last = q;
        D.setAttr(q.trigger, 'aria-expanded', q.open);
      };
      show(p);
      return { update: show, dispose: () => p.trigger.remove() };
    };
    const { host, ctx } = await menuAtDomains({ slots: { Menu: HostMenu } });
    (last as D.MenuProps | null)?.onOpenChange(true);
    expect(trigger(host).getAttribute('aria-expanded')).toBe('true');
    ctx.input.leave(ctx.layout.sectionIndex(S));
    await flush();
    expect(host.querySelector('.seqv-track-actions')).not.toBeNull();
    (last as D.MenuProps | null)?.onOpenChange(false);
    ctx.input.leave(ctx.layout.sectionIndex(S));
    pointer(ctx)('move', S, xOf(ctx, S, 12), at(ctx, S, 't.res').y);
    await flush();
    expect(host.querySelector('.seqv-track-actions')).toBeNull();
  });
});

describe('focus kept by key, and where it lands', () => {
  it('a switch pressed keeps focus on its button when the header is rebuilt, and so does a hover over other tracks', async () => {
    const { host, ctx } = await mounted();
    const sel = '[data-seqv-key="choice:section.rowSizing:natural"]';
    const before = host.querySelector(sel) as HTMLButtonElement;
    expect(before).not.toBeNull();
    before.focus();
    before.click();
    await flush();
    const after = host.querySelector(sel) as HTMLButtonElement;
    expect([after === before, after.getAttribute('aria-pressed'), active()]).toEqual([false, 'true', after]);
    for (const track of ['t.res', 't.dom', 't.res']) {
      const h = at(ctx, S, track);
      pointer(ctx)('move', S, xOf(ctx, S, 12), h.y);
    }
    expect(active()).toBe(after);
  });

  it("a switch with a menu in the header opens and closes like the track's; a choice rebuilds the header and focus stays on the trigger", async () => {
    const { host, ctx } = await mounted();
    const sel = '[data-seqv-key="switch:section.style"]';
    const t = host.querySelector(sel) as HTMLButtonElement;
    expect(t.getAttribute('aria-haspopup')).toBe('menu');
    t.click();
    press(active() as HTMLElement, 'Escape');
    expect([host.querySelector('[role=menu]'), active()]).toEqual([null, t]);
    t.click();
    const preset = items(host)[1];
    preset.focus();
    preset.click();
    await flush();
    expect(ctx.ui.style.value.get(S)).toBe(preset.textContent);
    const after = host.querySelector(sel) as HTMLButtonElement;
    expect([after === t, active(), host.querySelector('[role=menu]')]).toEqual([false, after, null]);
  });

  it('"Show X" moves focus to the next chip, then, with none left, to the surface', async () => {
    const { host, ctx } = await mounted();
    ctx.actions.run('track.hide', 't.res');
    ctx.actions.run('track.hide', 't.dom');
    await flush();
    const chips = () => [...host.querySelectorAll<HTMLButtonElement>('.seqv-chips button')];
    expect(chips().map((c) => c.textContent)).toEqual(['Show Residues', 'Show Domains']);
    chips()[0].focus();
    chips()[0].click();
    await flush();
    expect(chips().map((c) => c.textContent)).toEqual(['Show Domains']);
    expect(active()).toBe(chips()[0]);
    chips()[0].click();
    await flush();
    expect(host.querySelector('.seqv-chips')).toBeNull();
    expect(active()).toBe(surfaceOf(host));
  });

  it("a last chip's neighbor before it takes focus when none follows", async () => {
    const { host, ctx } = await mounted();
    ctx.actions.run('track.hide', 't.res');
    ctx.actions.run('track.hide', 't.dom');
    await flush();
    const chips = () => [...host.querySelectorAll<HTMLButtonElement>('.seqv-chips button')];
    chips()[1].focus();
    chips()[1].click();
    await flush();
    expect([chips().map((c) => c.textContent), active()]).toEqual([['Show Residues'], chips()[0]]);
  });

  it("Hide from the track's menu leaves focus on the section's surface, not on the page", async () => {
    const { host } = await menuAtDomains();
    trigger(host).click();
    const hide = items(host)[3];
    expect(hide.textContent).toBe('Hide track');
    hide.focus();
    hide.click();
    await flush();
    expect(host.querySelector('.seqv-track-actions')).toBeNull();
    expect(host.querySelector('.seqv-chips')?.textContent).toBe('Show Domains');
    expect(active()).toBe(surfaceOf(host));
  });

  it("the player's Next at the last snapshot and Previous at the first stay focusable, say so, and do nothing", async () => {
    const { host, ctx } = await mounted();
    const next = host.querySelector('[aria-label="Next snapshot"]') as HTMLButtonElement;
    const prev = host.querySelector('[aria-label="Previous snapshot"]') as HTMLButtonElement;
    expect([prev.getAttribute('aria-disabled'), prev.hasAttribute('disabled')]).toEqual(['true', false]);
    prev.focus();
    prev.click();
    expect([active(), ctx.snapshots.current.value]).toEqual([prev, null]);
    next.focus();
    next.click();
    await flush();
    expect([ctx.snapshots.current.value, active()]).toEqual(['s1', next]);
    next.click();
    await flush();
    expect([ctx.snapshots.current.value, next.getAttribute('aria-disabled'), next.hasAttribute('disabled')]).toEqual([
      's2',
      'true',
      false,
    ]);
    next.click();
    await flush();
    expect([ctx.snapshots.current.value, active()]).toEqual(['s2', next]);
  });

  it("focus inside a section that leaves the view lands on the viewer's root", async () => {
    const { host, ctx } = await mounted();
    surfaceOf(host).focus();
    expect(active()).toBe(surfaceOf(host));
    await ctx.actions.run('view.switch', 'msa');
    await flush();
    // The chains section's surface is gone with its frame; the alignment's has not taken focus.
    expect(active()).toBe(host.querySelector('.seqv-root'));
  });

  it("`spot(host).set` is for a host's parts too: the same key stays focused, a gone one lands on its group, then the root", () => {
    const host = D.h('div', { class: 'seqv-root', tabindex: -1 });
    document.body.append(host);
    cleanups.push(() => host.remove());
    const at = D.spot(host);
    const row = (...names: string[]) =>
      D.h('div', null, ...names.map((n) => D.button({ 'data-seqv-key': n, 'data-seqv-group': 'g' }, () => {}, n)));
    const one = row('a', 'b', 'c');
    at.set(one);
    (one.children[1] as HTMLElement).focus();
    at.set(row('a', 'b', 'c', 'd'));
    expect(active()?.dataset.seqvKey).toBe('b');
    at.set(row('a', 'c'));
    expect(active()?.dataset.seqvKey).toBe('c');
    (host.querySelector('[data-seqv-key=c]') as HTMLElement).focus();
    at.set(row('a'));
    expect(active()?.dataset.seqvKey).toBe('a');
    at.set(null);
    expect(active()).toBe(host);
    // A part disposed with focus inside lands it too.
    at.set(row('a'));
    (host.querySelector('[data-seqv-key=a]') as HTMLElement).focus();
    at.dispose();
    expect([host.querySelector('[data-seqv-key=a]'), active()]).toEqual([null, host]);
  });
});
