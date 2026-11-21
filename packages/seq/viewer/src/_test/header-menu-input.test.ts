// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An overlay in the mounted canvas owns its events (_spec/README.md "Input", Overlays): the header menu opened over the
// view, driven by the event sequences a browser sends (pointermove, pointerdown, pointerup, click; the container's
// pointer capture retargeting the release and the click to it, as browsers do). Each item applies its command, a click
// on an item selects nothing under it, the pointer over the menu leaves the view (no hover, no card), keys in the menu
// move no cursor, and `Color ▸` opens its submenu by click, hover and ArrowRight.

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type HeadersBehavior, HeadersProvider } from '../behaviors/headers.js';
import type { HoverCardBehavior } from '../behaviors/hover-card.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import type { HeadersParams } from '../header-entries.js';
import { headerMenuOverlay } from '../mount-header-menu.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/protein';
const Letters = `${Section}/letters`;
const Domains = `${Section}/domains`;

function spec() {
  const b = createSeqSpec({ title: 'header menu input' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters', { label: 'Letters' }).layer('l', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames;
let timers: Map<number, () => void>;
let timerId = 0;
const viewers: SeqViewer[] = [];
beforeEach(() => {
  frames = installFrames();
  timers = new Map();
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  document.body.textContent = '';
});

/** The rests (the card's, the announcer's, a submenu's) on a fake clock from here on: run when the test says. */
const fakeTimers = () => {
  vi.spyOn(DomHost, 'setTimeout').mockImplementation((fn) => {
    timers.set(++timerId, fn);
    return timerId;
  });
  vi.spyOn(DomHost, 'clearTimeout').mockImplementation((id) => {
    timers.delete(id as number);
  });
};
const runTimers = () => {
  const due = [...timers.values()];
  timers.clear();
  for (const fn of due) fn();
};

async function mount(params: HeadersParams = {}) {
  const target = document.createElement('div');
  document.body.appendChild(target);
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer(
    {
      ...base,
      ui: { canvas: true, history: uiHistory },
      config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
      behaviors: base.behaviors.map((b) =>
        b.provider === HeadersProvider
          ? { provider: HeadersProvider, params: { menu: headerMenuOverlay(), ...params } }
          : b,
      ),
    },
    target,
  );
  viewers.push(v);
  await v.load(spec());
  await v.settle();
  frames.flush();
  const container = target.querySelector<HTMLElement>('[role="application"]')!;
  // pointer capture as a browser keeps it: the release and the click after it go to the capturing element
  let captured: Element | undefined;
  container.setPointerCapture = () => {
    captured = container;
  };
  container.hasPointerCapture = () => captured === container;
  container.releasePointerCapture = () => {
    captured = undefined;
  };
  const at = (el: Element) => {
    const r = el.getBoundingClientRect();
    return { clientX: r.left + 1, clientY: r.top + 1 };
  };
  const pointer = (type: string, el: Element, xy = at(el)) =>
    el.dispatchEvent(
      new PointerEvent(type, { ...xy, bubbles: true, cancelable: true, button: 0, pointerId: 1, pointerType: 'mouse' }),
    );
  /** A click as the browser sends it: move, down, up and click (the last two at the capturing element if any). */
  const press = (el: Element, xy = at(el)) => {
    pointer('pointermove', el, xy);
    pointer('pointerdown', el, xy);
    const t = captured ?? el;
    pointer('pointerup', t, xy);
    captured = undefined;
    t.dispatchEvent(new MouseEvent('click', { ...xy, bubbles: true, cancelable: true, button: 0, detail: 1 }));
  };
  const headers = v.behaviors.get('headers') as HeadersBehavior;
  const card = v.behaviors.get('hoverCard') as HoverCardBehavior;
  return { v, container, headers, card, pointer, press };
}

const geo = (v: SeqViewer, track: string) => v.headers!.value.find((g) => g.target.track === track)!;
const tracks = (v: SeqViewer) => v.view.resolved.value!.sections[0].tracks.map((t) => t.address);

/** The menu of a track's header, opened at its label (as its "⋯" opens it). */
async function openMenu(v: SeqViewer, headers: HeadersBehavior, container: HTMLElement, track: string) {
  const g = geo(v, track);
  expect(headers.openMenu(g.target, { x: g.label.x, y: g.line.y + g.line.height })).toBe(true);
  await v.settle();
  frames.flush();
  const menu = container.querySelector<HTMLElement>('[role="menu"]')!;
  expect(menu).not.toBeNull();
  return menu;
}
const row = (menu: HTMLElement, label: RegExp) =>
  [...menu.querySelectorAll<HTMLElement>('[role^="menuitem"]')].find((b) => label.test(b.textContent ?? ''))!;
/** An item placed over the data (a menu opened over the view), at (x, y) of the canvas. */
const over = (el: HTMLElement, x: number, y: number) => {
  const r = el.closest('[role="application"]')!.getBoundingClientRect();
  vi.spyOn(el, 'getBoundingClientRect').mockReturnValue(new DOMRect(r.left + x, r.top + y, 100, 20));
  return el;
};

describe('the header menu over the view', () => {
  it('a click on a "⋯" with no move before it (the pointer elsewhere) opens its menu', async () => {
    const { v, container, headers, pointer } = await mount();
    const r = container.getBoundingClientRect();
    // the pointer on another track's line: its actions show, not these
    const g = geo(v, Letters);
    pointer('pointermove', container, {
      clientX: r.left + g.label.x + 2,
      clientY: r.top + g.line.y + g.line.height / 2,
    });
    frames.flush();
    const dots = geo(v, Domains).slots.find((s) => s.name === 'seq:menu')!.box;
    const xy = { clientX: r.left + dots.x + dots.width / 2, clientY: r.top + dots.y + dots.height / 2 };
    const canvas = container.querySelector('canvas')!;
    pointer('pointerdown', canvas, xy);
    const t = container.hasPointerCapture(1) ? container : canvas;
    pointer('pointerup', t, xy);
    t.dispatchEvent(new MouseEvent('click', { ...xy, bubbles: true, cancelable: true, button: 0, detail: 1 }));
    await v.settle();
    expect(headers.menu.value?.menu.target.track).toBe(Domains);
  });

  it('the pointer onto a header, then a click on its "⋯" as a browser clicks: the menu opens', async () => {
    const { v, container, headers, pointer, press } = await mount();
    const r = container.getBoundingClientRect();
    const g = geo(v, Domains);
    pointer('pointermove', container, {
      clientX: r.left + g.label.x + 2,
      clientY: r.top + g.line.y + g.line.height / 2,
    });
    frames.flush();
    const dots = geo(v, Domains).slots.find((s) => s.name === 'seq:menu')!.box;
    press(container.querySelector('canvas')!, {
      clientX: r.left + dots.x + dots.width / 2,
      clientY: r.top + dots.y + dots.height / 2,
    });
    await v.settle();
    expect(headers.menu.value?.menu.target.track).toBe(Domains);
    expect(container.querySelector('[role="menu"]')).not.toBeNull();
  });

  it('each item applies its command when clicked as a browser clicks, over the data: Move up, Reverse, Hide', async () => {
    const { v, container, headers, press } = await mount();
    // the menu's rows over the data, where a press would begin the view's drag (and capture the pointer)
    const item = async (label: RegExp) =>
      over(row(await openMenu(v, headers, container, Domains), label), v.dataX + 12, geo(v, Letters).line.y + 2);
    press(await item(/^Move up/));
    await v.settle();
    expect(tracks(v)).toEqual([Domains, Letters]);
    expect(container.querySelector('[role="menu"]')).toBeNull();
    const direction = () => {
      let n = v.spec.root.value;
      for (const id of Section.split('/').filter(Boolean)) n = n.children!.find((c) => c.id === id)!;
      return n.params?.direction;
    };
    press(await item(/^Reverse/));
    await v.settle();
    expect(direction()).toBe('reversed');
    press(await item(/^Hide track/));
    await v.settle();
    expect(tracks(v)).toEqual([Letters]);
  });

  it('a click on an item selects nothing in the view under it', async () => {
    const { v, container, headers, press } = await mount();
    const g = geo(v, Domains);
    const menu = await openMenu(v, headers, container, Letters);
    press(over(row(menu, /^Reverse/), v.dataX + 12, g.line.y + g.line.height / 2));
    await v.settle();
    expect(v.selection.value.value?.location.parts ?? []).toEqual([]);
    expect(v.selection.cursor.value).toBeUndefined();
  });

  it('the pointer over the menu leaves the view: its hover cleared, no card', async () => {
    const { v, container, headers, card, pointer } = await mount();
    const g = geo(v, Domains);
    const r = container.getBoundingClientRect();
    const menu = await openMenu(v, headers, container, Letters);
    fakeTimers();
    // on a residue: hovered
    pointer('pointermove', container, {
      clientX: r.left + v.dataX + 12,
      clientY: r.top + g.line.y + g.line.height / 2,
    });
    frames.flush();
    expect(v.hover.current.section).toBeDefined();
    // onto the menu, over the same residue
    pointer('pointermove', over(row(menu, /^Reverse/), v.dataX + 12, g.line.y + 2));
    frames.flush();
    runTimers();
    expect(v.hover.current.section).toBeUndefined();
    expect(card.shown.value).toBeUndefined();
  });

  it('keys in the menu are the menu’s: the arrows move no cursor in the view', async () => {
    const { v, container, headers } = await mount();
    const menu = await openMenu(v, headers, container, Letters);
    const before = v.selection.cursor.value;
    (document.activeElement as HTMLElement).dispatchEvent(
      new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, cancelable: true }),
    );
    expect(menu.contains(document.activeElement)).toBe(true);
    expect(v.selection.cursor.value).toBe(before);
  });
});

describe('an open menu and the view resolving again', () => {
  it('a load, or an edit elsewhere, keeps it; its track gone closes it', async () => {
    const { v, container, headers } = await mount();
    await openMenu(v, headers, container, Domains);
    // the same view loaded again (a source's answer, a reload): the menu's track is still there
    await v.load(spec());
    await v.settle();
    frames.flush();
    expect(headers.menu.value?.menu.target.track).toBe(Domains);
    expect(container.querySelector('[role="menu"]')).not.toBeNull();
    // another track hidden: still there
    await v.execute('track.hide', { address: Letters });
    await v.settle();
    frames.flush();
    expect(container.querySelector('[role="menu"]')).not.toBeNull();
    // its own track hidden: closed
    await v.execute('track.hide', { address: Domains });
    await v.settle();
    frames.flush();
    expect(headers.menu.value).toBeUndefined();
    expect(container.querySelector('[role="menu"]')).toBeNull();
  });
});

/** The open menu levels (the menu, then its flyouts), their rows' texts. */
const levels = (container: HTMLElement) =>
  [...container.querySelectorAll('[role="menu"]')].map((m) =>
    [...m.querySelectorAll('[role^="menuitem"]')].map((b) => b.textContent),
  );
const menusOf = (container: HTMLElement) => [...container.querySelectorAll<HTMLElement>('[role="menu"]')];
const key = (el: Element, k: string) =>
  el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));

describe('Color ▸ (a flyout beside its row: the menu stays)', () => {
  it('a click opens it beside the row, the parent still open and its row marked; a hovered scheme previews, a choice applies and closes both', async () => {
    const { v, container, headers, pointer } = await mount();
    const menu = await openMenu(v, headers, container, Letters);
    const color = row(menu, /^Color/);
    expect(color.getAttribute('aria-expanded')).toBe('false');
    color.click();
    await v.settle();
    frames.flush();
    expect(menusOf(container)).toHaveLength(2);
    // the parent is still there, whole, its Color row marked open
    expect(menu.isConnected).toBe(true);
    expect(levels(container)[0]).toContain('Hide track');
    expect(color.getAttribute('aria-expanded')).toBe('true');
    const flyout = menusOf(container)[1];
    expect(flyout.querySelectorAll('[role="menuitemradio"]').length).toBeGreaterThan(1);
    expect(headers.menu.value).toBeDefined();
    const choices = [...flyout.querySelectorAll<HTMLElement>('[role="menuitemradio"]')];
    const other = choices.find((c) => c.getAttribute('aria-checked') !== 'true')!;
    const name = other.textContent;
    pointer('pointerover', other);
    frames.flush();
    await v.settle();
    expect(v.spec.preview.value).toBeDefined();
    expect(menusOf(container)).toHaveLength(2);
    other.click();
    await v.settle();
    expect(menusOf(container)).toHaveLength(0);
    expect(headers.menu.value).toBeUndefined();
    const again = row(await openMenu(v, headers, container, Letters), /^Color/);
    again.click();
    await v.settle();
    expect(menusOf(container)[1].querySelector('[aria-checked="true"]')?.textContent).toBe(name);
  });

  it('is placed right of its row, left of the menu where there is no room, inside the window', async () => {
    const { v, container, headers } = await mount();
    const width = vi.spyOn(HTMLElement.prototype, 'offsetWidth', 'get').mockReturnValue(160);
    const height = vi.spyOn(HTMLElement.prototype, 'offsetHeight', 'get').mockReturnValue(100);
    try {
      const at = async (left: number, top: number) => {
        const menu = await openMenu(v, headers, container, Letters);
        const color = row(menu, /^Color/);
        vi.spyOn(menu, 'getBoundingClientRect').mockReturnValue(new DOMRect(left, top, 180, 300));
        vi.spyOn(color, 'getBoundingClientRect').mockReturnValue(new DOMRect(left, top + 30, 180, 24));
        color.click();
        await v.settle();
        const t = menusOf(container)[1].style.transform;
        headers.closeMenu();
        return t;
      };
      // innerWidth 1024, innerHeight 768 in jsdom: beside the row, its top a padding above it
      expect(await at(100, 50)).toBe('translate(282px, 74px)');
      // no room on the right: left of the menu
      expect(await at(900, 50)).toBe(`translate(${900 - 2 - 160}px, 74px)`);
      // low on the window: kept inside it
      const low = await at(100, 700);
      expect(low).toMatch(/^translate\(282px, \d+px\)$/);
      expect(Number(/, (\d+)px/.exec(low)![1]) + 100).toBeLessThanOrEqual(768 - 4);
    } finally {
      width.mockRestore();
      height.mockRestore();
    }
  });

  it('opens by ArrowRight; ← closes the flyout alone and focuses Color; Esc then closes the menu', async () => {
    const { v, container, headers } = await mount();
    const menu = await openMenu(v, headers, container, Letters);
    const color = row(menu, /^Color/);
    color.focus();
    key(color, 'ArrowRight');
    await v.settle();
    expect(menusOf(container)).toHaveLength(2);
    expect(levels(container)[0]).toContain('Hide track');
    expect(menusOf(container)[1].contains(document.activeElement)).toBe(true);
    key(document.activeElement!, 'ArrowLeft');
    expect(menusOf(container)).toHaveLength(1);
    expect(menu.isConnected).toBe(true);
    expect(document.activeElement).toBe(color);
    expect(color.getAttribute('aria-expanded')).toBe('false');
    expect(headers.menu.value).toBeDefined();
    // Esc in the flyout closes it alone, the focus on Color; Esc again closes the menu
    key(color, 'ArrowRight');
    await v.settle();
    key(document.activeElement!, 'Escape');
    expect(menusOf(container)).toHaveLength(1);
    expect(document.activeElement).toBe(color);
    key(color, 'Escape');
    await v.settle();
    expect(menusOf(container)).toHaveLength(0);
    expect(headers.menu.value).toBeUndefined();
  });

  it('opens by Enter and Space (the row’s own click) with the focus in the flyout', async () => {
    const { v, container, headers } = await mount();
    const color = row(await openMenu(v, headers, container, Letters), /^Color/);
    color.focus();
    // a keyboard click has detail 0: the focus goes into the flyout
    color.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true, detail: 0 }));
    await v.settle();
    expect(menusOf(container)).toHaveLength(2);
    expect(menusOf(container)[1].contains(document.activeElement)).toBe(true);
  });

  it('opens when the pointer rests on its row, the parent staying; another row’s rest closes it, arriving in it does not', async () => {
    const { v, container, headers, pointer } = await mount();
    const menu = await openMenu(v, headers, container, Letters);
    const color = row(menu, /^Color/);
    fakeTimers();
    pointer('pointerover', color);
    expect(menusOf(container)).toHaveLength(1);
    runTimers();
    await v.settle();
    expect(menusOf(container)).toHaveLength(2);
    expect(levels(container)[0]).toContain('Hide track');
    // onto another row of the menu: a rest, then the flyout closes
    pointer('pointerover', row(menu, /^Reverse/));
    expect(menusOf(container)).toHaveLength(2);
    // the pointer arrives in the flyout first (a diagonal move): nothing closes it
    menusOf(container)[1].dispatchEvent(new PointerEvent('pointerenter', { bubbles: false }));
    runTimers();
    expect(menusOf(container)).toHaveLength(2);
    pointer('pointerover', row(menu, /^Reverse/));
    runTimers();
    expect(menusOf(container)).toHaveLength(1);
    expect(color.getAttribute('aria-expanded')).toBe('false');
    expect(headers.menu.value).toBeDefined();
  });

  it('a menu dismissed while the pointer rests on it opens nothing after', async () => {
    const { v, container, headers, pointer } = await mount();
    for (const dismiss of ['Escape', 'outside', 'choice'] as const) {
      const menu = await openMenu(v, headers, container, Letters);
      const color = row(menu, /^Color/);
      fakeTimers();
      pointer('pointerover', color);
      if (dismiss === 'Escape') key(color, 'Escape');
      else if (dismiss === 'outside') pointer('pointerdown', document.body);
      else row(menu, /^Reverse/).click();
      await v.settle();
      expect(container.querySelector('[role="menu"]')).toBeNull();
      runTimers();
      await v.settle();
      expect(container.querySelector('[role="menu"]'), dismiss).toBeNull();
      expect(headers.menu.value).toBeUndefined();
      vi.mocked(DomHost.setTimeout).mockRestore();
      vi.mocked(DomHost.clearTimeout).mockRestore();
    }
  });

  it('a press outside closes the flyout and the menu; a press in the flyout closes nothing', async () => {
    const { v, container, headers, pointer } = await mount();
    const color = row(await openMenu(v, headers, container, Letters), /^Color/);
    color.click();
    await v.settle();
    pointer('pointerdown', menusOf(container)[1]);
    expect(menusOf(container)).toHaveLength(2);
    pointer('pointerdown', document.body);
    await v.settle();
    expect(menusOf(container)).toHaveLength(0);
    expect(headers.menu.value).toBeUndefined();
  });

  it('waits the `submenuRestMs` param; 0 turns resting off', async () => {
    const { v, container, headers, pointer } = await mount({ submenuRestMs: 120 });
    expect(headers.submenuRestMs).toBe(120);
    const color = row(await openMenu(v, headers, container, Letters), /^Color/);
    fakeTimers();
    pointer('pointerover', color);
    expect(vi.mocked(DomHost.setTimeout)).toHaveBeenLastCalledWith(expect.any(Function), 120);
    headers.closeMenu();
    vi.mocked(DomHost.setTimeout).mockClear();
    const off = await mount({ submenuRestMs: 0 });
    const c0 = row(await openMenu(off.v, off.headers, off.container, Letters), /^Color/);
    vi.mocked(DomHost.setTimeout).mockClear();
    off.pointer('pointerover', c0);
    expect(vi.mocked(DomHost.setTimeout)).not.toHaveBeenCalled();
  });

  it('a move within a submenu’s row keeps its rest: it opens at the first rest, not later', async () => {
    const { v, container, headers, pointer } = await mount();
    const color = row(await openMenu(v, headers, container, Letters), /^Color/);
    fakeTimers();
    pointer('pointerover', color);
    const n = [...timers.keys()];
    // onto the row's label (a child of the same row)
    pointer('pointerover', color.querySelector('span:not(:empty)') ?? color.lastElementChild!);
    expect([...timers.keys()]).toEqual(n);
  });
});

describe('Show hidden ▸ (the application’s way back to a hidden track), the same flyout', () => {
  it('is in the section’s menu only while a track is hidden; → lists each hidden track beside it, Enter shows one', async () => {
    const { v, container, headers } = await mount();
    const menu0 = await openMenu(v, headers, container, Letters);
    expect(menu0.textContent).not.toContain('Show hidden');
    headers.closeMenu();
    await v.execute('track.hide', { address: Domains });
    await v.settle();
    frames.flush();
    // a section with a hidden track: still no line, no chip; the menu has it
    expect(v.frame!.layout.sections[0].header).toBeUndefined();
    const menu = await openMenu(v, headers, container, Letters);
    const hidden = row(menu, /^Show hidden/);
    expect(hidden).toBeDefined();
    hidden.focus();
    key(hidden, 'ArrowRight');
    await v.settle();
    // the menu stays; one hidden track: no "Show all" beside it
    expect(menu.isConnected).toBe(true);
    expect(levels(container)[1]).toEqual(['Domains']);
    // Enter on a row is its click
    (document.activeElement as HTMLElement).click();
    await v.settle();
    expect(tracks(v)).toEqual([Letters, Domains]);
    expect(v.hiddenTracks.value).toEqual([]);
    expect(menusOf(container)).toHaveLength(0);
  });

  it('with several hidden: each in spec order, then Show all (one patch, one undo step)', async () => {
    const { v, container, headers } = await mount();
    // a third track stays shown: a section with every track hidden has no header to open a menu from
    const b = createSeqSpec({ title: 'three' });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
    b.views.view('main', (w) => {
      const s = w.section('protein', { of: p });
      for (const id of ['letters', 'domains', 'more']) s.track(id).layer('l', { type: 'letters' });
    });
    await v.load(b.build());
    await v.settle();
    frames.flush();
    const More = `${Section}/more`;
    await v.execute('track.hide', { address: Domains });
    await v.execute('track.hide', { address: Letters });
    await v.settle();
    frames.flush();
    expect(tracks(v)).toEqual([More]);
    expect(v.hiddenTracks.value.map((t) => t.address)).toEqual([Letters, Domains]);
    expect(headers.openMenu({ section: Section, track: More }, { x: 0, y: 0 })).toBe(true);
    await v.settle();
    frames.flush();
    const hidden = row(menusOf(container)[0], /^Show hidden/);
    hidden.focus();
    key(hidden, 'ArrowRight');
    await v.settle();
    expect(levels(container)[1]).toEqual(['letters', 'domains', 'Show all']);
    const steps = (v.spec.history as UiHistory).stack.steps.value.length;
    row(menusOf(container)[1], /^Show all/).click();
    await v.settle();
    expect(tracks(v)).toEqual([Letters, Domains, More]);
    expect((v.spec.history as UiHistory).stack.steps.value.length).toBe(steps + 1);
  });
});
