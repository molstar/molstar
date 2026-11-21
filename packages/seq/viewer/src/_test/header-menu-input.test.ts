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
// move no cursor, and `Color ▸` opens its submenu by a click.

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type HeadersBehavior, HeadersProvider } from '../behaviors/headers.js';
import type { HoverCardBehavior } from '../behaviors/hover-card.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { headerMenuOverlay } from '../mount-header-menu.js';
import { uiHistory } from '../ui-history.js';
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
let timers: (() => void)[];
const viewers: SeqViewer[] = [];
beforeEach(() => {
  frames = installFrames();
  timers = [];
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
  vi.spyOn(DomHost, 'setTimeout').mockImplementation((fn) => (timers.push(fn), timers.length));
  vi.spyOn(DomHost, 'clearTimeout').mockImplementation((id) => {
    if (typeof id === 'number' && id > 0) timers[id - 1] = () => {};
  });
};
const runTimers = () => {
  for (const fn of timers.splice(0)) fn();
};

async function mount() {
  const target = document.createElement('div');
  document.body.appendChild(target);
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer(
    {
      ...base,
      ui: { canvas: true, history: uiHistory },
      config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
      behaviors: base.behaviors.map((b) =>
        b.provider === HeadersProvider ? { provider: HeadersProvider, params: { menu: headerMenuOverlay() } } : b,
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

describe('Color ▸', () => {
  const labels = (container: HTMLElement) =>
    [...container.querySelectorAll('[role="menu"] [role^="menuitem"]')].map((b) => b.textContent);

  it('opens its submenu by a click, and a choice applies the scheme', async () => {
    const { v, container, headers, press } = await mount();
    const g = geo(v, Domains);
    press(over(row(await openMenu(v, headers, container, Letters), /^Color/), v.dataX + 12, g.line.y + 2));
    await v.settle();
    frames.flush();
    expect(labels(container)).not.toContain('Hide track');
    const choices = [...container.querySelectorAll<HTMLElement>('[role="menu"] [role="menuitemradio"]')];
    expect(choices.length).toBeGreaterThan(1);
    const other = choices.find((c) => c.getAttribute('aria-checked') !== 'true')!;
    const name = other.textContent;
    press(over(other, v.dataX + 12, g.line.y + 2));
    await v.settle();
    expect(container.querySelector('[role="menu"]')).toBeNull();
    press(row(await openMenu(v, headers, container, Letters), /^Color/));
    await v.settle();
    expect(container.querySelector('[role="menu"] [aria-checked="true"]')?.textContent).toBe(name);
  });
});
