// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { rondo, rondoRing } from '@molstar/seqx-ext-rondo';
import { partnerChromosomes, rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { type HeadersBehavior, HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { headerMenuOverlay } from '@molstar/seqx-viewer/mount-header-menu';
import { type UiHistory, uiHistory } from '@molstar/seqx-viewer/ui-history';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ribbonsOf, strongestStem } from '../views.js';

// The stem menu in a DOM (jsdom), driven by the event sequences a browser sends for a click (pointermove or
// pointerover, pointerdown, pointerup, click): "Show with a partner" has a submenu, and a click on its row opens it as
// the Right arrow does.

/** A canvas's 2D context that takes every call. */
function context(canvas: HTMLCanvasElement): unknown {
  const target: Record<string | symbol, unknown> = {
    canvas,
    measureText: (t: string) => ({ width: 6 * t.length, actualBoundingBoxAscent: 8, actualBoundingBoxDescent: 2 }),
    getImageData: (_x: number, _y: number, w: number, h: number) => ({
      data: new Uint8ClampedArray(4 * w * h),
      width: w,
      height: h,
    }),
    createImageData: (w: number, h: number) => ({ data: new Uint8ClampedArray(4 * w * h), width: w, height: h }),
    getTransform: () => ({ a: 1, b: 0, c: 0, d: 1, e: 0, f: 0 }),
    createLinearGradient: () => ({ addColorStop() {} }),
    createRadialGradient: () => ({ addColorStop() {} }),
    createPattern: () => ({}),
  };
  return new Proxy(target, {
    get: (t, k) => (k in t ? t[k] : () => {}),
    set: (t, k, v) => {
      t[k] = v;
      return true;
    },
  });
}

const raf = DomHost.requestAnimationFrame;
let viewer: SeqViewer | undefined;
let target: HTMLElement | undefined;

beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as unknown as { __ctx?: unknown }).__ctx ??= context(this)) as never;
  });
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  target?.remove();
  target = undefined;
  DomHost.requestAnimationFrame = raf;
  vi.useRealTimers();
  vi.restoreAllMocks();
});

/** Six chromosomes and 300 links between them (seeded by arithmetic: the pairs around the ring). */
function spec() {
  const chroms = Array.from({ length: 6 }, (_, i) => ({ name: `chr${i + 1}`, length: 8_000_000 - i * 700_000 }));
  const rows = {
    a_name: [] as string[],
    b_name: [] as string[],
    a_start: [] as number[],
    a_end: [] as number[],
    b_start: [] as number[],
    b_end: [] as number[],
  };
  for (let k = 0; k < 300; k++) {
    const a = k % 6;
    const b = (a + 1 + (k % 5)) % 6;
    rows.a_name.push(chroms[a].name);
    rows.b_name.push(chroms[b].name);
    const s = (k * 7919) % (chroms[a].length - 1000);
    const t = (k * 104729) % (chroms[b].length - 1000);
    rows.a_start.push(s);
    rows.a_end.push(s + 500);
    rows.b_start.push(t);
    rows.b_end.push(t + 500);
  }
  const b = createSeqSpec({ title: 'stem menu', extensions: [rondo] });
  const g = b.data.assembly('toy', assembly('Toy', chroms));
  const links = b.data.pairs('links', { a: g, rows } as never);
  b.views.view('main', (v) => {
    rondoRing(v, { assembly: g, links, square: false });
  });
  return b.build();
}

const menusOf = (root: ParentNode) => [...root.querySelectorAll<HTMLElement>('[role="menu"]')];

/** The viewer mounted over the toy ring (the host's history, the menu's port), settled, its parts landed. */
async function mounted() {
  target = document.createElement('div');
  document.body.append(target);
  const base = DefaultSeqViewerSpec();
  const v = (viewer = createSeqViewer(
    {
      ...base,
      ui: { canvas: true },
      extensions: [...(base.extensions ?? []), rondoViewer],
      config: [[SeqCoreConfig.render, { labelWidth: 140 }]],
      behaviors: base.behaviors.map((b) =>
        b.provider === HeadersProvider ? { provider: HeadersProvider, params: { menu: headerMenuOverlay() } } : b,
      ),
      history: uiHistory(),
    },
    target,
  ));
  void v.load(spec());
  v.resize(600);
  await Promise.resolve();
  // a lazy part the spec's model waits for: its import is real I/O, never the fake clock's
  for (
    let seen = new Set<unknown>(), m = v.spec.model.value;
    m.loading && !seen.has(m.loading);
    m = v.spec.model.value
  ) {
    seen.add(m.loading);
    await m.loading.catch(() => undefined);
  }
  for (let i = 0; i < 20; i++) {
    await vi.advanceTimersByTimeAsync(100);
    v.draw();
  }
  const headers = v.behavior<HeadersBehavior>('headers')!;
  await headers.loaded();
  const l = ribbonsOf(v)!;
  expect(l.ns).toBeGreaterThan(0);
  const container = target.querySelector<HTMLElement>('[role="application"]')!;
  const canvas = container.querySelector<HTMLCanvasElement>('canvas')!;
  return { v, headers, l, container, canvas };
}

/** The menu's DOM landed (a lazy chunk, loaded at its first open). */
async function menuLanded(headers: HeadersBehavior) {
  await vi.advanceTimersByTimeAsync(300);
  await headers.loaded();
  await vi.advanceTimersByTimeAsync(300);
}

/** A secondary click on the canvas at a point (jsdom's container at the page's origin): whether it was prevented. */
function secondaryClick(canvas: HTMLCanvasElement, x: number, y: number, shiftKey = false): boolean {
  const e = new MouseEvent('contextmenu', {
    clientX: x,
    clientY: y,
    button: 2,
    bubbles: true,
    cancelable: true,
    shiftKey,
  });
  canvas.dispatchEvent(e);
  return e.defaultPrevented;
}

/** A menu row by its text. */
const rowOf = (root: ParentNode, text: RegExp) =>
  [...root.querySelectorAll<HTMLElement>('[role^="menuitem"]')].find((b) => text.test(b.textContent ?? ''))!;

/** A key pressed where the focus is (the menu's rows take it). */
const press = (key: string) =>
  (document.activeElement ?? document.body).dispatchEvent(
    new KeyboardEvent('keydown', { key, bubbles: true, cancelable: true }),
  );

/** The viewer mounted over the toy ring, its stem menu open on the strongest stem: the container, the row, its events. */
async function openStemMenu() {
  const { headers, l, container } = await mounted();
  const { x, y } = strongestStem(l);
  expect(headers.context(x, y)).toBe(true);
  await menuLanded(headers);
  expect(menusOf(container)).toHaveLength(1);
  const row = rowOf(container, /Show with a partner/);
  expect(row.getAttribute('aria-haspopup')).toBe('menu');
  const r = row.getBoundingClientRect();
  const xy = { clientX: r.left + 2, clientY: r.top + 2 };
  const pointer = (type: string) =>
    row.dispatchEvent(
      new PointerEvent(type, { ...xy, bubbles: true, cancelable: true, button: 0, pointerId: 1, pointerType: 'mouse' }),
    );
  const click = () =>
    row.dispatchEvent(new MouseEvent('click', { ...xy, bubbles: true, cancelable: true, button: 0, detail: 1 }));
  return { container, row, pointer, click };
}

describe('The stem menu, by a pointer', () => {
  it('a click on "Show with a partner" opens its submenu, the menu staying, as the Right arrow does', async () => {
    const { container, row, pointer, click } = await openStemMenu();
    // the events of a click: over, move, down, up, click
    pointer('pointerover');
    pointer('pointermove');
    pointer('pointerdown');
    pointer('pointerup');
    click();
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(2);
    expect(row.getAttribute('aria-expanded')).toBe('true');
    expect(menusOf(container)[1].querySelectorAll('[role^="menuitem"]').length).toBeGreaterThan(0);
  });

  it('a click after the pointer rested on the row (the flyout already open) keeps it open', async () => {
    const { container, row, pointer, click } = await openStemMenu();
    pointer('pointerover');
    await vi.advanceTimersByTimeAsync(500);
    expect(menusOf(container)).toHaveLength(2);
    pointer('pointermove');
    pointer('pointerdown');
    pointer('pointerup');
    click();
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(2);
    expect(row.getAttribute('aria-expanded')).toBe('true');
  });
});

/** The ring section's spec node (its `include`, what a partner's row sets). */
const ringNode = (v: SeqViewer) =>
  v.spec.root.value.children!.find((c) => c.kind === 'views')!.children![0].children!.find((c) => c.id === 'ring')!;

// The stem menu's gestures (the ring's design, the stem menu's table), each by its name, in a DOM: the browser's events
// on the mounted canvas, the keys on the view and on the menu's rows.
describe('The stem menu’s gestures', () => {
  it('right-click a stem: the stem menu at the pointer', async () => {
    const { headers, l, canvas } = await mounted();
    const { x, y } = strongestStem(l);
    expect(secondaryClick(canvas, x, y)).toBe(true);
    const shown = headers.menu.value!;
    expect(shown.place).toEqual({ x, y });
    expect(shown.menu.target.element?.type).toBe('rondo:ribbons');
    expect(shown.menu.groups[0].map((i) => i.name)).toContain('rondo:show-partner');
    // Shift+secondary click stays the browser's
    headers.closeMenu();
    expect(secondaryClick(canvas, x, y, true)).toBe(false);
    expect(headers.menu.value).toBeUndefined();
  });

  it('right-click the ring’s empty space: the browser’s menu', async () => {
    const { v, headers, l, canvas } = await mounted();
    // a point of the ring's section that hits nothing (outside the ring, or the hole between ribbons)
    const box = v.frame!.layout.sections.find((s) => s.section.address === l.section)!.placement.box;
    let empty: [number, number] | undefined;
    for (let y = box.y + 30; y < box.y + box.height && !empty; y += 7)
      for (let x = box.x + 2; x < box.x + box.width && !empty; x += 7)
        if (!v.hit(x, y) && v.input.sectionAt(y) === l.section) empty = [x, y];
    expect(empty).toBeDefined();
    expect(secondaryClick(canvas, empty![0], empty![1])).toBe(false);
    expect(headers.menu.value).toBeUndefined();
  });

  it('focus a stem, Shift+F10 or the Menu key: the stem menu at the stem', async () => {
    const { v, headers, l } = await mounted();
    const { stem, x, y } = strongestStem(l);
    for (const key of [{ key: 'F10', shiftKey: true }, { key: 'ContextMenu' }]) {
      // the keyboard cursor on the stem (a click on it, or the cursor's keys, puts it there)
      v.selection.setCursor({
        section: l.section,
        system: l.systems[l.stemSystem[stem]],
        position: l.stemStart[stem],
        keyboard: true,
      });
      expect(v.input.key(key)).toBe(true);
      const shown = headers.menu.value!;
      expect(shown.menu.target.element).toEqual({ layer: l.key, key: l.stemKey[stem], type: 'rondo:ribbons' });
      expect(shown.place.x).toBeCloseTo(x, 6);
      expect(shown.place.y).toBeCloseTo(y, 6);
      headers.closeMenu();
      // past the time the browser's own secondary click after the Menu key is the same ask (EC-156)
      await vi.advanceTimersByTimeAsync(1100);
    }
    // on a track's header the menu key opens the header's menu, as before
    v.selection.setCursor({
      section: l.section,
      track: l.key.slice(0, l.key.lastIndexOf('/')),
      system: l.systems[l.stemSystem[stem]],
      position: l.stemStart[stem],
    });
    expect(v.input.key({ key: 'F10', shiftKey: true })).toBe(true);
    expect(headers.menu.value!.menu.target.element).toBeUndefined();
  });

  it('click a partner: shows the pair, one history step', async () => {
    const { v, headers, l, container, canvas } = await mounted();
    const { stem, x, y } = strongestStem(l);
    secondaryClick(canvas, x, y);
    await menuLanded(headers);
    expect(ringNode(v).params?.include).toBeUndefined();
    rowOf(container, /Show with a partner/).focus();
    press('ArrowRight');
    await vi.advanceTimersByTimeAsync(50);
    const partner = menusOf(container)[1].querySelector<HTMLElement>('[role^="menuitem"]')!;
    const name = l.names[l.stemSystem[stem]];
    const other = /(chr\d)/.exec(partner.textContent ?? '')![1];
    // a chromosome with pieces opens them: a piece shows the same pair
    if (partner.getAttribute('aria-haspopup') === 'menu') {
      partner.focus();
      press('ArrowRight');
      await vi.advanceTimersByTimeAsync(50);
      menusOf(container)[2].querySelector<HTMLElement>('[role^="menuitem"]')!.click();
    } else partner.click();
    await vi.advanceTimersByTimeAsync(300);
    await v.settle();
    expect(new Set(ringNode(v).params?.include as string[])).toEqual(new Set([name, other]));
    expect(menusOf(container)).toHaveLength(0);
    await (v.spec.history as UiHistory).stack.undo();
    await v.settle();
    expect(ringNode(v).params?.include).toBeUndefined();
  });

  it('hover or Right arrow on a chromosome with pieces: its submenu', async () => {
    const { headers, l, container, canvas } = await mounted();
    // a stem with a partner chromosome reached by more than one ribbon
    let stem = -1;
    for (let s = 0; s < l.ns && stem < 0; s++)
      if (l.stemPlaced[s] && l.stemPoly[s] >= 0 && partnerChromosomes(l, s).some((c) => c.pieces.length > 1)) stem = s;
    expect(stem).toBeGreaterThanOrEqual(0);
    const c = partnerChromosomes(l, stem).find((x) => x.pieces.length > 1)!;
    l.frame.point((l.stemT0[stem] + l.stemT1[stem]) / 2, 2);
    expect(secondaryClick(canvas, l.frame.px, l.frame.py)).toBe(true);
    await menuLanded(headers);
    rowOf(container, /Show with a partner/).focus();
    press('ArrowRight');
    await vi.advanceTimersByTimeAsync(50);
    const name = c.within ? `Within ${l.names[c.system]}` : l.names[c.system];
    const row = rowOf(menusOf(container)[1], new RegExp(`^${name} · `));
    expect(row.getAttribute('aria-haspopup')).toBe('menu');
    // the Right arrow on it: its pieces, each named by its range (the toy has no cytobands: in bases)
    row.focus();
    press('ArrowRight');
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(3);
    const pieces = [...menusOf(container)[2].querySelectorAll('[role^="menuitem"]')].map((p) => p.textContent);
    expect(pieces).toHaveLength(c.pieces.length);
    for (const p of pieces) expect(p).toMatch(/^[\d.]+–[\d.]+ (Mb|kb) · [\d,]+ links?/);
    // ← back to the row, its submenu closed; the pointer resting on the row opens it again
    press('ArrowLeft');
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(2);
    const r = row.getBoundingClientRect();
    row.dispatchEvent(
      new PointerEvent('pointerover', {
        clientX: r.left + 2,
        clientY: r.top + 2,
        bubbles: true,
        pointerId: 1,
        pointerType: 'mouse',
      }),
    );
    await vi.advanceTimersByTimeAsync(500);
    expect(menusOf(container)).toHaveLength(3);
  });

  it('Escape: closes the submenu, then the menu', async () => {
    const { container, row } = await openStemMenu();
    row.focus();
    press('ArrowRight');
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(2);
    press('Escape');
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(1);
    expect(document.activeElement).toBe(row);
    press('Escape');
    await vi.advanceTimersByTimeAsync(50);
    expect(menusOf(container)).toHaveLength(0);
  });
});
