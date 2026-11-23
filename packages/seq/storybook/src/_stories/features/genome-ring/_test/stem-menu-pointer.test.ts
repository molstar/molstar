// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { rondo, rondoRing } from '@molstar/seqx-ext-rondo';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { type HeadersBehavior, HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { headerMenuOverlay } from '@molstar/seqx-viewer/mount-header-menu';
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

/** The viewer mounted over the toy ring, its stem menu open on the strongest stem: the container, the row, its events. */
async function openStemMenu() {
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
  const { x, y } = strongestStem(l);
  expect(headers.context(x, y)).toBe(true);
  // the menu's DOM is a lazy chunk, loaded at its first open
  await vi.advanceTimersByTimeAsync(300);
  await headers.loaded();
  await vi.advanceTimersByTimeAsync(300);
  const container = target.querySelector<HTMLElement>('[role="application"]')!;
  expect(menusOf(container)).toHaveLength(1);
  const row = [...container.querySelectorAll<HTMLElement>('[role^="menuitem"]')].find((b) =>
    /Show with a partner/.test(b.textContent ?? ''),
  )!;
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
