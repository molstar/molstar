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
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { type UiHistory, uiHistory } from '@molstar/seqx-viewer/ui-history';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

// The square's hint beside the ring as the default UI places it (Features/Genome Ring, The Square Beside): one muted
// line under the square, its words wrapped to the square's width.

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
    rondoRing(v, { assembly: g, links });
  });
  return b.build();
}

async function mounted() {
  target = document.createElement('div');
  document.body.append(target);
  const base = DefaultSeqViewerSpec();
  const v = (viewer = createSeqViewer(
    { ...base, extensions: [...(base.extensions ?? []), rondoViewer], history: uiHistory() },
    target,
  ));
  void v.load(spec());
  v.resize(600);
  await Promise.resolve();
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
  return v;
}

describe('the square beside the ring', () => {
  it('the default placement draws the square with its hint under it', async () => {
    await mounted();
    const hints = [...target!.querySelectorAll<HTMLElement>('.seqx-navigator-hint')];
    expect(target!.querySelectorAll('[data-navigator]').length).toBe(1);
    expect(hints.length).toBe(1);
    expect(hints[0].textContent).toContain('double-click for the whole genome');
    // muted by the token (light and dark follow it), wrapped to the square's width
    expect(hints[0].style.color).toBe('var(--msui-text-muted, gray)');
    expect(hints[0].style.width).toBe(`${220}px`);
    // under the square, in its cell
    expect(hints[0].previousElementSibling?.hasAttribute('data-navigator')).toBe(true);
  });
});
