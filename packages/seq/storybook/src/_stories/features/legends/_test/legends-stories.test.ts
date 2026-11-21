// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { placementsOf, trackOf } from '../legends-compose.js';
import { type LegendPlacementsRoot, LegendPlacements } from '../legends.stories.js';

// The legend placements story in a DOM (jsdom) on a fake clock, mounted as Storybook mounts it: the same legends in
// three places by the host's options (the header's one track compact, the panel without it, the band compact in a row)
// and in the export's band below its figure.

let mounted: HTMLElement | undefined;

/** jsdom has no 2D context: one that takes every call and measures text by its length (nothing drawn is checked). */
function stubContext(canvas: HTMLCanvasElement): unknown {
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
beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as unknown as { __ctx?: unknown }).__ctx ??= stubContext(this)) as never;
  });
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  mounted?.remove();
  mounted = undefined;
  vi.useRealTimers();
});

const titles = (el: HTMLElement) => [...el.querySelectorAll('[role="group"]')].map((g) => g.getAttribute('aria-label'));

describe('Legends placements, headless', () => {
  it('a track’s address from a layer’s', () => {
    expect(trackOf('/views/main/p/domains/d')).toBe('/views/main/p/domains');
  });

  it('the same legends three ways and in the export', { tags: ['medium'] }, async () => {
    const el = (LegendPlacements.render as () => LegendPlacementsRoot)();
    document.body.append(el);
    mounted = el;
    await el.loaded;
    await vi.advanceTimersByTimeAsync(600);
    const s = el.session!;
    const hosts = el.hosts!;
    const all = s.viewer.legends.legends.value.map((l) => l.title);
    expect(all.length).toBeGreaterThan(1);
    const p = placementsOf(s.viewer)!;
    // the header's one track, compact (labels are names, not text)
    expect(titles(hosts.header)).toHaveLength(1);
    expect(hosts.header.querySelector<HTMLElement>('[data-form]')!.dataset.form).toBe('compact');
    // the panel has the others, full
    expect(titles(hosts.panel)).toHaveLength(all.length - 1);
    expect(titles(hosts.panel)).not.toContain(titles(hosts.header)[0]);
    expect(hosts.panel.querySelector<HTMLElement>('[data-form]')!.dataset.form).toBe('full');
    // the band has every legend in a row
    expect(titles(hosts.band)).toEqual(all);
    expect(hosts.band.querySelector<HTMLElement>('[data-orient]')!.dataset.orient).toBe('horizontal');
    expect(p.band.form).toBe('compact');
    // the export's band, below the figure
    const svg = (await s.viewer.renderSvg({ legend: 'below' })).svg;
    expect(svg).toMatch(/missense/);
  });
});
