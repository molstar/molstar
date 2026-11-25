// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { BuiltInSchemes } from '@molstar/seqx-core/theme/index';
import { LetterSchemes } from '@molstar/seqx-examples/p53/themes';
import { DomHost } from '@molstar/ui/dom';
import { loadLegendEyes } from '@molstar/seqx-viewer/mount-legend';
import { afterAll, afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import type { EditingSession } from '../../editing/compose.js';
import { EditingThemes, SchemesGallery, Themes } from '../themes.stories.js';
import { SchemesGalleryRoot, type ThemesEditingRoot } from '../themes-roots.js';
import {
  cancelScheme,
  Letters,
  pickerSchemes,
  pickScheme,
  previewScheme,
  setSignalDomain,
  signalDomains,
} from '../themes-compose.js';

// The p53 themes stories in a DOM (jsdom) on a fake clock, each mounted as Storybook mounts it (the story's `render`,
// its element connected to the page), the edits driven through the composition the story made: a scheme previewed
// while hovered (no step) and picked in one patch; the shared scale's domain set for both signals in one step and undone
// in one; the export's legends below the data.

const view = FeatureViews.find((v) => v.id === 'themes-editing')!;
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

/** A story's element as Storybook mounts it: rendered, then connected to the page. */
function mount<T extends HTMLElement>(story: { render?: unknown }): T {
  const el = (story.render as () => T)();
  document.body.append(el);
  mounted = el;
  return el;
}

async function make(): Promise<EditingSession> {
  const el = mount<ThemesEditingRoot>(EditingThemes);
  // the story renders once connected: a root with no model would render nothing (and so would this)
  expect(el.session).toBeDefined();
  await el.loaded;
  await vi.advanceTimersByTimeAsync(400);
  return el.session!;
}

/** The classes track's letters' color as shown (the preview included), printed. */
const lettersColor = (s: EditingSession) => {
  const layer = resolveAddress(s.viewer.spec.shown.value, Letters.split('/').filter(Boolean));
  const c = layer?.params?.color;
  return c === undefined ? '' : printSexpr(c as Node);
};

// The legend's eyes are a lazy chunk the first mount in a worker starts: the file waits for it, or it lands after the file ends.
afterAll(() => loadLegendEyes());

describe('p53 Themes, headless', () => {
  it('both views are listed, the same spec', () => {
    expect(P53Views.find((v) => v.id === 'themes')?.build).toBe(view.build);
  });

  it('a scheme previewed is no step, picked is one; undo restores the default', { tags: ['medium'] }, async () => {
    const s = await make();
    const offered = pickerSchemes(s);
    for (const name of LetterSchemes) expect(offered).toContain(name);
    expect(offered).not.toContain('nucleotide');

    const before = lettersColor(s);
    expect(previewScheme(s, 'zappo')).toBe(true);
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toContain('zappo');
    expect(s.stack.steps.value).toHaveLength(0);
    cancelScheme(s);
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toBe(before);

    previewScheme(s, 'taylor');
    await pickScheme(s, 'taylor');
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toContain('taylor');
    expect(s.stack.steps.value).toHaveLength(1);
    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toBe(before);
  });

  it('the shared scale’s domain: both signals follow in one step, undone in one', { tags: ['medium'] }, async () => {
    const s = await make();
    const [a, b] = signalDomains(s);
    // open at the top: the union of both signals (the accessibility's range is the larger)
    expect(a).toBeDefined();
    expect(a).toEqual(b);
    expect(a![0]).toBe(0);
    expect(a![1]).toBeGreaterThan(1);

    await setSignalDomain(s, [0, 0.5]);
    await vi.advanceTimersByTimeAsync(400);
    expect(signalDomains(s)).toEqual([
      [0, 0.5],
      [0, 0.5],
    ]);
    expect(s.stack.steps.value).toHaveLength(1);

    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(400);
    expect(signalDomains(s)).toEqual([a, b]);
  });

  it('the legends: one per vocabulary, the shared scale’s merged, and the export draws them below', {
    tags: ['medium'],
  }, async () => {
    const s = await make();
    const legends = s.viewer.legends.legends.value;
    const signal = legends.filter((l) => l.scale === '/views/signal');
    expect(signal).toHaveLength(1);
    expect(signal[0].layers).toHaveLength(2);
    expect(legends.some((l) => l.field === 'consequence')).toBe(true);
    const p = s.viewer.renderSvg({ idPrefix: 'th' });
    await vi.advanceTimersByTimeAsync(300);
    const svg = (await p).svg;
    expect(svg).toContain('data-seq-scale-');
    expect(svg).toMatch(/missense/);
  });

  it('the view mounts with content', () => {
    const t = mount<HTMLElement>(Themes);
    expect(t.querySelectorAll('*').length).toBeGreaterThan(3);
    t.remove();
  });

  // the gallery's quick representative: two schemes, each table a row of four visions' swatches
  it('the gallery mounts two schemes’ rows', () => {
    const el = document.createElement(SchemesGalleryRoot.tag) as SchemesGalleryRoot;
    el.model = BuiltInSchemes.slice(0, 2);
    document.body.append(el);
    mounted = el;
    const rows = el.querySelectorAll('tr');
    // a header, then a row a table a scheme has colors in
    expect(rows.length).toBeGreaterThan(2);
    expect(rows[1].querySelectorAll('td').length).toBe(2 + 4);
  });

  // every scheme in both tables under four visions, a swatch a color, in jsdom: a third of a second (a sweep)
  it('the gallery mounts with content', { tags: ['medium'] }, () => {
    const g = mount<SchemesGalleryRoot>(SchemesGallery);
    expect(g.querySelectorAll('tr').length).toBeGreaterThan(20);
  });

  // the editing root loads its session: slow by nature (the quick tests above each mount it too, and read its state)
  it('the editing root mounts with its picker', { tags: ['medium'] }, async () => {
    const e = mount<ThemesEditingRoot>(EditingThemes);
    await e.loaded;
    expect(e.querySelectorAll('[role="option"]').length).toBeGreaterThan(5);
  });
});
