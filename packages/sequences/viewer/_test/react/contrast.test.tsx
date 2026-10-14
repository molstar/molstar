/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The picked item's ink on its wash (seqspec:interaction#the-skin): the tab shown, the choice in effect and the snapshot
// shown read at WCAG AA (4.5:1, 11-12px text) in light and dark, measured as drawn: the wash over the strip's glass over
// the canvas's surface, composited in a canvas the way the page composites them.

import { act } from 'react';
import { afterEach, describe, expect, it } from 'vitest';
import '../../src/skin.css';
import { harness, settle } from './fixtures.js';

const h = harness();
afterEach(() => {
  delete document.documentElement.dataset.theme;
});

const g = new OffscreenCanvas(1, 1).getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
/** The colors over each other, bottom first, on white; the result as sRGB 0..255. */
function composite(colors: readonly string[]): [number, number, number] {
  g.globalCompositeOperation = 'source-over';
  g.fillStyle = '#fff';
  g.fillRect(0, 0, 1, 1);
  for (const c of colors) {
    g.fillStyle = c;
    g.fillRect(0, 0, 1, 1);
  }
  const [r, gr, b] = g.getImageData(0, 0, 1, 1).data;
  return [r, gr, b];
}
const lum = ([r, gr, b]: readonly number[]) => {
  const c = (v: number) => (v / 255 <= 0.04045 ? v / 255 / 12.92 : ((v / 255 + 0.055) / 1.055) ** 2.4);
  return 0.2126 * c(r) + 0.7152 * c(gr) + 0.0722 * c(b);
};
/** The text's contrast on its ground: every background from the viewer's root up to the element, composited. */
function contrast(e: HTMLElement): number {
  const grounds: string[] = [];
  for (let n: HTMLElement | null = e; n; n = n.classList.contains('seqv-root') ? null : n.parentElement)
    grounds.unshift(getComputedStyle(n).backgroundColor);
  const ground = composite(grounds);
  const ink = composite([...grounds, getComputedStyle(e).color]);
  const [a, b] = [lum(ink), lum(ground)].sort((x, y) => y - x);
  return (a + 0.05) / (b + 0.05);
}

describe("the picked item's contrast", () => {
  for (const theme of ['light', 'dark'])
    it(`${theme}: the tab shown, the choice in effect and the snapshot shown read at AA on their wash`, async () => {
      document.documentElement.dataset.theme = theme;
      const t = await h.opened({ chrome: 'full' });
      await act(async () => void (await t.ctx.snapshots.apply('s1')));
      await settle();
      const picked = [
        t.container.querySelector('.seqv-tab[aria-selected="true"]'),
        t.container.querySelector('.seqv-choice[aria-pressed="true"]'),
        t.container.querySelector('.seqv-dot[aria-current]'),
      ] as HTMLElement[];
      expect(picked.every((e) => e !== null)).toBe(true);
      for (const e of picked) expect([e.className, contrast(e)]).toEqual([e.className, expect.any(Number)]);
      for (const e of picked) expect(contrast(e), e.className).toBeGreaterThanOrEqual(4.5);
    });
});
