/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigator's wash (`seq-wash`, _spec/README.md "Render", Navigators): light enough that the picture under it stays
// legible in both themes (an overview whose window is a few floored rows of 10,000 is washed nearly whole: its shape
// must still read); the window stands out by its edge and grips. A bound on the wash's alpha, and every picture color
// (the residue classes, the ramp, the categories, a chromosome's palest stain, the centromere) keeping most of its
// contrast against the empty background once washed (or 3:1 against it, WCAG's for graphics).

import { describe, expect, it } from 'vitest';
import { DarkTokens, LightTokens, type SeqTokenName, type SeqTokens } from '../index.js';

/** The wash's alpha at most: a third of the picture's contrast lost, no more. */
const MaxWashAlpha = 0.3;
/** The share of a color's contrast (its WCAG ratio less 1) it keeps once washed, at least. */
const MinKept = 0.6;

const Picture: readonly SeqTokenName[] = [
  'seq-aa-hydrophobic',
  'seq-aa-aromatic',
  'seq-aa-polar',
  'seq-aa-positive',
  'seq-aa-negative',
  'seq-aa-special',
  'seq-ramp-0',
  'seq-ramp-2',
  'seq-ramp-4',
  'seq-cat-0',
  'seq-cat-3',
  'seq-cat-7',
  'seq-centromere',
];
/** A chromosome's palest stain (`gpos25` in the `giemsa` scheme), per theme. */
const PalestStain = { light: 0xc8c8c8ff, dark: 0x4a4f58ff };

const rgb = (c: number) => [(c >>> 24) & 255, (c >>> 16) & 255, (c >>> 8) & 255] as const;
const alpha = (c: number) => (c & 255) / 255;

function luminance([r, g, b]: readonly number[]): number {
  const f = (v: number) => {
    const s = v / 255;
    return s <= 0.03928 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  };
  return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
}

function ratio(a: readonly number[], b: readonly number[]): number {
  const x = luminance(a);
  const y = luminance(b);
  return (Math.max(x, y) + 0.05) / (Math.min(x, y) + 0.05);
}

/** A color under the wash, over it. */
function washed(c: number, wash: number): number[] {
  const a = alpha(wash);
  const w = rgb(wash);
  return rgb(c).map((v, i) => v * (1 - a) + w[i] * a);
}

describe('the navigator’s wash keeps the picture legible', () => {
  for (const [name, tokens] of [
    ['light', LightTokens],
    ['dark', DarkTokens],
  ] as [keyof typeof PalestStain, SeqTokens][])
    it(`${name}: its alpha bounded, every picture color keeping most of its contrast against the background`, () => {
      const wash = tokens.colors['seq-wash'];
      expect(alpha(wash)).toBeLessThanOrEqual(MaxWashAlpha);
      const bg = rgb(tokens.colors['msui-bg']);
      for (const c of [...Picture.map((t) => tokens.colors[t]), PalestStain[name]]) {
        const raw = ratio(rgb(c), bg) - 1;
        const kept = ratio(washed(c, wash), bg) - 1;
        // a strong color may lose more of a high ratio and still read (3:1, WCAG's for graphics)
        if (kept + 1 < 3) expect(kept / raw, c.toString(16)).toBeGreaterThanOrEqual(MinKept);
      }
    });
});
