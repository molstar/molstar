/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Black, contrastRgba, inkFor, White } from '../colors.js';

// the memoized ink equals the direct contrast rule, over fills that collide in the memo and across ink pairs
function direct(fill: number, dark: number, light: number): number {
  const cd = contrastRgba(fill, dark);
  const cl = contrastRgba(fill, light);
  if (Math.max(cd, cl) >= 4.5) return cd >= cl ? dark : light;
  return contrastRgba(fill, Black) >= contrastRgba(fill, White) ? Black : White;
}

describe('inkFor', () => {
  it('equals the contrast rule, first and memoized, for any pair of inks', () => {
    const pairs: [number, number][] = [
      [0x1a1a1aff, 0xf5f5f5ff],
      [0x555555ff, 0xaaaaaaff],
    ];
    for (let pass = 0; pass < 2; pass++)
      for (const [dark, light] of pairs)
        for (let i = 0; i < 5000; i++) {
          const fill = ((Math.imul(i, 0x9e3779b1) & 0xffffff00) | 0xff) >>> 0;
          expect(inkFor(fill, dark, light)).toBe(direct(fill, dark, light));
        }
  });

  it('ignores the fill alpha', () => {
    expect(inkFor(0x336699ff, Black, White)).toBe(inkFor(0x33669980, Black, White));
  });
});
