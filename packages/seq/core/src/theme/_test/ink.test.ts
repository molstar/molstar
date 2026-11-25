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

// the sweeps count mismatches and assert once: an `expect` a pair costs more than the pair
function sweepPairs(fills: number): number {
  const pairs: [number, number][] = [
    [0x1a1a1aff, 0xf5f5f5ff],
    [0x555555ff, 0xaaaaaaff],
    // one ink shared, the other different: a memo that ignores either ink returns the other pair's answer
    [0x1a1a1aff, 0x777777ff],
    [0x888888ff, 0xf5f5f5ff],
  ];
  let bad = 0;
  // a fill at a time, every pair: pairs whose slots collide evict each other, so the memo's ink compares are read
  for (let pass = 0; pass < 2; pass++)
    for (let i = 0; i < fills; i++) {
      const fill = ((Math.imul(i, 0x9e3779b1) & 0xffffff00) | 0xff) >>> 0;
      for (const [dark, light] of pairs) if (inkFor(fill, dark, light) !== direct(fill, dark, light)) bad++;
    }
  return bad;
}

function sweepSlots(fills: number): number {
  // 3,000 pairs over 1,024 slots collide: a slot is matched on both inks, so a hit never answers for another pair
  let bad = 0;
  for (let f = 0; f < fills; f++) {
    const fill = ((Math.imul(f + 1, 0x9e3779b1) & 0xffffff00) | 0xff) >>> 0;
    for (let k = 0; k < 3000; k++) {
      const ink = ((Math.imul(k + 7, 0x85ebca6b) & 0xffffff00) | 0xff) >>> 0;
      if (inkFor(fill, 0x1a1a1aff, ink) !== direct(fill, 0x1a1a1aff, ink)) bad++;
      if (inkFor(fill, ink, 0xf5f5f5ff) !== direct(fill, ink, 0xf5f5f5ff)) bad++;
    }
  }
  return bad;
}

describe('inkFor', () => {
  // the quick representatives: smaller sweeps; the full ones run in the medium lane
  it('equals the contrast rule, first and memoized, for any pair of inks', () => {
    expect(sweepPairs(1500)).toBe(0);
  });

  it('tells inks apart that share a memo slot: one ink shared, the other varied', () => {
    expect(sweepSlots(6)).toBe(0);
  });

  it('ignores the fill alpha', () => {
    expect(inkFor(0x336699ff, Black, White)).toBe(inkFor(0x33669980, Black, White));
  });
});

describe('inkFor, the full sweeps', { tags: ['medium'] }, () => {
  it('equals the contrast rule, first and memoized, for any pair of inks', () => {
    expect(sweepPairs(20000)).toBe(0);
  });

  it('tells inks apart that share a memo slot: one ink shared, the other varied', () => {
    expect(sweepSlots(60)).toBe(0);
  });
});
