/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The labels a section's layers offer at points (_spec/README.md "Render", Section layouts: Labels at points are the
// section's): settled together, the strongest kept, whichever layer offered it first.

import { describe, expect, it } from 'vitest';
import { TakenLabels, type TakenOwner } from '../sizing.js';

function owner(): TakenOwner & { kept: Map<number, boolean> } {
  const kept = new Map<number, boolean>();
  return { kept, placed: (id, k) => void kept.set(id, k) };
}

describe('TakenLabels', () => {
  it('the stronger of two labels that meet is kept, whichever was offered first', () => {
    const t = new TakenLabels();
    const first = owner();
    const second = owner();
    t.offer(first, 0, 0, 0, 40, 12, 20, 6, 24, 10);
    t.offer(second, 0, 10, 0, 50, 12, 30, 6, 299, 10);
    t.settle();
    expect(first.kept.get(0)).toBe(false);
    expect(second.kept.get(0)).toBe(true);
  });

  it('labels apart are all kept; tips within the clearance keep the stronger only', () => {
    const t = new TakenLabels();
    const o = owner();
    t.offer(o, 0, 0, 0, 20, 12, 10, 6, 1, 10);
    t.offer(o, 1, 100, 0, 120, 12, 110, 6, 2, 10);
    t.offer(o, 2, 100, 30, 120, 42, 105, 8, 3, 30);
    t.settle();
    expect([o.kept.get(0), o.kept.get(1), o.kept.get(2)]).toEqual([true, false, true]);
  });

  it('no label is kept over a wedge, the strongest aside; the lists are cleared by settling', () => {
    const t = new TakenLabels();
    const o = owner();
    t.block(0, 0, 100, 0, 50, 50);
    t.offer(o, 0, 40, 10, 60, 22, 50, 16, 9, 10);
    t.offer(o, 1, 200, 10, 220, 22, 210, 16, 5, 10);
    t.offer(o, 2, 45, 12, 65, 24, 55, 18, 1, 10);
    t.settle();
    // the strongest is kept whatever is under it; the one over it and over the wedge too is dropped
    expect([o.kept.get(0), o.kept.get(1), o.kept.get(2)]).toEqual([true, true, false]);
    const again = owner();
    t.offer(again, 0, 40, 10, 60, 22, 50, 16, 1, 10);
    t.settle();
    expect(again.kept.get(0)).toBe(true);
  });
});
