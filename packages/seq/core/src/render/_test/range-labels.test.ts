/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cut labels' cache (range-labels.ts): a cut is made once while its label is in use, however many labels are cut,
// so a warm layout makes no string.

import { describe, expect, it, vi } from 'vitest';
import { cutLabel } from '../range-labels.js';

describe('cut labels', () => {
  it('more labels cut than the cache’s first limit are made once while in use: a second pass makes none', () => {
    const labels = Array.from({ length: 6000 }, (_, i) => `a domain label ${i}`);
    for (const l of labels) expect(cutLabel(l, 5)).toBe(`${l.slice(0, 5)}…`);
    const slice = vi.spyOn(String.prototype, 'slice');
    for (const l of labels) cutLabel(l, 5);
    expect(slice).not.toHaveBeenCalled();
    slice.mockRestore();
  });

  it('labels no longer used are let go when new ones come', () => {
    const old = Array.from({ length: 5000 }, (_, i) => `an old label ${i}`);
    for (const l of old) cutLabel(l, 4);
    const fresh = Array.from({ length: 20_000 }, (_, i) => `a new label ${i}`);
    for (const l of fresh) cutLabel(l, 4);
    // the old ones were let go: asked again, they are made again
    const slice = vi.spyOn(String.prototype, 'slice');
    cutLabel(old[0], 4);
    expect(slice).toHaveBeenCalled();
    slice.mockRestore();
  });
});
