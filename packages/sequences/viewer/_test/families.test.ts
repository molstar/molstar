/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Vertical scrolling of a box section's rows (LY1; the viewport's clamp): limits after a resize and a row zoom, and the
// vertical overview's window and the scrollbar's thumb from the same numbers. On a synthetic family.

import { describe, expect, it } from 'vitest';
import { FAMILY_SECTION as F, syntheticFamily } from './families.js';
import { opened } from './layout-fixtures.js';

describe("a box section's rows scroll exactly over their content", () => {
  it('clamps to [0, rows − visible] after a resize and a row zoom; the windows match', async () => {
    const t = await opened(syntheticFamily(120), { size: { [F]: [1000, 300] } });
    const vp = t.vp(F);
    const s = vp.live.current;
    const l = () => t.layout(F);
    const rows = () => l().rows.count * l().rows.height;
    const visible = () => s.height - l().scrolledTop;
    expect([l().rows.count, s.contentHeight]).toEqual([120, rows()]);
    vp.scrollTo(0, 1e9, 'api');
    expect(s.scrollY).toBe(rows() - visible());
    vp.setSize(1000, 200, 1);
    l(); // the layout rebuilds on read
    vp.scrollTo(0, 1e9, 'api');
    expect(s.scrollY).toBe(rows() - visible());
    vp.setRowScale(0.5, 0);
    l();
    expect(s.contentHeight).toBe(rows());
    expect(s.scrollY).toBeLessThanOrEqual(rows() - visible());
    vp.scrollTo(0, 1e9, 'api');
    expect(s.scrollY).toBe(rows() - visible());
    // The overview's window and the thumb: the visible share of the rows, at the scroll's place.
    const v = t.ctx.layout.axis(F);
    if (!v) throw new Error('no axis view');
    const rh = l().height - l().scrolledTop;
    expect(v.vnavY1 - v.vnavY0).toBeCloseTo((rh * visible()) / rows(), 6);
    expect(v.vnavY1).toBeCloseTo(l().height, 6);
    vp.scrollTo(0, 0, 'api');
    const w = t.ctx.layout.axis(F);
    expect(w?.vnavY0).toBe(l().scrolledTop);
    expect(w?.vnavY1).toBeCloseTo(l().scrolledTop + (rh * visible()) / rows(), 6);
  });

  it('rows that fit do not scroll, and show no vertical navigation', async () => {
    const t = await opened(syntheticFamily(4, 4), { size: { [F]: [1000, 400] } });
    const vp = t.vp(F);
    vp.scrollTo(0, 50, 'api');
    expect([vp.live.current.scrollY, t.layout(F).vnavWidth]).toEqual([0, 0]);
  });
});
