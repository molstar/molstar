/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' levels across frames, through the rows layer's frames (performance.md, Levels): a moving frame at a new
// level pair builds a budget of it and draws the level drawn last; a settled frame builds the rest and draws it.

import { describe, expect, it } from 'vitest';
import type { Frame, FrameInput, RowsLay } from '../index.js';
import { pyramidLevelFor } from '../pyramid.js';
import { settle } from './helpers.js';
import { Columns, rowsSpec as spec, Section, Synthetic, Track } from './rows-fixtures.js';

const Rows = 1_000;
const Size = { width: 1440, height: 800 };
const at = (rowScale: number, end: number, settled: boolean): FrameInput => ({
  size: Size,
  dpr: 2,
  settled,
  viewports: { [Section]: { start: 0, end, bands: { [Track]: { bandTop: 0, rowScale } } } },
});

function rowsLay(f: Frame): RowsLay {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
}

/** The level pair a frame at `rowScale` over `end` columns wants (as the layer picks it). */
function wanted(lay: RowsLay, rowScale: number) {
  return [pyramidLevelFor(1 / (rowScale * 2), 4), pyramidLevelFor(1 / lay.axis.devicePx, 4)] as const;
}

describe('the rows’ levels across frames', () => {
  it('a moving frame at a new level pair draws the level drawn last while the new one builds; settled, the new one', async () => {
    const s = await settle(spec(Rows, true), [Synthetic], { rows: { buildBudget: 20_000 } } as never);
    // the fit, settled: its level built and drawn
    const fit = rowsLay(await s.frame(at(0.4, Columns, true)));
    expect(fit.images.length).toBeGreaterThan(0);
    // a zoom out of the rows, moving: a level pair not built yet, a frame builds a part of it
    const f = s.renderer.frame(at(0.05, Columns, false));
    const lay = rowsLay(f);
    const [r, c] = wanted(lay, 0.05);
    expect(lay.pyramid!.isBuilt(r, c)).toBe(false);
    expect(lay.images.length).toBeGreaterThan(0);
    // settled: built and drawn
    const g = rowsLay(s.renderer.frame(at(0.05, Columns, true)));
    expect(g.pyramid!.isBuilt(r, c)).toBe(true);
  });

  it('the budget is the config’s `rows.buildBudget`: a level a default frame would build whole waits for a smaller one', async () => {
    // 500 × 300 base cells: within the default budget, past 20,000
    const frame = async (budget?: number) => {
      const s = await settle(spec(500, true), [Synthetic], (budget ? { rows: { buildBudget: budget } } : {}) as never);
      await s.frame(at(0.4, Columns, true));
      const lay = rowsLay(s.renderer.frame(at(0.05, Columns, false)));
      const [r, c] = wanted(lay, 0.05);
      return lay.pyramid!.isBuilt(r, c);
    };
    expect(await frame()).toBe(true);
    expect(await frame(20_000)).toBe(false);
  });

  it('a zoom out of the letters (no level drawn yet) draws the rows pending while the level builds, never a whole build', async () => {
    const s = await settle(spec(Rows, true), [Synthetic], { rows: { buildBudget: 20_000 } } as never);
    // the letters, settled: no image level drawn
    expect(rowsLay(await s.frame(at(20, 40, true))).images.length).toBe(0);
    // out to the fit, moving: a part built, nothing drawn of it
    const f = s.renderer.frame(at(0.4, Columns, false));
    let lay = rowsLay(f);
    const [r, c] = wanted(lay, 0.4);
    expect(lay.pyramid!.isBuilt(r, c)).toBe(false);
    expect(lay.images.length).toBe(0);
    // the rows in view drawn pending meanwhile, never an empty band (EC-169)
    expect(lay.pending.length).toBe(2);
    expect(lay.pending[1]).toBeGreaterThan(0);
    const rects = f.lists.passes.data.rects;
    const ink = lay.tokens.colors['seq-pending'];
    expect(Array.from({ length: rects.n }, (_, k) => rects.fill[k]).filter((x) => x === ink).length).toBe(1);
    // the frames after it build the rest, then draw it
    let frames = 0;
    for (let w = Columns - 1; !lay.images.length && frames < 100; w--, frames++)
      lay = rowsLay(s.renderer.frame(at(0.4, w, false)));
    expect(frames).toBeGreaterThan(1);
    expect(lay.images.length).toBeGreaterThan(0);
  });
});
