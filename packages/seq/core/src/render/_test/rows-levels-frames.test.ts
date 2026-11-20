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
});
