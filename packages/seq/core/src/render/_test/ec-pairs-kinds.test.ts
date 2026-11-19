/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EC-63 (.seq/sources/seqstar/edge-cases.md): two pair kinds at one cell of a dot plot are sliced, side by side in the
// cell, each slice its own hit.

import { describe, expect, it } from 'vitest';
import { hitView, type GridLay, type PairsLay } from '../index.js';
import { settle } from './helpers.js';
import { gridSpec, pairRows } from './pairs-fixtures.js';

const Size = { width: 720, height: 900 };

describe('EC-63: two pair kinds at one cell', () => {
  it('the cell is cut into a slice a kind, in order, each slice hit as its own pair', async () => {
    // pairs (4, 12) of two kinds, and (6, 9) alone as the control
    const rows = {
      ...pairRows(
        [
          [4, 12],
          [4, 12],
          [6, 9],
        ],
        ['canonical', 'wobble', 'canonical'],
      ),
      id: ['wc', 'gu', 'other'],
    };
    const s = await settle(gridSpec({ n: 20, rows }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const sec = f.layout.sections[1];
    const lay = sec.lays[0][1].lay as unknown as PairsLay;
    const g = sec.lays[0][0].lay as GridLay;
    const a = g.axis!;
    const ay = g.axisY!;
    expect(lay.nd).toBe(3);
    a.cell(4);
    const cx0 = Math.min(a.sx0, a.sx1);
    const cw = Math.abs(a.sx1 - a.sx0);
    ay.cell(12);
    const cy = (ay.sx0 + ay.sx1) / 2;
    // the two dots at the cell: halves of its width, in the pairs' order, side by side; the control keeps its cell whole
    const dots = [0, 1, 2].map((i) => ({ x: lay.dots[4 * i], w: lay.dots[4 * i + 2] }));
    expect(dots[0].w).toBeCloseTo(cw / 2, 6);
    expect(dots[1].w).toBeCloseTo(cw / 2, 6);
    expect(dots[0].x).toBeCloseTo(cx0, 6);
    expect(dots[1].x).toBeCloseTo(cx0 + cw / 2, 6);
    expect(dots[2].w).toBeCloseTo(cw, 6);
    // each slice is its own hit
    expect(hitView(f.layout, cx0 + cw / 4, cy)?.element?.key).toBe('wc');
    expect(hitView(f.layout, cx0 + (3 * cw) / 4, cy)?.element?.key).toBe('gu');
    // and the two slices are two colors, one a kind
    expect(lay.dotColor[0]).not.toBe(lay.dotColor[1]);
    expect(lay.dotColor[2]).toBe(lay.dotColor[0]);
  });
});
