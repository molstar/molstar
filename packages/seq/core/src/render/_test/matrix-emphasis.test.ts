/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { location, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { type DrawLists, type Emphasis, type FrameInput, MatrixLay, recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { matrixSpec, Positions, States } from './matrix-fixtures.js';

// The matrix layer's emphasis (_spec/README.md "Render", The matrix layer, Emphasis): a member of its states system
// selected or hovered lights its row; contiguous rows are one band; a rect selected on its own cells lights no row.

const Input: FrameInput = { size: { width: 800, height: 320 }, dpr: 1, settled: true };

/** The pass's rects in the matrix's box, `[y, height, hasFill]`. */
function rows(lists: DrawLists, pass: 'selection' | 'hover', box: { x: number; width: number }) {
  const out: [number, number, boolean][] = [];
  for (const l of recordDrawLists(lists, [pass]).split('\n')) {
    if (!l.startsWith('rect ')) continue;
    const t = l.split(' ');
    const [x, y, w, h] = t.slice(1, 5).map(Number);
    if (x === box.x && w === box.width) out.push([y, h, t[5] !== 'none']);
  }
  return out;
}

async function open(emphasis: Emphasis) {
  const s = await settle(matrixSpec());
  const f = await s.frame({ ...Input, emphasis });
  const sec = f.layout.sections[0];
  const gt = sec.lays.flat().find((l) => l?.lay instanceof MatrixLay)!.lay as MatrixLay;
  return { f, gt };
}

describe('the matrix layer’s emphasis', () => {
  it('a member selected lights its row, contiguous members one band', async () => {
    const { f, gt } = await open({ selection: rangesLocation(States, [2, 3, 5, 8]) });
    const y = (s: number) => gt.box.y + s * gt.rowScale - gt.bandTop;
    expect(rows(f.lists, 'selection', gt.box).map(([a, h]) => [a, h])).toEqual([
      [y(2), gt.rowScale],
      [y(5), 3 * gt.rowScale],
    ]);
  });

  it('adjacent ranges are one band: no outline between them', async () => {
    const { f, gt } = await open({ selection: rangesLocation(States, [1, 2, 2, 3]) });
    expect(rows(f.lists, 'selection', gt.box)).toHaveLength(1);
  });

  it('a rect selected on the matrix’s own cells lights no row through the projected sides', async () => {
    const rect = location({
      rects: [{ a: { system: Positions, ranges: [3, 6] }, b: { system: States, ranges: [2, 3] } }],
    });
    const { f, gt } = await open({
      selection: rect,
      projected: { selection: rangesLocation(States, [2, 3]) },
    });
    expect(rows(f.lists, 'selection', gt.box)).toEqual([]);
  });

  it('a selection projected here lights its row lighter', async () => {
    const { f, gt } = await open({ projected: { selection: rangesLocation(States, [2, 3]) } });
    expect(rows(f.lists, 'selection', gt.box)).toHaveLength(1);
  });

  it('the hover’s tint is the config’s', async () => {
    const s = await settle(matrixSpec(), [], { emphasis: { hoverTint: 0.4 } });
    const base = await s.frame({ ...Input });
    const lay = base.layout.sections[0].lays.flat().find((l) => l?.lay instanceof MatrixLay)!.lay as MatrixLay;
    const hit = { layer: 'other', section: 'x', system: States, position: 3, target: 'row' } as const;
    const f = await s.frame({ ...Input, emphasis: { hover: hit } });
    const text = recordDrawLists(f.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect ') && l.split(' ')[1] === String(lay.box.x));
    expect(text).toHaveLength(1);
    // the fill's alpha: the hover ink's alpha at the config's tint
    const fill = text[0].split(' ')[5];
    expect(fill).not.toBe('none');
    const alpha = parseInt(fill.slice(-2), 16);
    expect(alpha).toBe(Math.round(0.4 * 255));
  });
});
