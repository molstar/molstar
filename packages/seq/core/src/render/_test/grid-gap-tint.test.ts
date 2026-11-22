/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DarkTokens, GridPlacement, LightTokens, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { twoSectionsSpec } from './grid-two-sections-fixtures.js';
import { letters } from './pairs-fixtures.js';

// The block gaps of a grid over two sections of blocks (_spec/README.md "Render", The grid layout): tinted a step off
// the cells' ground (`msui-line`) under the data, each exactly the space between two blocks, never a line; the grid over
// one section of blocks has none.

const Size = { width: 720, height: 900 };

const hex = (c: number) => `#${c.toString(16).padStart(8, '0').slice(0, 6)}`;

/** The background pass's rects in `ink` below `top`, as `[x, y, w, h]`. */
function tints(lists: Parameters<typeof recordDrawLists>[0], ink: string, top: number) {
  return recordDrawLists(lists, ['background'])
    .split('\n')
    .filter((l) => l.startsWith('rect ') && l.includes(ink))
    .map((l) => l.split(' ').slice(1, 5).map(Number))
    .filter((r) => r[1] + r[3] > top);
}

function single(): Node {
  const b = createSeqSpec();
  const ns = [10, 20, 15];
  const seqs = ns.map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'rna', letters: letters(n, k + 1) }));
  const grids: string[] = [];
  ns.forEach((nx, x) =>
    ns.forEach((ny, y) =>
      grids.push(b.data.grid(`g${x}${y}`, { x: seqs[x], y: seqs[y], values: new Float32Array(nx * ny) })),
    ),
  );
  b.views.view('v', (v) => {
    const s = v.section('chains', {});
    seqs.forEach((q, k) => s.block(`b${k}`, { of: q, label: `chain ${k}` }));
    s.track('letters').layer('letters', { type: 'letters' });
    v.section('plot', { layout: 'grid', follow: s.address }).track('plot', {}, (t) =>
      grids.forEach((d, k) => t.layer(`g${k}`, { type: 'grid', data: d })),
    );
  });
  return b.build();
}

describe('the gaps of a grid over two sections of blocks', () => {
  it('are tinted exactly between the blocks (x and y), in both themes, canvas and SVG alike', async () => {
    const spec = twoSectionsSpec({ xs: [10, 20, 15], ys: [12, 18] });
    const st = await settle(spec);
    for (const dpr of [1, 2]) {
      const f = await st.frame({ size: Size, settled: true, dpr });
      const plot = f.layout.sections[2].placement as GridPlacement;
      const xr = plot.xRow!;
      const yr = plot.yRow!;
      const ink = hex(LightTokens.colors['msui-line']);
      const rects = tints(f.lists, ink, plot.plotTop - 1);
      const down = rects.filter(([, , , h]) => h === plot.plotHeight);
      const across = rects.filter(([, , w]) => w === xr.width);
      // one gap between each pair of neighbouring blocks, a rect each, no line
      expect(
        down.map(([x, , w]) => [x, w]),
        `dpr ${dpr}`,
      ).toEqual([1, 2].map((b) => [xr.x[b - 1] + xr.w[b - 1], xr.x[b] - (xr.x[b - 1] + xr.w[b - 1])]));
      expect(
        across.map(([, y, , h]) => [y, h]),
        `dpr ${dpr}`,
      ).toEqual([[yr.x[0] + yr.w[0], yr.x[1] - (yr.x[0] + yr.w[0])]]);
      for (const [, , w, h] of rects) expect(Math.min(w, h)).toBeGreaterThan(1);
    }
    // both themes' exports draw the tint, in their own token
    for (const tokens of [LightTokens, DarkTokens]) {
      const r = await renderSvg(spec, { size: Size, tokens });
      expect(r.svg.toLowerCase()).toContain(hex(tokens.colors['msui-line']));
    }
  });

  it('the grid over one section of blocks has none', async () => {
    const st = await settle(single());
    const f = await st.frame({ size: Size, settled: true });
    const plot = f.layout.sections[1].placement as GridPlacement;
    expect(tints(f.lists, hex(LightTokens.colors['msui-line']), plot.plotTop - 1)).toEqual([]);
  });
});
