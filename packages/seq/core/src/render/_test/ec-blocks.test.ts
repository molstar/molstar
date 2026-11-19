/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Edge cases built but not tested before (.seq/sources/seqstar/edge-cases.md): EC-115 a stacked blocks section zoomed
// about a point keeps the block under the pointer; EC-99 a text section's ligand chips flow in a line with the words
// and wrap with them.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type BlockPlacement, PlacePoint, TextLayout, type TextPlacement, TracksLayout } from '../index.js';
import { seeded } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 600 };

function blocksSpec(o: { arrange?: 'row' | 'stack'; layout?: 'text'; ligands?: number; chains?: number; length?: number }) {
  const b = createSeqSpec({ title: 'blocks' });
  const r = seeded(11);
  const protein = (n: number) => {
    let s = '';
    for (let i = 0; i < n; i++) s += 'ACDEFGHIKLMNPQRSTVWY'[Math.floor(r() * 20)];
    return s;
  };
  const chains = Array.from({ length: o.chains ?? 2 }, (_, k) =>
    b.data.sequence(`c${k}`, { alphabet: 'protein', letters: protein(o.length ?? 90) }),
  );
  const ligands = Array.from({ length: o.ligands ?? 0 }, (_, k) =>
    b.data.sequence(`lig${k}`, { alphabet: 'protein', letters: 'X', name: `LIGAND${k}` }),
  );
  b.views.view('v', (v) => {
    const s = v.section('complex', {
      ...(o.arrange ? { arrange: o.arrange } : {}),
      ...(o.layout ? { layout: o.layout } : {}),
    });
    chains.forEach((c, k) => s.block(`b${k}`, { of: c, label: `chain ${k}` }));
    ligands.forEach((l, k) => s.block(`l${k}`, { of: l, label: `LIGAND${k}` }));
    s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

const placementOf = (f: { layout: { sections: readonly { placement?: unknown }[] } }) =>
  f.layout.sections[0].placement as BlockPlacement;

describe('EC-115: a stacked blocks section zoomed about a point', () => {
  it('the block under the pointer stays under it through successive zooms in', async () => {
    const s = await settle(blocksSpec({ arrange: 'stack', chains: 3 }));
    const f0 = await s.frame({ size: Size });
    const p0 = placementOf(f0);
    expect(p0.strips.axes.length).toBe(3);
    const axis = f0.layout.sections[0].axis!;
    // the pointer on the second block, a third of its way in
    const k = 1;
    const px = p0.strips.axes[k].x0 + p0.strips.axes[k].width / 3;
    const y = p0.strips.bandY[k] + 2;
    const at = p0.offsetAt(px, y);
    // a stack is one scale over the longest block's positions: the window is in a block's own units
    const total = axis.total;
    expect(total).toBe(90);
    let v = { start: 0, end: total };
    let p = p0;
    for (let i = 0; i < 3; i++) {
      const c = { total, width: axis.width, height: 600, cyclic: false, maxScale: 40, placement: p };
      v = TracksLayout.nav.zoom(c, v, 1.6, at) as typeof v;
      expect(v.end - v.start).toBeLessThan(total);
      const f = await s.frame({ size: Size, viewport: v });
      p = placementOf(f);
      // the same unit of the same block is at the same px, and the pointer there still reads it, in that block
      expect(p.strips.axes[k].xOfUnit(at)).toBeCloseTo(px, 1);
      const pt = new PlacePoint();
      expect(p.at([px, y], pt)).toBe(true);
      expect(pt.strip).toBe(k);
      expect(Math.floor(pt.offset)).toBe(Math.floor(at));
    }
  });
});

describe('EC-99: ligand chips in a text section', () => {
  const chipsOf = (p: TextPlacement) => p.strips.axes.map((_, i) => i).filter((i) => p.strips.chip[i] === 1);

  it('the chips share a line like words, and wrap onto the next line with them when the width runs out', async () => {
    const wide = await (await settle(blocksSpec({ layout: 'text', ligands: 6, chains: 1, length: 12 }))).frame({ size: Size });
    const pw = wide.layout.sections[0].placement as TextPlacement;
    expect(wide.layout.sections[0].layout).toBe(TextLayout);
    const cw = chipsOf(pw);
    expect(cw.length).toBe(6);
    // wide: one line, left to right, in order
    expect(new Set(cw.map((i) => pw.strips.top[i])).size).toBe(1);
    for (let j = 1; j < cw.length; j++)
      expect(pw.strips.axes[cw[j]].x0).toBeGreaterThan(pw.strips.axes[cw[j - 1]].x0 + pw.strips.axes[cw[j - 1]].width);

    const narrow = await (await settle(blocksSpec({ layout: 'text', ligands: 6, chains: 1, length: 12 }))).frame({
      size: { width: 260, height: 900 },
    });
    const pn = narrow.layout.sections[0].placement as TextPlacement;
    const cn = chipsOf(pn);
    expect(cn.length).toBe(6);
    const lines = [...new Set(cn.map((i) => pn.strips.top[i]))].sort((a, b) => a - b);
    expect(lines.length).toBeGreaterThan(1);
    // wrapped in order: a later chip is never on an earlier line, and none runs past the section's right edge
    for (let j = 1; j < cn.length; j++) expect(pn.strips.top[cn[j]]).toBeGreaterThanOrEqual(pn.strips.top[cn[j - 1]]);
    const right = pn.box.x + pn.box.width;
    for (const i of cn) expect(pn.strips.axes[i].x0 + pn.strips.axes[i].width).toBeLessThanOrEqual(right + 0.5);
    // the first line of chips starts after the words before it (the chains' lines), not at the section's top
    expect(pn.strips.top[cn[0]]).toBeGreaterThan(pn.strips.top[0]);
  });
});
