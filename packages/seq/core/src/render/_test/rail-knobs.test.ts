/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rail's knobs (_spec/README.md "Render", Navigators, The ruler rail): discs centred on the window's true ends at
// every size (the axis's ends and a wrapped window's second end too), a ring in the well's ground under each; hovered
// they grow about their centre by at least 4 px, dragged more.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type DrawLists, type FrameInput, hitView } from '../index.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 300 };
const Section = '/views/v/protein';
const rails = { [Section]: true };
const active = { [Section]: true };

function p53(cyclic = false) {
  const b = createSeqSpec({ title: 'knobs' });
  const p = b.data.sequence('p53', { alphabet: 'protein', letters: P53, ...(cyclic ? { cyclic: true } : {}) });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The window pass's discs, in draw order (a ring then a disc at each end): centre along x and y, and diameter. */
function discs(lists: DrawLists) {
  const l = lists.passes.window.paths;
  const out: { cx: number; cy: number; d: number }[] = [];
  for (let i = 0; i < l.n; i++) {
    let lo = Infinity;
    let hi = -Infinity;
    let ylo = Infinity;
    let yhi = -Infinity;
    const n = l.cmdCount[i] - 1;
    for (let k = 0; k < n; k++) {
      const x = l.points[2 * (l.pointStart[i] + k)];
      const y = l.points[2 * (l.pointStart[i] + k) + 1];
      lo = Math.min(lo, x);
      hi = Math.max(hi, x);
      ylo = Math.min(ylo, y);
      yhi = Math.max(yhi, y);
    }
    out.push({ cx: (lo + hi) / 2, cy: (ylo + yhi) / 2, d: hi - lo });
  }
  return out;
}

async function frame(viewport: { start: number; end: number }, cyclic = false) {
  const s = await settle(p53(cyclic));
  const input: FrameInput = { size: Size, viewport, rails, railsActive: active };
  const f = await s.frame(input);
  return { s, f, input, rail: f.layout.sections[0].rail! };
}

describe('the rail’s knobs', () => {
  it('are circles centred on the window’s ends, rings under them, at rest', async () => {
    const { f, rail } = await frame({ start: 100, end: 200 });
    const k = rail.width / rail.total;
    const d = discs(f.lists);
    expect(d).toHaveLength(4);
    // ring, disc, ring, disc: a disc a circle (its extent in y the same as in x), the ring 2 px wider
    for (const x of d) expect(x.d).toBeGreaterThan(0);
    expect(d[1].cx).toBeCloseTo(rail.x + 100 * k, 3);
    expect(d[3].cx).toBeCloseTo(rail.x + 200 * k, 3);
    expect(d[0].cx).toBeCloseTo(d[1].cx, 6);
    expect(d[0].d).toBeCloseTo(d[1].d + 2, 3);
    for (const x of d) expect(x.cy).toBeCloseTo(rail.y, 3);
    // round: the polygon's extent in y is its extent in x
    const l = f.lists.passes.window.paths;
    let yl = Infinity;
    let yh = -Infinity;
    for (let k2 = 0; k2 < l.cmdCount[1] - 1; k2++) {
      const y = l.points[2 * (l.pointStart[1] + k2) + 1];
      yl = Math.min(yl, y);
      yh = Math.max(yh, y);
    }
    expect(yh - yl).toBeCloseTo(d[1].d, 2);
  });

  it('at both ends of the axis stay centred on the window’s end, hovered too (no clamp)', async () => {
    for (const [viewport, at] of [
      [{ start: 0, end: 100 }, 'start'],
      [{ start: 293, end: 393 }, 'end'],
    ] as const) {
      const { s, f, input, rail } = await frame(viewport);
      const edge = at === 'start' ? rail.x : rail.x + rail.width;
      const i = at === 'start' ? 1 : 3;
      expect(discs(f.lists)[i].cx).toBeCloseTo(edge, 3);
      const h = { ...hitView(f.layout, rail.x + rail.width / 2, rail.bottom - 2)! };
      const g = discs(s.renderer.frame({ ...input, emphasis: { hover: h } }).lists);
      expect(g[i].cx).toBeCloseTo(edge, 3);
    }
  });

  it('a hovered knob grows about its centre by at least 4 px, a dragged one more', async () => {
    const { s, f, input, rail } = await frame({ start: 100, end: 200 });
    const rest = discs(f.lists);
    const h = { ...hitView(f.layout, rail.x + (150 / 393) * rail.width, rail.bottom - 2)! };
    const hover = discs(s.renderer.frame({ ...input, emphasis: { hover: h } }).lists);
    const drag = discs(s.renderer.frame({ ...input, emphasis: { drag: h } }).lists);
    for (const i of [1, 3]) {
      expect(hover[i].d - rest[i].d).toBeGreaterThanOrEqual(4);
      expect(hover[i].cx).toBeCloseTo(rest[i].cx, 3);
      expect(hover[i].cy).toBeCloseTo(rest[i].cy, 3);
      expect(drag[i].d).toBeGreaterThan(hover[i].d);
    }
  });

  it('a window over a cyclic origin has its second knob at the true end, where the hit grabs it', async () => {
    const { f, rail } = await frame({ start: 380, end: 410 }, true);
    expect(rail.end).toBeGreaterThan(rail.total);
    const k = rail.width / rail.total;
    const d = discs(f.lists);
    expect(d[1].cx).toBeCloseTo(rail.x + 380 * k, 3);
    const end = rail.x + (rail.end - rail.total) * k;
    expect(d[3].cx).toBeCloseTo(end, 3);
    expect(hitView(f.layout, end + 1, rail.y)?.part).toBe('end');
  });
});
