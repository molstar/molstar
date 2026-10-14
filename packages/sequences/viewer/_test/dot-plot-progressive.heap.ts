/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's progressive draw allocates nothing per frame (the memory rule; seqspec:viewer#performance): stacked
// lanes on a progressive backend, each lane's own grid, its fine tiles never built (the budget spent: none a frame
// once the coarse levels are), so every frame takes the budget, finds tiles missing, draws the coarse level instead,
// asks for another frame (`more`, the composer's `again`) and draws the plot's entries again alone, over two pans. The
// positive control: the same frames allocating an object each (`leaking`), through the same warm-up and measurement,
// are caught. Its own file: cases in one file share V8's feedback.

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

const N = 1_000;
/** The guards' measurement: warmed by count and by time, then `N` frames and `4N` measured. */
const measure = (frame: (i: number) => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (packages/sequences/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? '' : describeScaling(r, N);
};

/** Four stacked chains of 600, a grid on each, its coarse levels built and nothing finer. */
async function stacked() {
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { lit } = await import('./overlay-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { BackendCaps, Dirty } = await import('../src/draw/index.js');
  const { DotPlotLimits, DotPlotVisits } = await import('../src/draw/dot-plot.js');
  const L = 600;
  const b = createSeqViewSpec({ title: 'stacked' });
  const cs = Array.from({ length: 4 }, (_, k) => b.data.sequence({ ref: `c${k}`, value: 'M'.repeat(L) }));
  const value = cs.map(() => new Float32Array(L * L).map((_, x) => (x % 31) / 3));
  const g = b.data.annotation({
    ref: 'g',
    label: 'G',
    geometry: 'grid',
    value: { target: cs.map((_, k) => `c${k}`), b_target: cs.map((_, k) => `c${k}`), value },
  });
  const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', arrangement: 'stacked' });
  for (const c of cs) sec.block({ sequence: c });
  sec.track({ ref: 't', label: 'P' }).dot_plot({ ref: 'l-plot', data: g });
  const t = await lit(b.getDocument() as never, { width: 1200, height: 4000, dpr: 2 });
  const be = recording(BackendCaps.Progressive);
  // The first frame builds the coarse levels (forced); from then on no budget: the fine tiles never come.
  t.frame(be, Dirty.all);
  DotPlotLimits.samplesPerFrame = 0;
  const cell = t.vp?.live.current.cell ?? 1;
  let i = 0;
  const draw = () => {
    t.vp?.setView(cell * 2, i++ & 1 ? 40 : 0, 'api');
    t.frame(be, Dirty.all);
  };
  return { t, be, draw, visits: DotPlotVisits, dispose: () => t.ctx.dispose() };
}

describe("the dot plot's progressive draw allocates nothing per frame", () => {
  it('stacked lanes, tiles missing past the budget: the coarse level drawn, another frame asked, the entry alone', async () => {
    vi.resetModules();
    const g = await stacked();
    g.draw();
    const tiles = g.visits.tiles;
    g.draw();
    // Not vacuous: tiles drawn (the coarse ones), more asked, nothing built.
    const { Op } = await import('../src/testing/index.js');
    let images = 0;
    for (let k = 0; k < g.be.opCount; k++) if (g.be.ops[4 * k] === Op.image) images++;
    expect(images).toBeGreaterThan(0);
    expect(g.t.composer.again).toBe(true);
    expect(guard(g.draw)).toBe('');
    expect(g.visits.tiles).toBe(tiles);
    g.dispose();
  }, 180_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    vi.resetModules();
    const g = await stacked();
    expect(controlCaught(() => measure(leaking(g.draw)))).toBe(true);
    g.dispose();
  }, 180_000);
});
