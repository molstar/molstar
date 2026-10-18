/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's progressive case, for its heap guard and its quick twin (dot-plot-progressive.heap.ts,
// dot-plot-progressive.quick.ts): four stacked chains of 600 on a progressive backend, a grid on each, its coarse
// levels built by the first frame and nothing finer (`DotPlotLimits.samplesPerFrame` 0 from then on: the fine tiles
// never come); `draw` sets one of two views and draws a whole frame; `dispose` puts the budget back. Each call loads
// the viewer afresh (`vi.resetModules`), so the limit is this module copy's.

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/** The quick counts: a frame of the four lanes is about 0.5 ms (M4), so windows of 250 steps. */
export const DRAW: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

export async function stacked() {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { lit } = await import('./overlay-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { BackendCaps, Dirty } = await import('../src/draw/index.js');
  const { DotPlotLimits, DotPlotVisits } = await import('../src/draw/dot-plot.js');
  const { Op } = await import('../src/testing/index.js');
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
  const budget = DotPlotLimits.samplesPerFrame;
  DotPlotLimits.samplesPerFrame = 0;
  const cell = t.vp?.live.current.cell ?? 1;
  let i = 0;
  const draw = () => {
    t.vp?.setView(cell * 2, i++ & 1 ? 40 : 0, 'api');
    t.frame(be, Dirty.all);
  };
  const images = () => {
    let n = 0;
    for (let k = 0; k < be.opCount; k++) if (be.ops[4 * k] === Op.image) n++;
    return n;
  };
  return {
    t,
    be,
    draw,
    images,
    visits: DotPlotVisits,
    dispose: () => {
      DotPlotLimits.samplesPerFrame = budget;
      t.ctx.dispose();
    },
  };
}
