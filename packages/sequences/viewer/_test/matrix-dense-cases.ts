/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS2's dense matrix cases, for their heap guard and its quick twin (matrix-dense.heap.ts, matrix-dense.quick.ts): a
// 3,000-state genotype matrix drawn through its per-device-row tiles at a device pixel ratio, or two stacked blocks of
// one state, each step a pointer move and a whole frame, a fresh viewer for each (`vi.resetModules`).

import { vi } from 'vitest';

/** The cases: a device pixel ratio, and two blocks stacked (height·dpr not whole at 4/3). */
export const DENSE = [
  [1.5, false],
  [1, false],
  [4 / 3, true],
] as const;

/** The genotype matrix composed, and a step: a pointer move and a whole frame, the zoom switched every 64 steps. */
export async function setup(dpr: number, stacked = false) {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { genotypeDoc } = await import('./genotype-fixtures.js');
  const { entryOf } = await import('./draw-fixtures.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const { MatrixVisits } = await import('../src/draw/matrix-tiles.js');
  const doc = genotypeDoc({ samples: 3000, positions: 2000, stacked });
  const t = await composed(doc, { width: 1000, height: 400, dpr });
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  const be = recording();
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const cells = [vp.live.current.cell * 40, vp.live.current.cell];
  const e = entryOf(t.ctx, 'G');
  const step = (i: number) => {
    // Every 64 steps the other zoom; the pointer down and across the matrix's rows between.
    if ((i & 63) === 0) vp.setView(cells[(i >> 6) & 1], 0, 'api');
    p.x = e.left + 3.5 + ((i * 7) % 600);
    p.y = e.top + 0.25 + ((i * 0.37) % (e.height - 0.5));
    t.ctx.input.pointer(p);
    t.composer.redraw();
    t.frame(be, 7);
  };
  return { step, visits: MatrixVisits, hover: () => t.ctx.hover.live.current.state };
}
