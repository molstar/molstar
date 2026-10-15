/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A split track's row guides' case, for its heap guard and its quick twin (split-guides.heap.ts,
// split-guides.quick.ts): case 13's variants grid (22 rows) at DPR 2, a frame a pan and a whole redraw, counting the
// frames that drew the 21 guides (the grid tone at the guide's alpha, 2 device px tall). A fresh viewer
// (`vi.resetModules`).

import { vi } from 'vitest';

/** Case 13 composed at DPR 2, and a frame: a pan and a whole frame, the guides counted. */
export async function guidesFrame() {
  vi.resetModules();
  const { Dirty, SkinToken } = await import('../src/draw/index.js');
  const { COPY_BACK, composed, recording } = await import('./composer-fixtures.js');
  const { case13 } = await import('../../core/_test/format/cases.js');
  const t = await composed(case13(), { width: 1000, height: 900, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(3, 200, 'api');
  t.frame(be, Dirty.all);
  const GRID = SkinToken.grid;
  let guides = 0;
  let full = 0;
  const frame = (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
    t.frame(be, Dirty.all);
    const r = be.rectItems;
    let n = 0;
    for (let j = 0; j < r.n; j++) if (r.color[j] === GRID && r.alpha[j] < 255 && r.h[j] === 2) n++;
    guides += n;
    if (n >= 21) full++;
  };
  return { frame, full: () => full, guides: () => guides };
}
