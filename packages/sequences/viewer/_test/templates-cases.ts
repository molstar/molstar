/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The template's 1,000 instance tracks as a step, for the heap guard and its quick twin (templates.heap.ts,
// templates.quick.ts): the template document (template-fixtures.ts) with 1,000 instances (ranges on the protein),
// composed at 1,000 x 400, DPR 2; a step is a pan's pass of the streams and a composer frame. The viewer is loaded
// afresh per call (`vi.resetModules`).

import { vi } from 'vitest';
import type { Streams } from '../src/streams.js';

/** The instances the template document is given. */
export const INSTANCES = 1_000;

/** The template document with 1,000 instances (ranges on the protein), composed at 1,000 x 400, DPR 2. */
export async function templated() {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/composer.js');
  const { range, templateDoc, templateLoaders, P } = await import('./template-fixtures.js');
  const genes = Array.from({ length: INSTANCES }, (_, k) => range(`g${k}`, 1 + (k % (P.length - 8))));
  const { loaders } = templateLoaders({ genes, notes: [] });
  const t = await composed(templateDoc(), { width: 1000, height: 400, dpr: 2, loaders });
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  vp.zoomTo(0, 20, 120);
  const be = recording();
  for (let k = 0; k < 8; k++) {
    vp.panBy(k & 1 ? -9 : 9, 0, 'pointer');
    t.frame(be, Dirty.all);
  }
  let i = 0;
  return {
    t,
    frame: () => {
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      (t.ctx.streams as Streams).pass();
      t.frame(be, Dirty.all);
    },
  };
}
