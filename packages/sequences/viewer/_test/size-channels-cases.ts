/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of size-channels.heap.ts and its quick twin (size-channels.quick.ts): frames of layers drawn through size
// themes (pairs stroked at a width theme's, sorted by color then width, with an end radius; ranges at a thickness;
// points at a radius; values at a thickness), through pans: a frame is a pan and a whole redraw. A fresh viewer
// (`vi.resetModules`).

import { vi } from 'vitest';

async function setup() {
  vi.resetModules();
  const { composed, recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(100) });
  const n = 300;
  const scale = (field: string, range: [number, number]) => ({ theme: 'field_scale', params: { field, range } });
  const start = Array.from({ length: n }, (_, i) => i * 6);
  const pairs = b.data.annotation({
    ref: 'pp',
    geometry: 'pairs',
    target: a,
    value: {
      start,
      b_target: start.map(() => 'a'),
      b_start: start.map((s) => s + 3 + (s % 5)),
      category: start.map((_, i) => `k${i % 4}`),
      count: start.map((_, i) => 1 + (i % 7)),
    },
  } as never);
  const items = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: {
      start,
      stop: start.map((s) => s + 4),
      conf: start.map((_, i) => 0.2 + (i % 5) / 5),
      category: start.map((_, i) => `c${i % 3}`),
    },
  } as never);
  const d = b.data.annotation({
    ref: 'd',
    geometry: 'values',
    target: a,
    value: { index: Array.from({ length: 1500 }, (_, i) => i), value: Array.from({ length: 1500 }, (_, i) => i % 11) },
  });
  const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
  v.track({ label: 'Pairs' }).pairs({
    data: pairs,
    width: scale('count', [1, 6]),
    radius: scale('count', [2, 5]),
  } as never);
  v.track({ label: 'Bars' }).ranges({ data: items, thickness: scale('conf', [0.2, 1]) } as never);
  v.track({ label: 'Dots' }).points({ data: items, radius: scale('conf', [2, 5]) } as never);
  v.track({ label: 'Values' }).values({ data: d, thickness: scale('value', [0.3, 1]) } as never);
  const t = await composed(b.getDocument(), { width: 1000, height: 400, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(4, 0, 'api');
  t.frame(be, 7);
  return { t, be };
}

/** The case: a frame (it gets a running index), a pan and a whole redraw. */
export async function sizeFrame() {
  const { t, be } = await setup();
  return (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, 'pointer');
    t.frame(be, 7);
  };
}
